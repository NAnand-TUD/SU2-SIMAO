/*!
 * \file driver_structure.cpp
 * \brief The main subroutines for driving multi-zone problems.
 * \author R. Sanchez, O. Burghardt
 * \version 6.0.1 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/driver_structure.hpp"
#include "../include/definition_structure.hpp"


CMultizoneDriver::CMultizoneDriver(char* confFile,
                       unsigned short val_nZone,
                       unsigned short val_nDim,
                       bool val_periodic,
                       SU2_Comm MPICommunicator) : CDriver(confFile,
                                                          val_nZone,
                                                          val_nDim,
                                                          val_periodic,
                                                          MPICommunicator) {

  /*--- Initialize the counter for TimeIter ---*/
  TimeIter = 0;

  /*--- Initialize some useful booleans ---*/
  fsi = false; cht = false;

  /*--- Structure for multizone convergence ---*/
  init_res     = new su2double*[nZone];
  residual     = new su2double*[nZone];
  residual_rel = new su2double*[nZone];
  nVarZone     = new unsigned short[nZone];

  for (iZone = 0; iZone < nZone; iZone++){
    nVarZone[iZone] = 0;
    /*--- Account for all the solvers ---*/
    for (unsigned short iSol = 0; iSol < MAX_SOLS; iSol++){
      if (solver_container[iZone][INST_0][MESH_0][iSol] != NULL)
        nVarZone[iZone] += solver_container[iZone][INST_0][MESH_0][iSol]->GetnVar();
    }
    init_res[iZone]     = new su2double[nVarZone[iZone]];
    residual[iZone]     = new su2double[nVarZone[iZone]];
    residual_rel[iZone] = new su2double[nVarZone[iZone]];
    /*--- Initialize the residual containers to 0.0 ---*/
    for (unsigned short iVar = 0; iVar < nVarZone[iZone]; iVar++){
      init_res[iZone][iVar]     = 0.0;
      residual[iZone][iVar]     = 0.0;
      residual_rel[iZone][iVar] = 0.0;
    }
  }

  /*----------------------------------------------------*/
  /*------ Determine the properties of the problem -----*/
  /*----------------------------------------------------*/

  bool structural_zone = false;
  bool fluid_zone = false;
  bool heat_zone = false;

  /*--- If there is at least a fluid and a structural zone ---*/
  for (iZone = 0; iZone < nZone; iZone++){
    switch (config_container[iZone]->GetKind_Solver()) {
    case EULER: case NAVIER_STOKES: case RANS:
      fluid_zone = true;
      break;
    case FEM_ELASTICITY: case FEM_MODAL:
      structural_zone = true;
      break;
    case HEAT_EQUATION_FVM:
      heat_zone = true;
      break;
    }
  }

  /*--- If the problem has FSI properties ---*/
  if (fluid_zone && structural_zone) fsi = true;
  /*--- If the problem has CHT properties ---*/
  if (fluid_zone && heat_zone) cht = true;

  /*----------------------------------------------------*/
  /*- Define if a prefixed motion is imposed in a zone -*/
  /*----------------------------------------------------*/

  prefixed_motion = new bool[nZone];
  for (iZone = 0; iZone < nZone; iZone++){
    switch (config_container[iZone]->GetKind_GridMovement()){
      case RIGID_MOTION: case DEFORMING:
      case EXTERNAL: case EXTERNAL_ROTATION:
      case AEROELASTIC: case AEROELASTIC_RIGID_MOTION:
      case ELASTICITY:
        prefixed_motion[iZone] = true; break;
      case FLUID_STRUCTURE: case FLUID_STRUCTURE_STATIC:
      case STEADY_TRANSLATION: case MOVING_WALL: case ROTATING_FRAME:
      case NO_MOVEMENT: case GUST: default:
        prefixed_motion[iZone] = false; break;
    }
  }

  if(config_container[ZONE_0]->GetUnsteady_Simulation() == HARMONIC_BALANCE) {
      D = NULL;
      /*--- allocate dynamic memory for the Harmonic Balance operator ---*/
      D = new su2double *[nInst[ZONE_0]];
      for (iInst = 0; iInst < nInst[ZONE_0]; iInst++) D[iInst] = new su2double[nInst[ZONE_0]];
  }

}

CMultizoneDriver::~CMultizoneDriver(void) {

  for (iZone = 0; iZone < nZone; iZone++){
    delete [] init_res[iZone];
    delete [] residual[iZone];
    delete [] residual_rel[iZone];
  }

  delete [] init_res;
  delete [] residual;
  delete [] residual_rel;

  delete [] prefixed_motion;

}

void CMultizoneDriver::StartSolver() {

  /*--- Main external loop of the solver. Runs for the number of time steps required. ---*/

//   iteration_container[0][INST_0]->Solve(output, integration_container, geometry_container, solver_container,numerics_container, config_container, surface_movement, grid_movement, FFDBox, 0, INST_0);
  
  if (rank == MASTER_NODE)
    cout << endl <<"------------------------------ Begin Solver -----------------------------" << endl;

  if (rank == MASTER_NODE){
    cout << endl <<"Simulation Run using the Multizone Driver" << endl;
    if (driver_config->GetTime_Domain())
      cout << "The simulation will run for " << driver_config->GetnTime_Iter() << " time steps." << endl;
  }

  /*--- Set the initial time iteration to the restart iteration. ---*/
  if (driver_config->GetRestart()) TimeIter = driver_config->GetRestart_Iter();

  /*--- Run the problem until the number of time iterations required is reached. ---*/
  while ( TimeIter < driver_config->GetnTime_Iter() ) {

    /*--- Perform some preprocessing before starting the time-step simulation. ---*/
    cout << "\n\nMZ driver timeiter loop -- preprocess\n";
    Preprocess(TimeIter);

    /*--- Run a block iteration of the multizone problem. ---*/

    switch (driver_config->GetKind_MZSolver()){
      case MZ_BLOCK_GAUSS_SEIDEL: Run_GaussSeidel(); break;  // Block Gauss-Seidel iteration
      case MZ_BLOCK_JACOBI: Run_Jacobi(); break;             // Block-Jacobi iteration
      default: Run_GaussSeidel(); break;
    }

    /*--- Update the solution for dual time stepping strategy ---*/

    Update();

    /*--- Monitor the computations after each iteration. ---*/

    Monitor(TimeIter);

    /*--- Output the solution in files. ---*/
    /* COmmented Out temp*/
    Output(TimeIter);

    /*--- If the convergence criteria has been met, terminate the simulation. ---*/

    if (StopCalc) break;

    TimeIter++;

  }

}

void CMultizoneDriver::Preprocess(unsigned long TimeIter) {

  bool unsteady = driver_config->GetTime_Domain();

  cout << "\n\nmz solver; unsteady= " << unsteady << endl;
  for (iZone = 0; iZone < nZone; iZone++){

    /*--- Set the value of the external iteration to TimeIter. -------------------------------------*/
    /*--- TODO: This should be generalised for an homogeneous criteria throughout the code. --------*/
    config_container[iZone]->SetExtIter(TimeIter);

    /*--- Read the target pressure for inverse design. ---------------------------------------------*/
    /*--- TODO: This routine should be taken out of output, and made general for multiple zones. ---*/
    if (config_container[iZone]->GetInvDesign_Cp() == YES)
      output->SetCp_InverseDesign(solver_container[iZone][INST_0][MESH_0][FLOW_SOL],
          geometry_container[iZone][INST_0][MESH_0], config_container[iZone], TimeIter);

    /*--- Read the target heat flux ----------------------------------------------------------------*/
    /*--- TODO: This routine should be taken out of output, and made general for multiple zones. ---*/
    if (config_container[iZone]->GetInvDesign_HeatFlux() == YES)
      output->SetHeatFlux_InverseDesign(solver_container[iZone][INST_0][MESH_0][FLOW_SOL],
          geometry_container[iZone][INST_0][MESH_0], config_container[iZone], TimeIter);

    /*--- Set the initial condition for EULER/N-S/RANS ---------------------------------------------*/
    /*--- For FSI, this is set after the mesh has been moved. --------------------------------------*/
    if ((config_container[iZone]->GetKind_Solver() ==  EULER) ||
        (config_container[iZone]->GetKind_Solver() ==  NAVIER_STOKES) ||
        (config_container[iZone]->GetKind_Solver() ==  RANS) ) {
        if(!fsi) solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container[iZone][INST_0], solver_container[iZone][INST_0], config_container[iZone], TimeIter);
    }

  }

#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif

  /*--- Run a predictor step ---*/
  for (iZone = 0; iZone < nZone; iZone++){
    if (config_container[iZone]->GetPredictor())
      iteration_container[iZone][INST_0]->Predictor(output, integration_container, geometry_container, solver_container,numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone, INST_0);
  }

  /*--- Perform a dynamic mesh update if required. ---*/

  DynamicMeshUpdate(TimeIter);

  cout<<"After dynamic Mesh Update Line 249 driver_direct_multizone \n";

  cout << "unsteady flag1: " << unsteady << endl;
//   cout << "unsteady flag2: " << config_container[iZone]->GetUnsteady_Simulation() << endl;
  /*--- Updating zone interface communication patterns for unsteady problems with pre-fixed motion in the config file ---*/
  if ( unsteady /*|| config_container[iZone]->GetUnsteady_Simulation()==HARMONIC_BALANCE*/) {
    cout << "here mz driver 260\n";
    for (iZone = 0; iZone < nZone; iZone++) {
      for (unsigned short jZone = 0; jZone < nZone; jZone++){
          cout << "zones: " << iZone << "\t" << jZone << "\n";
          for(iInst=0; iInst<nInst[iZone]; iInst++)
                cout << "isnts: " << iInst << "\n";
                if(jZone != iZone && interpolator_container[iZone][jZone][iInst] != NULL && prefixed_motion[iZone]){
                    cout << "here mz driver 264\n";
                    interpolator_container[iZone][jZone][iInst]->Set_TransferCoeff(config_container);
                }
      }
    }
  }

}

void CMultizoneDriver::Run_GaussSeidel() {

  unsigned long iOuter_Iter;
  unsigned short jZone, UpdateMesh;
  bool DeformMesh = false;
  unsigned long ExtIter = 0;
  bool Convergence = false;

  cout << "\n\n run gauss-seidel\n\n";
  unsigned long OuterIter = 0; 
  for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetOuterIter(OuterIter);


  /*--- Loop over the number of outer iterations ---*/
  for (iOuter_Iter = 0; iOuter_Iter < driver_config->GetnOuter_Iter(); iOuter_Iter++){
    cout << "\n\n run GS main outer iteration " << iOuter_Iter <<"\n\n";
    /*--- Loop over the number of zones (IZONE) ---*/
    for (iZone = 0; iZone < nZone; iZone++){

      /*--- In principle, the mesh does not need to be updated ---*/
      UpdateMesh = 0;
    cout << "\n\n run GS: zone --> " << iZone <<"\n(Fluid: zone 0; Structure: zone 1)\n\n";

      /*--- Set the OuterIter ---*/
      config_container[iZone]->SetOuterIter(iOuter_Iter);
        cout << "completed SetOuterIter\n";
//         if(iOuter_Iter > 1) exit(0);
        
      /*--- Transfer from all the remaining zones ---*/
      for (jZone = 0; jZone < nZone; jZone++){
        /*--- The target zone is iZone ---*/
        if (jZone != iZone){
            cout << "\n\n run GS trasfer data:";
            if(jZone == 0) cout << " donor=fluid; target=structure\n\n";
            if(jZone == 1) cout << " donor=structure; target=fluid\n\n";
            DeformMesh = Transfer_Data(jZone, iZone);
            cout << "finish data transfer\n";
            if (DeformMesh) UpdateMesh+=1;
        }
      }
      /*--- If a mesh update is required due to the transfer of data ---*/
    cout << "\n\n run GS dynupdatemesh --> " << iZone <<" updatemesh= " << UpdateMesh <<"\n";
      if (UpdateMesh > 0) DynamicMeshUpdate(iZone, ExtIter);


      /*--- Iterate the zone as a block, either to convergence or to a max number of iterations ---*/
      for (iInst = 0; iInst < nInst[iZone]; iInst++) {
          cout << "\n\n run GS iteration->solve --> " << iZone << "\n(Fluid: zone 0; Structure: zone 1)\n";
          iteration_container[iZone][iInst]->Solve(output, integration_container, geometry_container, solver_container,
                                                   numerics_container, config_container, surface_movement,
                                                   grid_movement, FFDBox, iZone, iInst);
          cout << "\n\n completed iteration->solve\n";

          /*--- A corrector step can help preventing numerical instabilities ---*/
          if (config_container[0]->GetKind_Solver() != FEM_MODAL &&
              config_container[1]->GetKind_Solver() != FEM_MODAL)
              Corrector(iZone);
      }
    cout<<"iZone :: "<<iZone<<" Modal :: "<<(config_container[iZone][INST_0].GetDynamic_Method() == MODAL_HARMONIC_BALANCE)<<
            "  Fluid :: "<<(config_container[iZone][INST_0].GetUnsteady_Simulation() == HARMONIC_BALANCE)<<endl;
    // Solver Update
    if (config_container[iZone][INST_0].GetUnsteady_Simulation() == HARMONIC_BALANCE)
        FluidHBUpdate(iZone,FLOW_SOL);
    else if (config_container[iZone][INST_0].GetDynamic_Method() == MODAL_HARMONIC_BALANCE){
        cout<<"I was called HB Modal Update \n";
        FluidHBUpdate(iZone,MODAL_SOL);
    }
    else
        continue;
    }


    /*--- This is temporary. Each zone has to be monitored independently. Right now, fixes CHT output. ---*/

    Monitor(iOuter_Iter);

    Output(iOuter_Iter);

    Convergence = OuterConvergence(iOuter_Iter);
    cout << "convergence: " << Convergence << endl;
//     if (Convergence) break;

  }

}

void CMultizoneDriver::Run_Jacobi() {

  unsigned long iOuter_Iter;
  unsigned short jZone, UpdateMesh;
  bool DeformMesh = false;
  unsigned long ExtIter = 0;
  bool Convergence = false;

    cout << "\n\n run jacobi\n\n";
    
  unsigned long OuterIter = 0; for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetOuterIter(OuterIter);

  /*--- Loop over the number of outer iterations ---*/
  for (iOuter_Iter = 0; iOuter_Iter < driver_config->GetnOuter_Iter(); iOuter_Iter++){

          cout << "\n\n run jacobi iOuter --> " << iOuter_Iter <<"\n\n";

    /*--- Transfer from all zones ---*/
    for (iZone = 0; iZone < nZone; iZone++){
          cout << "\n\n run jacobi zone --> " << iZone <<"\n\n";

      /*--- In principle, the mesh does not need to be updated ---*/
      UpdateMesh = 0;

      /*--- Set the OuterIter ---*/
      config_container[iZone]->SetOuterIter(iOuter_Iter);

      /*--- Transfer from all the remaining zones ---*/
      for (jZone = 0; jZone < nZone; jZone++){
        /*--- The target zone is iZone ---*/
        if (jZone != iZone && transfer_container[iZone][jZone][INST_0] != NULL){
                      cout << "\n\n run jacobi deform zone --> " << jZone <<"\n\n";
          DeformMesh = Transfer_Data(jZone, iZone);
          if (DeformMesh) UpdateMesh+=1;
        }
      }
      /*--- If a mesh update is required due to the transfer of data ---*/
    cout << "\n\n run jacobi DynamicMeshUpdate --> " << iZone <<"\n\n";
      if (UpdateMesh > 0) DynamicMeshUpdate(iZone, ExtIter);

    }

      /*--- Loop over the number of zones (IZONE) ---*/
    for (iZone = 0; iZone < nZone; iZone++){
    cout << "\n\n run jacobi zone2 --> " << iZone <<"\n\n";
      /*--- Set the OuterIter ---*/
      config_container[iZone]->SetOuterIter(iOuter_Iter);
        cout << "\n\n run jacobi iteration->solve --> " << iZone <<"\n\n";

      /*--- Iterate the zone as a block, either to convergence or to a max number of iterations ---*/
      iteration_container[iZone][INST_0]->Solve(output, integration_container, geometry_container, solver_container,numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone, INST_0);

      /*--- A corrector step can help preventing numerical instabilities ---*/
        cout << "\n\n run jacobi iteration->corrector --> " << iZone <<"\n\n";
      Corrector(iZone);

    }
    cout << "\n\n run jacobi iteration->convergence\n\n";
    Convergence = OuterConvergence(iOuter_Iter);

    if (Convergence) break;

  }

}

void CMultizoneDriver::Corrector(unsigned short val_iZone) {

    if (config_container[val_iZone]->GetRelaxation())
      iteration_container[val_iZone][INST_0]->Relaxation(output, integration_container, geometry_container, solver_container,
          numerics_container, config_container, surface_movement, grid_movement, FFDBox, val_iZone, INST_0);

}

bool CMultizoneDriver::OuterConvergence(unsigned long OuterIter) {

  bool Convergence = false;
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  unsigned short iRes, iVar;
  unsigned short nVarSol;

  /*--- Compute the residual for the all the zones ---*/
  for (iZone = 0; iZone < nZone; iZone++){
    iVar = 0; // Initialize the variable index for each zone
//     if(config_container[iZone]->GetKind_Solver() == FEM_MODAL) return 1;
    
    /*--- Account for all the solvers ---*/
    for (unsigned short iSol = 0; iSol < MAX_SOLS; iSol++){
      /*-- If the solver position iSol is enabled --*/
      if (solver_container[iZone][INST_0][MESH_0][iSol] != NULL){
        nVarSol = solver_container[iZone][INST_0][MESH_0][iSol]->GetnVar();
//         cout<<"Compute Residul multizone being called"<<endl;
        /*--- Compute the block residual on each solver ---*/
        solver_container[iZone][INST_0][MESH_0][iSol]->ComputeResidual_Multizone(geometry_container[iZone][INST_0][MESH_0],config_container[iZone]);

        /*--- Loop over all the variables in the solver ---*/
        for (iRes = 0; iRes < nVarSol; iRes++){
          /*--- Store the log10 of the residual value ---*/
//          cout<< "iZone" << iZone << "iVar" << iVar <<"iRes"<<iRes<<endl;
//          cout<<log10(solver_container[iZone][INST_0][MESH_0][iSol]->GetRes_BGS(0))<<endl;
//          cout<<log10(solver_container[iZone][INST_0][MESH_0][iSol]->GetRes_BGS(1))<<endl;
//          cout<<log10(solver_container[iZone][INST_0][MESH_0][iSol]->GetRes_BGS(2))<<endl;
          residual[iZone][iVar] = log10(solver_container[iZone][INST_0][MESH_0][iSol]->GetRes_BGS(iRes));

          /*--- If it is the first iteration, the init_res is the current residual ---*/
          if (OuterIter == 0) init_res[iZone][iVar] = residual[iZone][iVar];
          /*--- residual_rel checks the difference in order of magnitude between the current and the initial residual ---*/
          residual_rel[iZone][iVar] = fabs(residual[iZone][iRes] - init_res[iZone][iRes]);
          /*--- Move the position of the container ---*/
          iVar++;
        }
      }
    }
  }

  /*--- This still has to be generalised ---*/
  if (fsi && (config_container[0]->GetKind_Solver() == FEM_ELASTICITY ||
      config_container[1]->GetKind_Solver() == FEM_ELASTICITY)){

    /*--- This should be possible to change dinamically in the config ---*/
    if (rank == MASTER_NODE){

      cout << endl << "-------------------------------------------------------------------------" << endl;
      cout << endl;
      cout << "Convergence summary for BGS iteration ";
      cout << OuterIter << endl;
      cout << endl;
      /*--- TODO: This is a workaround until the TestCases.py script incorporates new classes for nested loops. ---*/
      cout << "Iter[ID]" << "    BGSRes[Rho]" << "   BGSRes[RhoE]" << "     BGSRes[Ux]" << "     BGSRes[Uy]" << endl;
      cout.precision(6); cout.setf(ios::fixed, ios::floatfield);
      cout.width(8); cout << OuterIter*1000;
      cout.width(15); cout << residual[ZONE_0][0];
      cout.width(15); cout << residual[ZONE_0][3];
      cout.width(15); cout << residual[ZONE_1][0];
      cout.width(15); cout << residual[ZONE_1][1];
      cout << endl;
    }
    for (iZone = 0; iZone < nZone; iZone++){
      if (config_container[iZone]->GetKind_Solver() == FEM_ELASTICITY){
        integration_container[iZone][INST_0][FEA_SOL]->Convergence_Monitoring_FSI(geometry_container[iZone][INST_0][MESH_0], config_container[iZone], solver_container[iZone][INST_0][MESH_0][FEA_SOL], OuterIter);
        Convergence = integration_container[iZone][INST_0][FEA_SOL]->GetConvergence_FSI();
      }
    }
  }
  /*--- Update the residual for the all the zones ---*/
  for (iZone = 0; iZone < nZone; iZone++){
    if(config_container[iZone]->GetKind_Solver() == FEM_MODAL) break;
    /*--- Accounting for all the solvers ---*/
    for (unsigned short iSol = 0; iSol < MAX_SOLS; iSol++){
      /*-- If the solver position iSol is enabled --*/
      if (solver_container[iZone][INST_0][MESH_0][iSol] != NULL){
        solver_container[iZone][INST_0][MESH_0][iSol]->UpdateSolution_BGS(geometry_container[iZone][INST_0][MESH_0],
            config_container[iZone]);}
    }
  }

  if (rank == MASTER_NODE) cout.setf(ios::scientific, ios::floatfield);

  /*-----------------------------------------------------------------*/
  /*-------------------- Output FSI history -------------------------*/
  /*-----------------------------------------------------------------*/
  if (fsi && (config_container[0]->GetKind_Solver() == FEM_ELASTICITY ||
      config_container[1]->GetKind_Solver() == FEM_ELASTICITY)){
    bool ZONE_FLOW=0, ZONE_FEA=1;
    /*--- This is a hack to test it works. ---*/
    for (iZone = 0; iZone < nZone; iZone++){
      if (config_container[iZone]->GetKind_Solver() == FEM_ELASTICITY) ZONE_FEA = iZone;
      if (config_container[iZone]->GetKind_Solver() == NAVIER_STOKES) ZONE_FLOW = iZone;
    }
    output->SpecialOutput_FSI(&FSIHist_file, geometry_container, solver_container,
        config_container, integration_container, 0,
        ZONE_FLOW, ZONE_FEA, false);
  }
//     cout << "end convergence\n";
  return Convergence;

}

void CMultizoneDriver::Update() {

  unsigned short jZone, UpdateMesh;
  bool DeformMesh = false;
    cout << "\n\nMZ driver update\n\n";
  /*--- For enabling a consistent restart, we need to update the mesh with the interface information that introduces displacements --*/
  /*--- Loop over the number of zones (IZONE) ---*/
  for (iZone = 0; iZone < nZone; iZone++){

    UpdateMesh = 0;

      /*--- Transfer from all the remaining zones (JZONE != IZONE)---*/
      for (jZone = 0; jZone < nZone; jZone++){
        /*--- The target zone is iZone ---*/
        if (jZone != iZone){
          DeformMesh = Transfer_Data(jZone, iZone);
          if (DeformMesh) UpdateMesh += 1;
        }
      }
    /*--- If a mesh update is required due to the transfer of data ---*/
    if (UpdateMesh > 0) DynamicMeshUpdate(iZone, ExtIter);

    // Update harmonic Balance
    cout<<"Modal :: "<<(config_container[iZone][INST_0].GetDynamic_Method() == MODAL_HARMONIC_BALANCE)<<
    "Fluid :: "<<(config_container[iZone][INST_0].GetUnsteady_Simulation() == HARMONIC_BALANCE)<<endl;
    if (config_container[iZone][INST_0].GetUnsteady_Simulation() == HARMONIC_BALANCE)
        FluidHBUpdate(iZone,FLOW_SOL);
    else if (config_container[iZone][INST_0].GetDynamic_Method() == MODAL_HARMONIC_BALANCE)
        FluidHBUpdate(iZone,MODAL_SOL);
    else
        continue;

    for(iInst=0; iInst<nInst[iZone]; iInst++)
    iteration_container[iZone][iInst]->Update(output, integration_container, geometry_container,
        solver_container, numerics_container, config_container,
        surface_movement, grid_movement, FFDBox, iZone, iInst);

    /*--- Set the Convergence_FSI boolean to false for the next time step ---*/
    for (unsigned short iSol = 0; iSol < MAX_SOLS; iSol++){
      if (solver_container[iZone][INST_0][MESH_0][iSol] != NULL)
        integration_container[iZone][INST_0][iSol]->SetConvergence_FSI(false);
    }
  }
  
  cout << "MZ driver end update\n\n";

}

void CMultizoneDriver::Output(unsigned long TimeIter) {

  bool output_files = false;

  /*--- Determine whether a solution needs to be written
   after the current iteration ---*/
  if (

      /*--- General if statements to print output statements ---*/

      (TimeIter+1 >= config_container[ZONE_0]->GetnTime_Iter()) || (StopCalc) ||

      /*--- Unsteady problems ---*/

      (((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
        (config_container[ZONE_0]->GetUnsteady_Simulation() == TIME_STEPPING)) &&
       ((TimeIter == 0) || (ExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0))) ||

      ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND) &&
       ((TimeIter == 0) || ((TimeIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0) ||
                           ((TimeIter-1) % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0)))) ||

      ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND) &&
       ((TimeIter == 0) || ((TimeIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0)))) ||

      ((config_container[ZONE_0]->GetDynamic_Analysis() == DYNAMIC) &&
       ((TimeIter == 0) || (TimeIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0))) ||

      /*--- No inlet profile file found. Print template. ---*/

      (config_container[ZONE_0]->GetWrt_InletFile())

      ) {

    output_files = true;

  }

  /*--- Determine whether a solution doesn't need to be written
   after the current iteration ---*/

  if (config_container[ZONE_0]->GetFixed_CL_Mode()) {
    if (config_container[ZONE_0]->GetnExtIter()-config_container[ZONE_0]->GetIter_dCL_dAlpha() - 1 < ExtIter) output_files = false;
    if (config_container[ZONE_0]->GetnExtIter() - 1 == ExtIter) output_files = true;
  }

  /*--- write the solution ---*/

  if (output_files) {

    /*--- Time the output for performance benchmarking. ---*/
#ifndef HAVE_MPI
    StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
    StopTime = MPI_Wtime();
#endif
    UsedTimeCompute += StopTime-StartTime;
#ifndef HAVE_MPI
    StartTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
    StartTime = MPI_Wtime();
#endif

    if (rank == MASTER_NODE) cout << endl << "-------------------------- File Output Summary --------------------------\n";

    /*--- Execute the routine for writing restart, volume solution,
     surface solution, and surface comma-separated value files. ---*/
    cout << "driver MZ" << endl;
     output->SetResult_Files_Parallel(solver_container, geometry_container, config_container, TimeIter, nZone);


    /*--- Execute the routine for writing special output. ---*/
    output->SetSpecial_Output(solver_container, geometry_container, config_container, TimeIter, nZone);


    if (rank == MASTER_NODE) cout << "-------------------------------------------------------------------------" << endl << endl;

    /*--- Store output time and restart the timer for the compute phase. ---*/
#ifndef HAVE_MPI
    StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
    StopTime = MPI_Wtime();
#endif
    UsedTimeOutput += StopTime-StartTime;
    OutputCount++;
    BandwidthSum = config_container[ZONE_0]->GetRestart_Bandwidth_Agg();
#ifndef HAVE_MPI
    StartTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
    StartTime = MPI_Wtime();
#endif

  }

}

void CMultizoneDriver::DynamicMeshUpdate(unsigned long ExtIter) {

  bool harmonic_balance;
  for (iZone = 0; iZone < nZone; iZone++) {
   harmonic_balance = (config_container[iZone]->GetUnsteady_Simulation() == HARMONIC_BALANCE);
    /*--- Dynamic mesh update ---*/
    if ((config_container[iZone]->GetGrid_Movement()) && (!harmonic_balance) && (!fsi)) {
      iteration_container[iZone][INST_0]->SetGrid_Movement(geometry_container, surface_movement, grid_movement, FFDBox, solver_container, config_container, iZone, INST_0, 0, ExtIter );
    }

  }
}

void CMultizoneDriver::DynamicMeshUpdate(unsigned short val_iZone, unsigned long ExtIter) {
    for(iInst = 0; iInst < nInst[val_iZone]; iInst++) {
        iteration_container[val_iZone][iInst]->SetGrid_Movement(geometry_container, surface_movement, grid_movement, FFDBox,
                                                                solver_container, config_container, val_iZone, iInst, 0,
                                                                ExtIter);
    }
}

bool CMultizoneDriver::Transfer_Data(unsigned short donorZone, unsigned short targetZone) {

  bool UpdateMesh = false;

  /*--- Select the transfer method according to the magnitudes being transferred ---*/
    cout << "\n\ntransfer data: updatemesh = " << UpdateMesh << "\t" <<
    " solver: " << config_container[donorZone]->GetKind_Solver() << "\t" <<
    " transfer: " << transfer_types[donorZone][targetZone] << endl;

  if (transfer_types[donorZone][targetZone] == SLIDING_INTERFACE) {

    transfer_container[donorZone][targetZone][INST_0]->Broadcast_InterfaceData(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
                                                                       geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
                                                                       config_container[donorZone], config_container[targetZone]);
    if (config_container[targetZone]->GetKind_Solver() == RANS)
      transfer_container[donorZone][targetZone][INST_0]->Broadcast_InterfaceData(solver_container[donorZone][INST_0][MESH_0][TURB_SOL],solver_container[targetZone][INST_0][MESH_0][TURB_SOL],
                                                                         geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
                                                                         config_container[donorZone], config_container[targetZone]);
  }
  else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_FS) {
    transfer_container[donorZone][targetZone][INST_0]->Broadcast_InterfaceData(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][HEAT_SOL],
                                                                       geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
                                                                       config_container[donorZone], config_container[targetZone]);
  }
  else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_WEAKLY_FS) {
    transfer_container[donorZone][targetZone][INST_0]->Broadcast_InterfaceData(solver_container[donorZone][INST_0][MESH_0][HEAT_SOL],solver_container[targetZone][INST_0][MESH_0][HEAT_SOL],
                                                                       geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
                                                                       config_container[donorZone], config_container[targetZone]);
  }
  else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_SF) {
    transfer_container[donorZone][targetZone][INST_0]->Broadcast_InterfaceData(solver_container[donorZone][INST_0][MESH_0][HEAT_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
                                                                       geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
                                                                       config_container[donorZone], config_container[targetZone]);
  }
  else if (transfer_types[donorZone][targetZone] == CONJUGATE_HEAT_WEAKLY_SF) {
    transfer_container[donorZone][targetZone][INST_0]->Broadcast_InterfaceData(solver_container[donorZone][INST_0][MESH_0][HEAT_SOL],solver_container[targetZone][INST_0][MESH_0][HEAT_SOL],
                                                                       geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
                                                                       config_container[donorZone], config_container[targetZone]);
  }
  else if (transfer_types[donorZone][targetZone] == STRUCTURAL_DISPLACEMENTS){
      if(config_container[donorZone]->GetKind_Solver() == FEM_MODAL) {
        cout << "\n\ntransfer data: CSD structural displacements\n Disp \t\t Vels\n" << endl;
        for (unsigned short iMode=0; iMode < 4; ++iMode){
            cout << solver_container[donorZone][INST_0][MESH_0][MODAL_SOL]->getGeneralizedDisplacement(iMode) <<
            "\t" << solver_container[donorZone][INST_0][MESH_0][MODAL_SOL]->getGeneralizedVelocity(iMode) << endl;
        }

        for (iInst = 0; iInst<nInst[donorZone]; iInst++){
        transfer_container[donorZone][targetZone][iInst]->Broadcast_InterfaceData(solver_container[donorZone][iInst][MESH_0][MODAL_SOL],solver_container[targetZone][iInst][MESH_0][FLOW_SOL],
                geometry_container[donorZone][iInst][MESH_0],geometry_container[targetZone][iInst][MESH_0],
                config_container[donorZone], config_container[targetZone]);
        }
          UpdateMesh = true;
    }
    else{
          cout << "\n\ntransfer data: FEM structural displacements" << endl;
        transfer_container[donorZone][targetZone][INST_0]->Broadcast_InterfaceData(solver_container[donorZone][INST_0][MESH_0][FEA_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
        geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
        config_container[donorZone], config_container[targetZone]);
      UpdateMesh = true;  
    }
  }
  
  else if (transfer_types[donorZone][targetZone] == STRUCTURAL_DISPLACEMENTS && config_container[donorZone]->GetKind_Solver() != FEM_ELASTICITY) {
    transfer_container[donorZone][targetZone][INST_0]->Broadcast_InterfaceData(solver_container[donorZone][INST_0][MESH_0][FEA_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],config_container[donorZone], config_container[targetZone]);
    UpdateMesh = true;
  }
  else if (transfer_types[donorZone][targetZone] == FLOW_TRACTION) {
                cout << "\n\ntransfer data: fluid forces" << endl;
        for (iInst=0; iInst<nInst[donorZone]; iInst++) {
            cout<<transfer_container[donorZone][targetZone][iInst]<<endl;
            transfer_container[donorZone][targetZone][iInst]->Broadcast_InterfaceData(
                    solver_container[donorZone][iInst][MESH_0][FLOW_SOL],
                    solver_container[targetZone][iInst][MESH_0][FEA_SOL],
                    geometry_container[donorZone][iInst][MESH_0],
                    geometry_container[targetZone][iInst][MESH_0],
                    config_container[donorZone], config_container[targetZone]);
        }
  }
  else if ((transfer_types[donorZone][targetZone] == NO_TRANSFER)
           || (transfer_types[donorZone][targetZone] == ZONES_ARE_EQUAL)
           || (transfer_types[donorZone][targetZone] == NO_COMMON_INTERFACE)) { }
  else {
    cout << "WARNING: One of the intended interface transfer routines is not known to the chosen driver and has not been executed." << endl;
  }

  return UpdateMesh;
}

void CMultizoneDriver::FluidHBUpdate(unsigned short val_iZone, unsigned short val_Sol){
    if (val_Sol == FLOW_SOL)
        cout<<"Fluid_Harmonic_Balance Update \n";
    else if (val_Sol == MODAL_SOL)
        cout<<"Modal Hamronic Balance Update \n";

    for (iInst = 0; iInst < nInst[val_iZone]; iInst++) {
        /*--- Compute the harmonic balance terms across all zones ---*/
        SetHarmonicBalance(val_iZone, iInst, val_Sol);

    }

    /*--- Precondition the harmonic balance source terms ---*/
    if (config_container[val_iZone]->GetHB_Precondition() == YES) {
        StabilizeHarmonicBalance(val_iZone, val_Sol);

    }

    for (iInst = 0; iInst < nInst[val_iZone]; iInst++) {

        /*--- Update the harmonic balance terms across all zones ---*/
        iteration_container[val_iZone][iInst]->Update(output, integration_container, geometry_container,
                                                   solver_container, numerics_container, config_container,
                                                   surface_movement, grid_movement, FFDBox, val_iZone, iInst);

    }


}

void CMultizoneDriver::ModalHBUpdate(unsigned short val_iZone) {
    cout<<"Modal_Harmonic_Balance Update \n";
}

void CMultizoneDriver::SetHarmonicBalance(unsigned short val_iZone, unsigned short iInst, unsigned short val_Sol) {

    cout<<"++++++ Setting Harmonic Balance +++++ \n";

    unsigned short iVar, jInst, iMGlevel, val_Sol_Adj;

    if (val_Sol == FLOW_SOL)
        val_Sol_Adj = ADJFLOW_SOL;
    else if (val_Sol == MODAL_SOL)
        val_Sol_Adj = ADJFLOW_SOL;
    else
        exit(1);

    unsigned short nVar = solver_container[val_iZone][INST_0][MESH_0][val_Sol]->GetnVar();
    unsigned long iPoint;

    bool implicit = (config_container[val_iZone]->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    bool adjoint = (config_container[val_iZone]->GetContinuous_Adjoint());
    if (adjoint) {
        implicit = (config_container[val_iZone]->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
    }

    unsigned long ExtIter = config_container[val_iZone]->GetExtIter();

    /*--- Retrieve values from the config file ---*/
    su2double *U = new su2double[nVar];
    su2double *U_old = new su2double[nVar];
    su2double *Psi = new su2double[nVar];
    su2double *Psi_old = new su2double[nVar];
    su2double *Source = new su2double[nVar];
    su2double deltaU, deltaPsi;

    /*--- Compute period of oscillation ---*/
    su2double period = config_container[val_iZone]->GetHarmonicBalance_Period();

    /*--- Non-dimensionalize the input period, if necessary.  */
    period /= config_container[val_iZone]->GetTime_Ref();

    if (ExtIter == 0)
        ComputeHB_Operator(val_iZone);

    /*--- Compute various source terms for explicit direct, implicit direct, and adjoint problems ---*/
    /*--- Loop over all grid levels ---*/
    for (iMGlevel = 0; iMGlevel <= config_container[val_iZone]->GetnMGLevels(); iMGlevel++) {

        /*--- Loop over each node in the volume mesh ---*/
        for (iPoint = 0; iPoint < geometry_container[val_iZone][iInst][iMGlevel]->GetnPoint(); iPoint++) {

            for (iVar = 0; iVar < nVar; iVar++) {
                Source[iVar] = 0.0;
            }

            /*--- Step across the columns ---*/
            for (jInst = 0; jInst < nInst[val_iZone]; jInst++) {

                /*--- Retrieve solution at this node in current zone ---*/
                for (iVar = 0; iVar < nVar; iVar++) {

                    if (!adjoint) {
                        U[iVar] = solver_container[val_iZone][jInst][iMGlevel][val_Sol]->node[iPoint]->GetSolution(iVar);
                        Source[iVar] += U[iVar]*D[iInst][jInst];

                        if (implicit) {
                            U_old[iVar] = solver_container[val_iZone][jInst][iMGlevel][val_Sol]->node[iPoint]->GetSolution_Old(iVar);
                            deltaU = U[iVar] - U_old[iVar];
                            Source[iVar] += deltaU*D[iInst][jInst];
                        }

                    }

                    else {
                        Psi[iVar] = solver_container[val_iZone][jInst][iMGlevel][val_Sol_Adj]->node[iPoint]->GetSolution(iVar);
                        Source[iVar] += Psi[iVar]*D[jInst][iInst];

                        if (implicit) {
                            Psi_old[iVar] = solver_container[val_iZone][jInst][iMGlevel][val_Sol_Adj]->node[iPoint]->GetSolution_Old(iVar);
                            deltaPsi = Psi[iVar] - Psi_old[iVar];
                            Source[iVar] += deltaPsi*D[jInst][iInst];
                        }
                    }
//                    cout<<"HB Source Term :: "<<U[iVar]<<" "<<D[iInst][jInst]<<" "<<Source[iVar]<<endl;
                }

                /*--- Store sources for current row ---*/
                for (iVar = 0; iVar < nVar; iVar++) {
                    if (!adjoint) {
                        solver_container[val_iZone][iInst][iMGlevel][val_Sol]->node[iPoint]->SetHarmonicBalance_Source(iVar, Source[iVar]);
                        //cout<<solver_container[val_iZone][iInst][iMGlevel][val_Sol]->node[iPoint]->GetHarmonicBalance_Source(iVar);
                    }
                    else {
                        solver_container[val_iZone][iInst][iMGlevel][val_Sol_Adj]->node[iPoint]->SetHarmonicBalance_Source(iVar, Source[iVar]);
                    }
                }

            }
        }
    }

    /*--- Source term for a turbulence model ---*/
    if (config_container[val_iZone]->GetKind_Solver() == RANS) {

        /*--- Extra variables needed if we have a turbulence model. ---*/
        unsigned short nVar_Turb = solver_container[val_iZone][INST_0][MESH_0][TURB_SOL]->GetnVar();
        su2double *U_Turb = new su2double[nVar_Turb];
        su2double *Source_Turb = new su2double[nVar_Turb];

        /*--- Loop over only the finest mesh level (turbulence is always solved
         on the original grid only). ---*/
        for (iPoint = 0; iPoint < geometry_container[val_iZone][INST_0][MESH_0]->GetnPoint(); iPoint++) {
            for (iVar = 0; iVar < nVar_Turb; iVar++) Source_Turb[iVar] = 0.0;
            for (jInst = 0; jInst < nInst[val_iZone]; jInst++) {

                /*--- Retrieve solution at this node in current zone ---*/
                for (iVar = 0; iVar < nVar_Turb; iVar++) {
                    U_Turb[iVar] = solver_container[val_iZone][jInst][MESH_0][TURB_SOL]->node[iPoint]->GetSolution(iVar);
                    Source_Turb[iVar] += U_Turb[iVar]*D[iInst][jInst];
                }
            }

            /*--- Store sources for current iZone ---*/
            for (iVar = 0; iVar < nVar_Turb; iVar++)
                solver_container[val_iZone][iInst][MESH_0][TURB_SOL]->node[iPoint]->SetHarmonicBalance_Source(iVar, Source_Turb[iVar]);
        }

        delete [] U_Turb;
        delete [] Source_Turb;
    }

    delete [] Source;
    delete [] U;
    delete [] U_old;
    delete [] Psi;
    delete [] Psi_old;

}

void CMultizoneDriver::SetHarmonicBalance_Modal(unsigned short val_iZone, unsigned short iInst, unsigned short val_Sol) {

    cout<<"++++++ Setting Harmonic Balance +++++ \n";

    unsigned short iVar, jInst, iMGlevel, val_Sol_Adj;

    if (val_Sol == FLOW_SOL)
        val_Sol_Adj = ADJFLOW_SOL;
    else if (val_Sol == MODAL_SOL)
        val_Sol_Adj = ADJFLOW_SOL;
    else
        exit(1);

    unsigned short nVar = solver_container[val_iZone][INST_0][MESH_0][val_Sol]->GetnVar();
    unsigned long iPoint;

    bool implicit = (config_container[val_iZone]->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    bool adjoint = (config_container[val_iZone]->GetContinuous_Adjoint());
    if (adjoint) {
        implicit = (config_container[val_iZone]->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
    }

    unsigned long ExtIter = config_container[val_iZone]->GetExtIter();

    /*--- Retrieve values from the config file ---*/
    su2double *U = new su2double[nVar];
    su2double *U_old = new su2double[nVar];
    su2double *Psi = new su2double[nVar];
    su2double *Psi_old = new su2double[nVar];
    su2double *Source = new su2double[nVar];
    su2double deltaU, deltaPsi;

    /*--- Compute period of oscillation ---*/
    su2double period = config_container[val_iZone]->GetHarmonicBalance_Period();

    /*--- Non-dimensionalize the input period, if necessary.  */
    period /= config_container[val_iZone]->GetTime_Ref();

    if (ExtIter == 0)
        ComputeHB_Operator(val_iZone);

    /*--- Compute various source terms for explicit direct, implicit direct, and adjoint problems ---*/
    /*--- Loop over all grid levels ---*/
    for (iMGlevel = 0; iMGlevel <= config_container[val_iZone]->GetnMGLevels(); iMGlevel++) {

            for (iVar = 0; iVar < nVar; iVar++) {
                Source[iVar] = 0.0;
            }

            /*--- Step across the columns ---*/
            for (jInst = 0; jInst < nInst[val_iZone]; jInst++) {

                /*--- Retrieve solution at this node in current zone ---*/
                for (iVar = 0; iVar < nVar; iVar++) {

                    if (!adjoint) {
                        U[iVar] = solver_container[val_iZone][jInst][iMGlevel][val_Sol]->Get_QSol(iVar);
                        Source[iVar] += U[iVar]*D[iInst][jInst];

                        if (implicit) {
                            U_old[iVar] = solver_container[val_iZone][jInst][iMGlevel][val_Sol]->node[iPoint]->GetSolution_Old(iVar);
                            deltaU = U[iVar] - U_old[iVar];
                            Source[iVar] += deltaU*D[iInst][jInst];
                        }

                    }

                    else {
                        Psi[iVar] = solver_container[val_iZone][jInst][iMGlevel][val_Sol_Adj]->node[iPoint]->GetSolution(iVar);
                        Source[iVar] += Psi[iVar]*D[jInst][iInst];

                        if (implicit) {
                            Psi_old[iVar] = solver_container[val_iZone][jInst][iMGlevel][val_Sol_Adj]->node[iPoint]->GetSolution_Old(iVar);
                            deltaPsi = Psi[iVar] - Psi_old[iVar];
                            Source[iVar] += deltaPsi*D[jInst][iInst];
                        }
                    }
//                    cout<<"HB Source Term :: "<<U[iVar]<<" "<<D[iInst][jInst]<<" "<<Source[iVar]<<endl;
                }

                /*--- Store sources for current row ---*/
                for (iVar = 0; iVar < nVar; iVar++) {
                    if (!adjoint) {
                        solver_container[val_iZone][iInst][iMGlevel][val_Sol]->node[iPoint]->SetHarmonicBalance_Source(iVar, Source[iVar]);
                        //cout<<solver_container[val_iZone][iInst][iMGlevel][val_Sol]->node[iPoint]->GetHarmonicBalance_Source(iVar);
                    }
                    else {
                        solver_container[val_iZone][iInst][iMGlevel][val_Sol_Adj]->node[iPoint]->SetHarmonicBalance_Source(iVar, Source[iVar]);
                    }
                }

            }
    }

    /*--- Source term for a turbulence model ---*/
    if (config_container[val_iZone]->GetKind_Solver() == RANS) {

        /*--- Extra variables needed if we have a turbulence model. ---*/
        unsigned short nVar_Turb = solver_container[val_iZone][INST_0][MESH_0][TURB_SOL]->GetnVar();
        su2double *U_Turb = new su2double[nVar_Turb];
        su2double *Source_Turb = new su2double[nVar_Turb];

        /*--- Loop over only the finest mesh level (turbulence is always solved
         on the original grid only). ---*/
        for (iPoint = 0; iPoint < geometry_container[val_iZone][INST_0][MESH_0]->GetnPoint(); iPoint++) {
            for (iVar = 0; iVar < nVar_Turb; iVar++) Source_Turb[iVar] = 0.0;
            for (jInst = 0; jInst < nInst[val_iZone]; jInst++) {

                /*--- Retrieve solution at this node in current zone ---*/
                for (iVar = 0; iVar < nVar_Turb; iVar++) {
                    U_Turb[iVar] = solver_container[val_iZone][jInst][MESH_0][TURB_SOL]->node[iPoint]->GetSolution(iVar);
                    Source_Turb[iVar] += U_Turb[iVar]*D[iInst][jInst];
                }
            }

            /*--- Store sources for current iZone ---*/
            for (iVar = 0; iVar < nVar_Turb; iVar++)
                solver_container[val_iZone][iInst][MESH_0][TURB_SOL]->node[iPoint]->SetHarmonicBalance_Source(iVar, Source_Turb[iVar]);
        }

        delete [] U_Turb;
        delete [] Source_Turb;
    }

    delete [] Source;
    delete [] U;
    delete [] U_old;
    delete [] Psi;
    delete [] Psi_old;

}

void CMultizoneDriver::ComputeHB_Operator(unsigned short val_iZone) {

    const   complex<su2double> J(0.0,1.0);
    unsigned short i, j, k, iInst;

    su2double *Omega_HB       = new su2double[nInst[val_iZone]];
    complex<su2double> **E    = new complex<su2double>*[nInst[val_iZone]];
    complex<su2double> **Einv = new complex<su2double>*[nInst[val_iZone]];
    complex<su2double> **DD   = new complex<su2double>*[nInst[val_iZone]];
    for (iInst = 0; iInst < nInst[val_iZone]; iInst++) {
        E[iInst]    = new complex<su2double>[nInst[val_iZone]];
        Einv[iInst] = new complex<su2double>[nInst[val_iZone]];
        DD[iInst]   = new complex<su2double>[nInst[val_iZone]];
    }

    /*--- Get simualation period from config file ---*/
    su2double Period = config_container[val_iZone]->GetHarmonicBalance_Period();

    /*--- Non-dimensionalize the input period, if necessary.      */
//    Period /= config_container[val_iZone]->GetTime_Ref();

    /*--- Build the array containing the selected frequencies to solve ---*/
    for (iInst = 0; iInst < nInst[val_iZone]; iInst++) {
        Omega_HB[iInst]  = config_container[val_iZone]->GetOmega_HB()[iInst];
//        Omega_HB[iInst] /= config_container[val_iZone]->GetOmega_Ref(); //TODO: check
        cout<< "Omega HB = "<<Omega_HB[iInst]<<endl;
//        cout<< "Omega Ref = "<<config_container[val_iZone]->GetOmega_Ref()<<endl;
    }

    /*--- Build the diagonal matrix of the frequencies DD ---*/
    for (i = 0; i < nInst[val_iZone]; i++) {
        for (k = 0; k < nInst[val_iZone]; k++) {
            if (k == i ) {
                DD[i][k] = J*Omega_HB[k];
            }
        }
    }


    /*--- Build the harmonic balance inverse matrix ---*/
    for (i = 0; i < nInst[val_iZone]; i++) {
        for (k = 0; k < nInst[val_iZone]; k++) {
            Einv[i][k] = complex<su2double>(cos(Omega_HB[k]*(i*Period/nInst[val_iZone]))) + J*complex<su2double>(sin(Omega_HB[k]*(i*Period/nInst[val_iZone])));
        }
    }

    /*---  Invert inverse harmonic balance Einv with Gauss elimination ---*/

    /*--  A temporary matrix to hold the inverse, dynamically allocated ---*/
    complex<su2double> **temp = new complex<su2double>*[nInst[val_iZone]];
    for (i = 0; i < nInst[val_iZone]; i++) {
        temp[i] = new complex<su2double>[2 * nInst[val_iZone]];
    }

    /*---  Copy the desired matrix into the temporary matrix ---*/
    for (i = 0; i < nInst[val_iZone]; i++) {
        for (j = 0; j < nInst[val_iZone]; j++) {
            temp[i][j] = Einv[i][j];
            temp[i][nInst[val_iZone] + j] = 0;
        }
        temp[i][nInst[val_iZone] + i] = 1;
    }

    su2double max_val;
    unsigned short max_idx;

    /*---  Pivot each column such that the largest number possible divides the other rows  ---*/
    for (k = 0; k < nInst[val_iZone] - 1; k++) {
        max_idx = k;
        max_val = abs(temp[k][k]);
        /*---  Find the largest value (pivot) in the column  ---*/
        for (j = k; j < nInst[val_iZone]; j++) {
            if (abs(temp[j][k]) > max_val) {
                max_idx = j;
                max_val = abs(temp[j][k]);
            }
        }
        /*---  Move the row with the highest value up  ---*/
        for (j = 0; j < (nInst[val_iZone] * 2); j++) {
            complex<su2double> d = temp[k][j];
            temp[k][j] = temp[max_idx][j];
            temp[max_idx][j] = d;
        }
        /*---  Subtract the moved row from all other rows ---*/
        for (i = k + 1; i < nInst[val_iZone]; i++) {
            complex<su2double> c = temp[i][k] / temp[k][k];
            for (j = 0; j < (nInst[val_iZone] * 2); j++) {
                temp[i][j] = temp[i][j] - temp[k][j] * c;
            }
        }
    }
    /*---  Back-substitution  ---*/
    for (k = nInst[val_iZone] - 1; k > 0; k--) {
        if (temp[k][k] != complex<su2double>(0.0)) {
            for (int i = k - 1; i > -1; i--) {
                complex<su2double> c = temp[i][k] / temp[k][k];
                for (j = 0; j < (nInst[val_iZone] * 2); j++) {
                    temp[i][j] = temp[i][j] - temp[k][j] * c;
                }
            }
        }
    }
    /*---  Normalize the inverse  ---*/
    for (i = 0; i < nInst[val_iZone]; i++) {
        complex<su2double> c = temp[i][i];
        for (j = 0; j < nInst[val_iZone]; j++) {
            temp[i][j + nInst[val_iZone]] = temp[i][j + nInst[val_iZone]] / c;
        }
    }
    /*---  Copy the inverse back to the main program flow ---*/
    for (i = 0; i < nInst[val_iZone]; i++) {
        for (j = 0; j < nInst[val_iZone]; j++) {
            E[i][j] = temp[i][j + nInst[val_iZone]];
        }
    }
    /*---  Delete dynamic template  ---*/
    for (i = 0; i < nInst[val_iZone]; i++) {
        delete[] temp[i];
    }
    delete[] temp;


    /*---  Temporary matrix for performing product  ---*/
    complex<su2double> **Temp    = new complex<su2double>*[nInst[val_iZone]];

    /*---  Temporary complex HB operator  ---*/
    complex<su2double> **Dcpx    = new complex<su2double>*[nInst[val_iZone]];

    for (iInst = 0; iInst < nInst[val_iZone]; iInst++){
        Temp[iInst]    = new complex<su2double>[nInst[val_iZone]];
        Dcpx[iInst]   = new complex<su2double>[nInst[val_iZone]];
    }


    /*---  Calculation of the HB operator matrix ---*/
    for (int row = 0; row < nInst[val_iZone]; row++) {
        for (int col = 0; col < nInst[val_iZone]; col++) {
            for (int inner = 0; inner < nInst[val_iZone]; inner++) {
                Temp[row][col] += Einv[row][inner] * DD[inner][col];
            }
        }
    }

    unsigned short row, col, inner;

    for (row = 0; row < nInst[val_iZone]; row++) {
        for (col = 0; col < nInst[val_iZone]; col++) {
            for (inner = 0; inner < nInst[val_iZone]; inner++) {
                Dcpx[row][col] += Temp[row][inner] * E[inner][col];
            }
        }
    }

    /*---  Take just the real part of the HB operator matrix ---*/
    for (i = 0; i < nInst[val_iZone]; i++) {
        for (k = 0; k < nInst[val_iZone]; k++) {
            D[i][k] = real(Dcpx[i][k]);
        }
    }

    //    cout<<" +++ D Matrix +++\n";
    //    for(i=0; i<nInst[val_iZone]; i++) {
    //        for (k = 0; k < nInst[val_iZone]; k++)
    //            cout << D[i][k] << "\t";
    //        cout<<"\n";
    //    }

    /*--- Deallocate dynamic memory ---*/
    for (iInst = 0; iInst < nInst[val_iZone]; iInst++){
        delete [] E[iInst];
        delete [] Einv[iInst];
        delete [] DD[iInst];
        delete [] Temp[iInst];
        delete [] Dcpx[iInst];
    }
    delete [] E;
    delete [] Einv;
    delete [] DD;
    delete [] Temp;
    delete [] Dcpx;
    delete [] Omega_HB;

}

void CMultizoneDriver::StabilizeHarmonicBalance(unsigned short val_iZone, unsigned short val_Sol){

    unsigned short nInstHB = nInst[val_iZone];
    unsigned short i, j, k, iVar, iInst, jInst, iMGlevel, val_Sol_Adj;

    if (val_Sol == FLOW_SOL)
        val_Sol_Adj = ADJFLOW_SOL;
    else if (val_Sol == MODAL_SOL)
        val_Sol_Adj = ADJFLOW_SOL;
    else
        exit(1);

    unsigned short nVar = solver_container[val_iZone][INST_0][MESH_0][val_Sol]->GetnVar();
    unsigned long iPoint;
    bool adjoint = (config_container[val_iZone]->GetContinuous_Adjoint());

    /*--- Retrieve values from the config file ---*/
    su2double *Source     = new su2double[nInstHB];
    su2double *Source_old = new su2double[nInstHB];
    su2double Delta;

    su2double **Pinv     = new su2double*[nInstHB];
    su2double **P        = new su2double*[nInstHB];
    for (iInst = 0; iInst < nInstHB; iInst++) {
        Pinv[iInst]       = new su2double[nInstHB];
        P[iInst]          = new su2double[nInstHB];
    }

    /*--- Loop over all grid levels ---*/
    for (iMGlevel = 0; iMGlevel <= config_container[val_iZone]->GetnMGLevels(); iMGlevel++) {

        /*--- Loop over each node in the volume mesh ---*/
        for (iPoint = 0; iPoint < geometry_container[val_iZone][INST_0][iMGlevel]->GetnPoint(); iPoint++) {

            /*--- Get time step for current node ---*/
            Delta = solver_container[val_iZone][INST_0][iMGlevel][val_Sol]->node[iPoint]->GetDelta_Time();

            /*--- Setup stabilization matrix for this node ---*/
            for (iInst = 0; iInst < nInstHB; iInst++) {
                for (jInst = 0; jInst < nInstHB; jInst++) {
                    if (jInst == iInst ) {
                        Pinv[iInst][jInst] = 1.0 + Delta*D[iInst][jInst];
                    }
                    else {
                        Pinv[iInst][jInst] = Delta*D[iInst][jInst];
                    }
                }
            }

            /*--- Invert stabilization matrix Pinv with Gauss elimination---*/

            /*--  A temporary matrix to hold the inverse, dynamically allocated ---*/
            su2double **temp = new su2double*[nInstHB];
            for (i = 0; i < nInstHB; i++) {
                temp[i] = new su2double[2 * nInstHB];
            }

            /*---  Copy the desired matrix into the temporary matrix ---*/
            for (i = 0; i < nInstHB; i++) {
                for (j = 0; j < nInstHB; j++) {
                    temp[i][j] = Pinv[i][j];
                    temp[i][nInstHB + j] = 0;
                }
                temp[i][nInstHB + i] = 1;
            }

            su2double max_val;
            unsigned short max_idx;

            /*---  Pivot each column such that the largest number possible divides the other rows  ---*/
            for (k = 0; k < nInstHB - 1; k++) {
                max_idx = k;
                max_val = abs(temp[k][k]);
                /*---  Find the largest value (pivot) in the column  ---*/
                for (j = k; j < nInstHB; j++) {
                    if (abs(temp[j][k]) > max_val) {
                        max_idx = j;
                        max_val = abs(temp[j][k]);
                    }
                }

                /*---  Move the row with the highest value up  ---*/
                for (j = 0; j < (nInstHB * 2); j++) {
                    su2double d = temp[k][j];
                    temp[k][j] = temp[max_idx][j];
                    temp[max_idx][j] = d;
                }
                /*---  Subtract the moved row from all other rows ---*/
                for (i = k + 1; i < nInstHB; i++) {
                    su2double c = temp[i][k] / temp[k][k];
                    for (j = 0; j < (nInstHB * 2); j++) {
                        temp[i][j] = temp[i][j] - temp[k][j] * c;
                    }
                }
            }

            /*---  Back-substitution  ---*/
            for (k = nInstHB - 1; k > 0; k--) {
                if (temp[k][k] != su2double(0.0)) {
                    for (int i = k - 1; i > -1; i--) {
                        su2double c = temp[i][k] / temp[k][k];
                        for (j = 0; j < (nInstHB * 2); j++) {
                            temp[i][j] = temp[i][j] - temp[k][j] * c;
                        }
                    }
                }
            }

            /*---  Normalize the inverse  ---*/
            for (i = 0; i < nInstHB; i++) {
                su2double c = temp[i][i];
                for (j = 0; j < nInstHB; j++) {
                    temp[i][j + nInstHB] = temp[i][j + nInstHB] / c;
                }
            }

            /*---  Copy the inverse back to the main program flow ---*/
            for (i = 0; i < nInstHB; i++) {
                for (j = 0; j < nInstHB; j++) {
                    P[i][j] = temp[i][j + nInstHB];
                }
            }

            /*---  Delete dynamic template  ---*/
            for (iInst = 0; iInst < nInstHB; iInst++) {
                delete[] temp[iInst];
            }
            delete[] temp;

            /*--- Loop through variables to precondition ---*/
            for (iVar = 0; iVar < nVar; iVar++) {

                /*--- Get current source terms (not yet preconditioned) and zero source array to prepare preconditioning ---*/
                for (iInst = 0; iInst < nInstHB; iInst++) {
                    Source_old[iInst] = solver_container[ZONE_0][iInst][iMGlevel][val_Sol]->node[iPoint]->GetHarmonicBalance_Source(iVar);
                    Source[iInst] = 0;
                }

                /*--- Step through columns ---*/
                for (iInst = 0; iInst < nInstHB; iInst++) {
                    for (jInst = 0; jInst < nInstHB; jInst++) {
                        Source[iInst] += P[iInst][jInst]*Source_old[jInst];
                    }

                    /*--- Store updated source terms for current node ---*/
                    if (!adjoint) {
                        solver_container[val_iZone][iInst][iMGlevel][val_Sol]->node[iPoint]->SetHarmonicBalance_Source(iVar, Source[iInst]);
                    }
                    else {
                        solver_container[val_iZone][iInst][iMGlevel][val_Sol_Adj]->node[iPoint]->SetHarmonicBalance_Source(iVar, Source[iInst]);
                    }
                }

            }
        }
    }

    /*--- Deallocate dynamic memory ---*/
    for (iInst = 0; iInst < nInstHB; iInst++){
        delete [] P[iInst];
        delete [] Pinv[iInst];
    }
    delete [] P;
    delete [] Pinv;
    delete [] Source;
    delete [] Source_old;
}
