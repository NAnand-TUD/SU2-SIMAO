/*!
 * \file solver_direct_elasticity.cpp
 * \brief Main subroutines for solving direct FEM elasticity problems.
 * \author R. Sanchez
 * \version 6.2.0 "Falcon"
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

#include "../include/solver_structure.hpp"
#include <algorithm>

CFEASolver::CFEASolver(void) : CSolver() {
  
  nElement = 0;
  nDim = 0;
  nMarker = 0;
  
  nFEA_Terms = 1;
  
  nPoint = 0;
  nPointDomain = 0;
  
  Total_CFEA = 0.0;
  WAitken_Dyn = 0.0;
  WAitken_Dyn_tn1 = 0.0;
  loadIncrement = 1.0;
  
  element_container = NULL;
  node = NULL;
  
  element_properties = NULL;
  elProperties = NULL;
  
  GradN_X = NULL;
  GradN_x = NULL;
  
  Jacobian_c_ij = NULL;
  Jacobian_s_ij = NULL;
  Jacobian_k_ij = NULL;
  
  MassMatrix_ij = NULL;
  
  mZeros_Aux = NULL;
  mId_Aux = NULL;
  
  Res_Stress_i = NULL;
  Res_Ext_Surf = NULL;
  Res_Time_Cont = NULL;
  Res_FSI_Cont = NULL;
  
  Res_Dead_Load = NULL;
  
  nodeReactions = NULL;
  
  solutionPredictor = NULL;
  
  SolRest = NULL;
  
  normalVertex = NULL;
  stressTensor = NULL;
  
  Solution_Interm = NULL;
  
  iElem_iDe   = NULL;
  
}

CFEASolver::CFEASolver(CGeometry *geometry, CConfig *config) : CSolver() {
  
  unsigned long iPoint;
  unsigned short iVar, jVar, iDim, jDim;
  unsigned short iTerm, iKind;

  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
  bool fsi = config->GetFSI_Simulation();                        // FSI simulation
  bool gen_alpha = (config->GetKind_TimeIntScheme_FEA() == GENERALIZED_ALPHA);  // Generalized alpha method requires residual at previous time step.
  
  bool de_effects = config->GetDE_Effects();                      // Test whether we consider dielectric elastomers
  
  bool body_forces = config->GetDeadLoad();  // Body forces (dead loads).
  bool incompressible = (config->GetMaterialCompressibility() == INCOMPRESSIBLE_MAT);
  
  element_based = false;          // A priori we don't have an element-based input file (most of the applications will be like this)
  
  nElement      = geometry->GetnElem();
  nDim          = geometry->GetnDim();
  nMarker       = geometry->GetnMarker();
  
  nPoint        = geometry->GetnPoint();
  nPointDomain  = geometry->GetnPointDomain();
  
  /*--- Number of different terms for FEA ---*/
  nFEA_Terms = 1;
  if (de_effects) nFEA_Terms++;       // The DE term is DE_TERM = 1
  if (incompressible) nFEA_Terms = 3; // The incompressible term is INC_TERM = 2
  
  /*--- Here is where we assign the kind of each element ---*/
  
  /*--- First level: different possible terms of the equations ---*/
  element_container = new CElement** [MAX_TERMS];
  for (iTerm = 0; iTerm < MAX_TERMS; iTerm++)
    element_container[iTerm] = new CElement* [MAX_FE_KINDS];
  
  for (iTerm = 0; iTerm < MAX_TERMS; iTerm++) {
    for (iKind = 0; iKind < MAX_FE_KINDS; iKind++) {
      element_container[iTerm][iKind] = NULL;
    }
  }
  
  if (nDim == 2) {
    
      /*--- Basic terms ---*/
      element_container[FEA_TERM][EL_TRIA] = new CTRIA1(nDim, config);
      element_container[FEA_TERM][EL_QUAD] = new CQUAD4(nDim, config);
    
      if (de_effects){
        element_container[DE_TERM][EL_TRIA] = new CTRIA1(nDim, config);
        element_container[DE_TERM][EL_QUAD] = new CQUAD4(nDim, config);
      }

      if (incompressible){
        element_container[INC_TERM][EL_TRIA] = new CTRIA1(nDim, config);
        element_container[INC_TERM][EL_QUAD] = new CQUAD1(nDim, config);
      }

  }
  else if (nDim == 3) {

      element_container[FEA_TERM][EL_TETRA] = new CTETRA1(nDim, config);
      element_container[FEA_TERM][EL_HEXA] = new CHEXA8(nDim, config);

      if (de_effects){
        element_container[DE_TERM][EL_TETRA] = new CTETRA1(nDim, config);
        element_container[DE_TERM][EL_HEXA] = new CHEXA8(nDim, config);
      }

      if (incompressible) {
        element_container[INC_TERM][EL_TETRA] = new CTETRA1(nDim, config);
        element_container[INC_TERM][EL_HEXA] = new CHEXA1(nDim, config);
      }


  }
  
  node              = new CVariable*[nPoint];
  
  /*--- Set element properties ---*/
  elProperties = new unsigned long[4];
  for (iVar = 0; iVar < 4; iVar++)
    elProperties[iVar] = 0;
  Set_ElementProperties(geometry, config);
  
  GradN_X = new su2double [nDim];
  GradN_x = new su2double [nDim];
  
  Total_CFEA      = 0.0;
  WAitken_Dyn       = 0.0;
  WAitken_Dyn_tn1   = 0.0;
  loadIncrement     = 0.0;
  
  SetFSI_ConvValue(0,0.0);
  SetFSI_ConvValue(1,0.0);
  
  nVar = nDim;
  
  /*--- Define some auxiliary vectors related to the residual ---*/
  
  Residual = new su2double[nVar];          for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS = new su2double[nVar];      for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max = new su2double[nVar];      for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Point_Max = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  
  /*--- Residual i and residual j for the Matrix-Vector Product for the Mass Residual ---*/
  if (dynamic){
    Residual_i = new su2double[nVar];        for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
    Residual_j = new su2double[nVar];        for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
  }
  
  /*--- Define some auxiliary vectors related to the solution ---*/
  
  Solution   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
  
  Solution_Interm = NULL;
  if (gen_alpha) {
    Solution_Interm = new su2double[nVar];
    for (iVar = 0; iVar < nVar; iVar++) Solution_Interm[iVar] = 0.0;
  }
  
  nodeReactions = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) nodeReactions[iVar]   = 0.0;
  
  /*--- The length of the solution vector depends on whether the problem is static or dynamic ---*/
  
  unsigned short nSolVar;
  string text_line, filename;
  ifstream restart_file;

  if (dynamic) nSolVar = 3 * nVar;
  else nSolVar = nVar;
  
  SolRest = new su2double[nSolVar];

  /*--- Initialize from zero everywhere. ---*/

  for (iVar = 0; iVar < nSolVar; iVar++) SolRest[iVar] = 0.0;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint] = new CFEAVariable(SolRest, nDim, nVar, config);
  }
  
  bool reference_geometry = config->GetRefGeom();
  if (reference_geometry) Set_ReferenceGeometry(geometry, config);
  
  bool prestretch_fem = config->GetPrestretch();
  if (prestretch_fem) Set_Prestretch(geometry, config);  
  
  /*--- Term ij of the Jacobian ---*/
  
  Jacobian_ij = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_ij[iVar] = new su2double [nVar];
    for (jVar = 0; jVar < nVar; jVar++) {
      Jacobian_ij[iVar][jVar] = 0.0;
    }
  }
  
  /*--- Term ij of the Mass Matrix (only if dynamic analysis) ---*/
  MassMatrix_ij = NULL;
  if (dynamic) {
    MassMatrix_ij = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      MassMatrix_ij[iVar] = new su2double [nVar];
      for (jVar = 0; jVar < nVar; jVar++) {
        MassMatrix_ij[iVar][jVar] = 0.0;
      }
    }
  }
  
  Jacobian_c_ij = NULL;
  Jacobian_s_ij = NULL;
  if (nonlinear_analysis) {
    
    /*--- Term ij of the Jacobian (constitutive contribution) ---*/
    
    Jacobian_c_ij = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_c_ij[iVar] = new su2double [nVar];
      for (jVar = 0; jVar < nVar; jVar++) {
        Jacobian_c_ij[iVar][jVar] = 0.0;
      }
    }
    
    /*--- Term ij of the Jacobian (stress contribution) ---*/
    
    Jacobian_s_ij = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_s_ij[iVar] = new su2double [nVar];
      for (jVar = 0; jVar < nVar; jVar++) {
        Jacobian_s_ij[iVar][jVar] = 0.0;
      }
    }
    
  }
  
  /*--- Term ij of the Jacobian (incompressibility term) ---*/
  Jacobian_k_ij = NULL;
  if (incompressible) {
    Jacobian_k_ij = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_k_ij[iVar] = new su2double [nVar];
      for (jVar = 0; jVar < nVar; jVar++) {
        Jacobian_k_ij[iVar][jVar] = 0.0;
      }
    }
  }
  
  /*--- Stress contribution to the node i ---*/
  Res_Stress_i = new su2double[nVar];
  
  /*--- Contribution of the external surface forces to the residual (auxiliary vector) ---*/
  Res_Ext_Surf = new su2double[nVar];
  
  /*--- Contribution of the body forces to the residual (auxiliary vector) ---*/
  Res_Dead_Load = NULL;
  if (body_forces) {
    Res_Dead_Load = new su2double[nVar];
  }
  
  /*--- Contribution of the fluid tractions to the residual (auxiliary vector) ---*/
  Res_FSI_Cont = NULL;
  if (fsi) {
    Res_FSI_Cont = new su2double[nVar];
  }
  
  /*--- Time integration contribution to the residual ---*/
  Res_Time_Cont = NULL;
  if (dynamic) {
    Res_Time_Cont = new su2double [nVar];
  }
  
  /*--- Matrices to impose clamped boundary conditions ---*/
  
  mZeros_Aux = new su2double *[nDim];
  for(iDim = 0; iDim < nDim; iDim++)
    mZeros_Aux[iDim] = new su2double[nDim];
  
  mId_Aux = new su2double *[nDim];
  for(iDim = 0; iDim < nDim; iDim++)
    mId_Aux[iDim] = new su2double[nDim];
  
  for(iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      mZeros_Aux[iDim][jDim] = 0.0;
      mId_Aux[iDim][jDim] = 0.0;
    }
    mId_Aux[iDim][iDim] = 1.0;
  }
  
  
  /*--- Initialization of matrix structures ---*/
  if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Non-Linear Elasticity)." << endl;
  
  Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);
  
  if (dynamic) {
    MassMatrix.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);
    TimeRes_Aux.Initialize(nPoint, nPointDomain, nVar, 0.0);
    TimeRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  }
  
  /*--- Initialization of linear solver structures ---*/
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  LinSysAux.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  LinSysReact.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Initialize the auxiliary vector and matrix for the computation of the nodal Reactions ---*/
  
  normalVertex = new su2double [nDim];
  
  stressTensor = new su2double* [nDim];
  for (iVar = 0; iVar < nVar; iVar++) {
    stressTensor[iVar] = new su2double [nDim];
  }
  
  /*---- Initialize the auxiliary vector for the solution predictor ---*/
  
  solutionPredictor = new su2double [nVar];
  
  iElem_iDe = NULL;

  /*--- Initialize the value of the total objective function ---*/
   Total_OFRefGeom = 0.0;
   Total_OFRefNode = 0.0;
   Total_OFVolFrac = 0.0;

   /*--- Initialize the value of the global objective function ---*/
   Global_OFRefGeom = 0.0;
   Global_OFRefNode = 0.0;

   /*--- Initialize the value of the total gradient for the forward mode ---*/
   Total_ForwardGradient = 0.0;

   if (config->GetDirectDiff() == D_YOUNG ||
       config->GetDirectDiff() == D_POISSON ||
       config->GetDirectDiff() == D_RHO ||
       config->GetDirectDiff() == D_RHO_DL ||
       config->GetDirectDiff() == D_EFIELD ||
       config->GetDirectDiff() == D_MACH ||
       config->GetDirectDiff() == D_PRESSURE){

     /*--- Header of the temporary output file ---*/
     ofstream myfile_res;

     if (config->GetDirectDiff() == D_YOUNG) myfile_res.open ("Output_Direct_Diff_E.txt");
     if (config->GetDirectDiff() == D_POISSON) myfile_res.open ("Output_Direct_Diff_Nu.txt");
     if (config->GetDirectDiff() == D_RHO) myfile_res.open ("Output_Direct_Diff_Rho.txt");
     if (config->GetDirectDiff() == D_RHO_DL) myfile_res.open ("Output_Direct_Diff_Rho_DL.txt");
     if (config->GetDirectDiff() == D_EFIELD) myfile_res.open ("Output_Direct_Diff_EField.txt");

     if (config->GetDirectDiff() == D_MACH) myfile_res.open ("Output_Direct_Diff_Mach.txt");
     if (config->GetDirectDiff() == D_PRESSURE) myfile_res.open ("Output_Direct_Diff_Pressure.txt");

     myfile_res << "Objective Function " << "\t";

     if (dynamic) myfile_res << "O. Function Averaged " << "\t";

     myfile_res << "Sensitivity Local" << "\t";

     myfile_res << "Sensitivity Averaged" << "\t";

     myfile_res << endl;

     myfile_res.close();

   }

  /*--- Initialize the BGS residuals in FSI problems. ---*/
  if (config->GetMultizone_Residual()){

    FSI_Residual      = 0.0;
    RelaxCoeff        = 1.0;
    ForceCoeff        = 1.0;

    Residual_BGS      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_BGS[iVar]  = 0.0;
    Residual_Max_BGS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max_BGS[iVar]  = 0.0;

    /*--- Define some structures for locating max residuals ---*/

    Point_Max_BGS       = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max_BGS[iVar]  = 0;
    Point_Max_Coord_BGS = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord_BGS[iVar] = new su2double[nDim];
      for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord_BGS[iVar][iDim] = 0.0;
    }
  }
  else{
    ForceCoeff        = 1.0;
  }

  /*--- Penalty value - to maintain constant the stiffness in optimization problems - TODO: this has to be improved ---*/
  PenaltyValue = 0.0;

  /*--- Perform the MPI communication of the solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*--- If dynamic, we also need to communicate the old solution ---*/
  
  if(dynamic) Set_MPI_Solution_Old(geometry, config);
  
}

CFEASolver::~CFEASolver(void) {
  
  unsigned short iVar, jVar;
  unsigned long iElem;
  
  if (element_container != NULL) {
    for (iVar = 0; iVar < MAX_TERMS; iVar++) {
      for (jVar = 0; jVar < MAX_FE_KINDS; jVar++) {
        if (element_container[iVar][jVar] != NULL) delete element_container[iVar][jVar];
      }
      delete [] element_container[iVar];
    }
    delete [] element_container;
  }
  
  if (element_properties != NULL){
    for (iElem = 0; iElem < nElement; iElem++)
      if (element_properties[iElem] != NULL) delete element_properties[iElem];
    delete [] element_properties;
  }
  
  for (iVar = 0; iVar < nVar; iVar++) {
    if (Jacobian_s_ij != NULL) delete [] Jacobian_s_ij[iVar];
    if (Jacobian_c_ij != NULL) delete [] Jacobian_c_ij[iVar];
    if (Jacobian_k_ij != NULL) delete[] Jacobian_k_ij[iVar];
    delete [] mZeros_Aux[iVar];
    delete [] mId_Aux[iVar];
    delete [] stressTensor[iVar];
  }
  
  if (Jacobian_s_ij != NULL) delete [] Jacobian_s_ij;
  if (Jacobian_c_ij != NULL) delete [] Jacobian_c_ij;
  if (Jacobian_k_ij != NULL) delete [] Jacobian_k_ij;
  delete [] Res_Stress_i;
  delete [] Res_Ext_Surf;
  if (Res_Time_Cont != NULL) delete[] Res_Time_Cont;
  if (Res_Dead_Load != NULL) delete[] Res_Dead_Load;
  delete [] SolRest;
  delete [] GradN_X;
  delete [] GradN_x;
  
  delete [] mZeros_Aux;
  delete [] mId_Aux;
  
  delete [] nodeReactions;
  
  delete [] normalVertex;
  delete [] stressTensor;
  
  if (iElem_iDe != NULL) delete [] iElem_iDe;
}

void CFEASolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
  
  
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.
  
  unsigned short nSolVar;
  
  if (dynamic) nSolVar = 3 * nVar;
  else nSolVar = nVar;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nSolVar;     nBufferR_Vector = nVertexR*nSolVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sent ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
        if (dynamic) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Buffer_Send_U[(iVar+nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Vel(iVar);
            Buffer_Send_U[(iVar+2*nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Accel(iVar);
          }
        }
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
        if (dynamic) {
          for (iVar = nVar; iVar < 3*nVar; iVar++)
            Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
        }
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy solution variables. ---*/
        for (iVar = 0; iVar < nSolVar; iVar++)
          SolRest[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Store received values back into the variable. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, SolRest[iVar]);
        
        if (dynamic) {
          
          for (iVar = 0; iVar < nVar; iVar++) {
            node[iPoint]->SetSolution_Vel(iVar, SolRest[iVar+nVar]);
            node[iPoint]->SetSolution_Accel(iVar, SolRest[iVar+2*nVar]);
          }
          
        }
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}

void CFEASolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
  
  
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
  unsigned short nSolVar;
  
  nSolVar = 3 * nVar;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nSolVar;     nBufferR_Vector = nVertexR*nSolVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sent ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++) {
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_time_n(iVar);
          Buffer_Send_U[(iVar+nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Vel_time_n(iVar);
          Buffer_Send_U[(iVar+2*nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Accel_time_n(iVar);
        }
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nSolVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy solution variables. ---*/
        for (iVar = 0; iVar < nSolVar; iVar++)
          SolRest[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Store received values back into the variable. ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          node[iPoint]->SetSolution_time_n(iVar, SolRest[iVar]);
          node[iPoint]->SetSolution_Vel_time_n(iVar, SolRest[iVar+nVar]);
          node[iPoint]->SetSolution_Accel_time_n(iVar, SolRest[iVar+2*nVar]);
        }
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}

void CFEASolver::Set_MPI_Solution_DispOnly(CGeometry *geometry, CConfig *config) {
  

  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;         nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sent ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy solution variables. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          SolRest[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Store received values back into the variable. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, SolRest[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}

void CFEASolver::Set_MPI_Solution_Pred(CGeometry *geometry, CConfig *config) {
  

  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;     nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sent ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Pred(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy predicted solution variables back into the variables. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution_Pred(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}

void CFEASolver::Set_MPI_Solution_Pred_Old(CGeometry *geometry, CConfig *config) {
  
  /*--- We are communicating the solution predicted, current and old, and the old solution ---*/
  /*--- necessary for the Aitken relaxation ---*/
  
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
  /*--- Analogous to the dynamic solution, in this case we need 3 * nVar variables per node ---*/
  unsigned short nSolVar;
  nSolVar = 3 * nVar;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nSolVar;     nBufferR_Vector = nVertexR*nSolVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sent ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++) {
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Old(iVar);
          Buffer_Send_U[(iVar+nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Pred(iVar);
          Buffer_Send_U[(iVar+2*nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Pred_Old(iVar);
        }
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nSolVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Store received values back into the variable. ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          node[iPoint]->SetSolution_Old(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);
          node[iPoint]->SetSolution_Pred(iVar, Buffer_Receive_U[(iVar+nVar)*nVertexR+iVertex]);
          node[iPoint]->SetSolution_Pred_Old(iVar, Buffer_Receive_U[(iVar+2*nVar)*nVertexR+iVertex]);
        }
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}

void CFEASolver::Set_ElementProperties(CGeometry *geometry, CConfig *config) {

  unsigned long iElem;
  unsigned long index;

  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();
  
  bool topology_mode = config->GetTopology_Optimization();

  string filename;
  ifstream properties_file;

  element_properties = new CElementProperty*[nElement];

  /*--- Restart the solution from file information ---*/

  filename = config->GetFEA_FileName();

  /*--- If multizone, append zone name ---*/
  if (nZone > 1)
    filename = config->GetMultizone_FileName(filename, iZone);

  if (rank == MASTER_NODE) cout << "Filename: " << filename << "." << endl;

  properties_file.open(filename.data(), ios::in);

  /*--- In case there is no file, all elements get the same property (0) ---*/

  if (properties_file.fail()) {
    if (rank == MASTER_NODE){
      cout << "There is no element-based properties file." << endl;
      cout << "The structural domain has uniform properties." << endl;
      
      if (topology_mode)
        SU2_MPI::Error("Topology mode requires an element-based properties file.",CURRENT_FUNCTION);
    }

    for (iElem = 0; iElem < nElement; iElem++){
      element_properties[iElem] = new CElementProperty(0, 0, 0, 0);
    }

    element_based = false;

  }
  else{

    element_based = true;

    /*--- In case this is a parallel simulation, we need to perform the
       Global2Local index transformation first. ---*/

    long *Global2Local = new long[geometry->GetGlobal_nElemDomain()];

    /*--- First, set all indices to a negative value by default ---*/

    for (iElem = 0; iElem < geometry->GetGlobal_nElemDomain(); iElem++)
      Global2Local[iElem] = -1;

    /*--- Now fill array with the transform values only for the points in the rank (including halos) ---*/

    for (iElem = 0; iElem < nElement; iElem++)
      Global2Local[geometry->elem[iElem]->GetGlobalIndex()] = iElem;

    /*--- Read all lines in the restart file ---*/

    long iElem_Local;
    unsigned long iElem_Global_Local = 0, iElem_Global = 0; string text_line;
    unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

    /*--- The first line is the header ---*/

    getline (properties_file, text_line);

    for (iElem_Global = 0; iElem_Global < geometry->GetGlobal_nElemDomain(); iElem_Global++ ) {

      getline (properties_file, text_line);

      istringstream point_line(text_line);

      /*--- Retrieve local index. If this element from the restart file lives
         only on a different processor, the value of iPoint_Local will be -1.
         Otherwise, the local index for this node on the current processor
         will be returned and used to instantiate the vars. ---*/

      iElem_Local = Global2Local[iElem_Global];

      if (iElem_Local >= 0) {

        if (config->GetDE_Effects())
          point_line >> index >> elProperties[0] >> elProperties[1] >> elProperties[2] >> elProperties[3];
        else
          point_line >> index >> elProperties[0] >> elProperties[1] >> elProperties[2] >> elProperties[3];

        element_properties[iElem_Local] = new CElementProperty(elProperties[0],
                                                         elProperties[1],
                                                         elProperties[2],
                                                         elProperties[3]);

        /*--- For backwards compatibility we only read a fifth column in topology mode ---*/
        if (topology_mode) {
          su2double elDensity;
          point_line >> elDensity;
          element_properties[iElem_Local]->SetDesignDensity(elDensity);
        }

        iElem_Global_Local++;
      }

    }

    /*--- Detect a wrong solution file ---*/

    if (iElem_Global_Local < nElement) { sbuf_NotMatching = 1; }

#ifndef HAVE_MPI
    rbuf_NotMatching = sbuf_NotMatching;
#else
    SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (rbuf_NotMatching != 0) {
      SU2_MPI::Error(string("The properties file ") + filename + string(" doesn't match with the mesh file!\n")  + 
                     string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
    }

    /*--- Close the restart file ---*/

    properties_file.close();

    /*--- Free memory needed for the transformation ---*/

    delete [] Global2Local;


  }


}


void CFEASolver::Set_Prestretch(CGeometry *geometry, CConfig *config) {
  
  unsigned long iPoint;
  unsigned long index;
  
  unsigned short iVar;
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();
  
  string filename;
  ifstream prestretch_file;
  
  
  /*--- Restart the solution from file information ---*/
  
  filename = config->GetPrestretch_FEMFileName();
  
  /*--- If multizone, append zone name ---*/
  if (nZone > 1)
    filename = config->GetMultizone_FileName(filename, iZone);
  
  if (rank == MASTER_NODE) cout << "Filename: " << filename << "." << endl;
  
  prestretch_file.open(filename.data(), ios::in);
  
  /*--- In case there is no file ---*/
  
  if (prestretch_file.fail()) {
    SU2_MPI::Error(string("There is no FEM prestretch reference file ") + filename, CURRENT_FUNCTION);
  }
  /*--- In case this is a parallel simulation, we need to perform the
   Global2Local index transformation first. ---*/
  
  map<unsigned long,unsigned long> Global2Local;
  map<unsigned long,unsigned long>::const_iterator MI;
  
  /*--- Now fill array with the transform values only for local points ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
  
  /*--- Read all lines in the restart file ---*/
  
  long iPoint_Local;
  unsigned long iPoint_Global_Local = 0, iPoint_Global = 0; string text_line;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;
  
  /*--- The first line is the header ---*/
  
  getline (prestretch_file, text_line);
  
  while (getline (prestretch_file, text_line)) {
    istringstream point_line(text_line);
    
    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/
    
    MI = Global2Local.find(iPoint_Global);
    if (MI != Global2Local.end()) {
    
      iPoint_Local = Global2Local[iPoint_Global];
      
      if (nDim == 2) point_line >> Solution[0] >> Solution[1] >> index;
      if (nDim == 3) point_line >> Solution[0] >> Solution[1] >> Solution[2] >> index;
      
      for (iVar = 0; iVar < nVar; iVar++) node[iPoint_Local]->SetPrestretch(iVar, Solution[iVar]);
      
      iPoint_Global_Local++;
    }
    iPoint_Global++;
  }
  
  /*--- Detect a wrong solution file ---*/
  
  if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }
  
#ifndef HAVE_MPI
  rbuf_NotMatching = sbuf_NotMatching;
#else
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (rbuf_NotMatching != 0) {
      SU2_MPI::Error(string("The solution file ") + filename + string(" doesn't match with the mesh file!\n") +
                     string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }
  
  /*--- Close the restart file ---*/
  
  prestretch_file.close();
  
  /*--- We need to communicate here the prestretched geometry for the halo nodes. ---*/
  /*--- We avoid creating a new function as this may be reformatted.              ---*/

  unsigned short iMarker, MarkerS, MarkerR;
  unsigned long iVertex, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;

#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];

      /*--- Copy the solution that should be sent ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetPrestretch(iVar);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }

#endif

      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;

      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();

        /*--- Copy solution variables. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          SolRest[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];

        /*--- Store received values back into the variable. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetPrestretch(iVar, SolRest[iVar]);

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;

    }

  }


}

void CFEASolver::Set_ReferenceGeometry(CGeometry *geometry, CConfig *config) {

  unsigned long iPoint;
  unsigned long index;

  unsigned short iVar;
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();
  unsigned short file_format = config->GetRefGeom_FileFormat();

  string filename;
  su2double dull_val;
  ifstream reference_file;


  /*--- Restart the solution from file information ---*/

  filename = config->GetRefGeom_FEMFileName();

  /*--- If multizone, append zone name ---*/
  if (nZone > 1)
    filename = config->GetMultizone_FileName(filename, iZone);

  reference_file.open(filename.data(), ios::in);

  /*--- In case there is no file ---*/

  if (reference_file.fail()) {
    SU2_MPI::Error( "There is no FEM reference geometry file!!", CURRENT_FUNCTION);
  }

  if (rank == MASTER_NODE) cout << "Filename: " << filename << " and format " << file_format << "." << endl;

  /*--- In case this is a parallel simulation, we need to perform the
   Global2Local index transformation first. ---*/

  long *Global2Local = new long[geometry->GetGlobal_nPointDomain()];

  /*--- First, set all indices to a negative value by default ---*/

  for (iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++)
    Global2Local[iPoint] = -1;

  /*--- Now fill array with the transform values only for local points ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;

  /*--- Read all lines in the restart file ---*/

  long iPoint_Local;
  unsigned long iPoint_Global_Local = 0, iPoint_Global = 0; string text_line;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

  /*--- The first line is the header ---*/

  getline (reference_file, text_line);

  while (getline (reference_file, text_line)) {
    istringstream point_line(text_line);

    /*--- Retrieve local index. If this node from the restart file lives
       on a different processor, the value of iPoint_Local will be -1.
       Otherwise, the local index for this node on the current processor
       will be returned and used to instantiate the vars. ---*/

    iPoint_Local = Global2Local[iPoint_Global];

    if (iPoint_Local >= 0) {

      if (nDim == 2) point_line >> index >> dull_val >> dull_val >> Solution[0] >> Solution[1];
      if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2];

      for (iVar = 0; iVar < nVar; iVar++) node[iPoint_Local]->SetReference_Geometry(iVar, Solution[iVar]);

      iPoint_Global_Local++;
    }
    iPoint_Global++;
  }

  /*--- Detect a wrong solution file ---*/

  if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }

#ifndef HAVE_MPI
  rbuf_NotMatching = sbuf_NotMatching;
#else
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
  
  if (rbuf_NotMatching != 0) {
    SU2_MPI::Error(string("The solution file ") + filename + string(" doesn't match with the mesh file!\n")  + 
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- I don't think we need to communicate ---*/

  /*--- Close the restart file ---*/

  reference_file.close();

  /*--- Free memory needed for the transformation ---*/

  delete [] Global2Local;

}



void CFEASolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics, unsigned short iMesh, unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output) {
  
  
  unsigned long iPoint;
  bool initial_calc = (config->GetExtIter() == 0);                  // Checks if it is the first calculation.
  bool first_iter = (config->GetIntIter() == 0);                          // Checks if it is the first iteration
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.
  bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
  bool newton_raphson = (config->GetKind_SpaceIteScheme_FEA() == NEWTON_RAPHSON);    // Newton-Raphson method
  bool restart = config->GetRestart();                        // Restart analysis
  bool initial_calc_restart = (SU2_TYPE::Int(config->GetExtIter()) == config->GetDyn_RestartIter()); // Initial calculation for restart
  
  bool disc_adj_fem = (config->GetKind_Solver() == DISC_ADJ_FEM);     // Discrete adjoint FEM solver
  
  bool incremental_load = config->GetIncrementalLoad();               // If an incremental load is applied
  
  bool body_forces = config->GetDeadLoad();                     // Body forces (dead loads).
  
  bool fsi = config->GetFSI_Simulation();
  bool consistent_interpolation = (!config->GetConservativeInterpolation() ||
                                  (config->GetKindInterpolation() == WEIGHTED_AVERAGE));
  
  bool topology_mode = config->GetTopology_Optimization();  // Density-based topology optimization
  
  
  /*--- Set vector entries to zero ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
    LinSysAux.SetBlock_Zero(iPoint);
    LinSysRes.SetBlock_Zero(iPoint);
    LinSysSol.SetBlock_Zero(iPoint);
  }
  
  /*--- Set matrix entries to zero ---*/
  
  /*
   * If the problem is linear, we only need one Jacobian matrix in the problem, because
   * it is going to be constant along the calculations. Therefore, we only initialize
   * the Jacobian matrix once, at the beginning of the simulation.
   *
   * We don't need first_iter, because there is only one iteration per time step in linear analysis.
   * For dynamic problems we also need to recompute the Jacobian as that is where the RHS is computed
   * as a residual (and we need that for AD).
   */
  if ((initial_calc && linear_analysis)||
      (restart && initial_calc_restart && linear_analysis) ||
      (dynamic && disc_adj_fem) ||
      (dynamic && linear_analysis)) {
    Jacobian.SetValZero();
  }
  
  /*
   * For topology optimization we apply a filter on the design density field to avoid
   * numerical issues (checkerboards), ensure mesh independence, and impose a length scale.
   * This has to be done before computing the mass matrix and the dead load terms.
   * This filter, and the volume fraction objective function, require the element volumes,
   * so we ask "geometry" to compute them.
   * This only needs to be done for the undeformed (initial) shape.
   */
  if (topology_mode && (restart || initial_calc || initial_calc_restart || disc_adj_fem)) {
    geometry->SetElemVolume(config);
    FilterElementDensities(geometry,config);
  }
  
  /*
   * If the problem is dynamic, we need a mass matrix, which will be constant along the calculation
   * both for linear and nonlinear analysis. Only initialized once, at the first time step.
   *
   * The same with the integration constants, as for now we consider the time step to be constant.
   *
   * We need first_iter, because in nonlinear problems there are more than one subiterations in the first time step.
   */
  if ((dynamic && initial_calc && first_iter) ||
      (dynamic && restart && initial_calc_restart && first_iter) ||
      (dynamic && disc_adj_fem)) {
    MassMatrix.SetValZero();
    Compute_IntegrationConstants(config);
    Compute_MassMatrix(geometry, solver_container, numerics, config);
  }
  
  /*
   * If body forces are taken into account, we need to compute the term that goes into the residual,
   * which will be constant along the calculation both for linear and nonlinear analysis.
   *
   * Only initialized once, at the first iteration or the beginning of the calculation after a restart.
   *
   * We need first_iter, because in nonlinear problems there are more than one subiterations in the first time step.
   */
  
  if ((body_forces && initial_calc && first_iter) ||
      (body_forces && restart && initial_calc_restart && first_iter)) {
    // If the load is incremental, we have to reset the variable to avoid adding up over the increments
    if (incremental_load) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) node[iPoint]->Clear_BodyForces_Res();
    }
    // Compute the dead load term
    Compute_DeadLoad(geometry, solver_container, numerics, config);
  }
  
  /*
   * If the problem is nonlinear, we need to initialize the Jacobian and the stiffness matrix at least at the beginning
   * of each time step. If the solution method is Newton Rapshon, we initialize it also at the beginning of each
   * iteration.
   */
  
  if ((nonlinear_analysis) && ((newton_raphson) || (first_iter)))  {
    Jacobian.SetValZero();
    //    StiffMatrix.SetValZero();
  }
  
  /*
   * Some external forces may be considered constant over the time step.
   */
  if (first_iter)  {
    for (iPoint = 0; iPoint < nPoint; iPoint++) node[iPoint]->Clear_SurfaceLoad_Res();
  }
  
  /*
   * If we apply nonlinear forces, we need to clear the residual on each iteration
   */
  unsigned short iMarker;
  unsigned long iVertex;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
    switch (config->GetMarker_All_KindBC(iMarker)) {
      case LOAD_BOUNDARY:
        /*--- Only if the load is nonzero - reduces computational cost ---*/
        if(config->GetLoad_Value(config->GetMarker_All_TagBound(iMarker)) != 0 ) {
          /*--- For all the vertices in the marker iMarker ---*/
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            /*--- Retrieve the point ID ---*/
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            /*--- Clear the residual of the node, to avoid adding on previous values ---*/
            node[iPoint]->Clear_SurfaceLoad_Res();
          }
        }
        break;
      case DAMPER_BOUNDARY:
        /*--- For all the vertices in the marker iMarker ---*/
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          /*--- Retrieve the point ID ---*/
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          /*--- Clear the residual of the node, to avoid adding on previous values ---*/
          node[iPoint]->Clear_SurfaceLoad_Res();
        }
        break;
    }
  }
  
  /*
   * FSI loads (computed upstream) need to be integrated if a nonconservative interpolation scheme is in use
   */
  if (fsi && first_iter && consistent_interpolation) Integrate_FSI_Loads(geometry,config);

}

void CFEASolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) { }

void CFEASolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {
  
  unsigned long iPoint, nPoint;
  bool incremental_load = config->GetIncrementalLoad();              // If an incremental load is applied
  
  nPoint = geometry[MESH_0]->GetnPoint();
  
  /*--- We store the current solution as "Solution Old", for the case that we need to retrieve it ---*/
  
  if (incremental_load) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) node[iPoint]->Set_OldSolution();
  }
  
  
}

void CFEASolver::ResetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {
  
  unsigned long iPoint, nPoint;
  bool incremental_load = config->GetIncrementalLoad();              // If an incremental load is applied
  
  nPoint = geometry[MESH_0]->GetnPoint();
  
  /*--- We store the current solution as "Solution Old", for the case that we need to retrieve it ---*/
  
  if (incremental_load) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) node[iPoint]->Set_Solution();
  }
  
}

void CFEASolver::Compute_StiffMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {
  
  unsigned long iElem, iVar, jVar;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord, val_Sol;
  int EL_KIND = 0;
  
  su2double *Kab = NULL, *Ta  = NULL;
  unsigned short NelNodes, jNode;
  
  bool topology_mode = config->GetTopology_Optimization();
  su2double simp_exponent = config->GetSIMP_Exponent();
  su2double simp_minstiff = config->GetSIMP_MinStiffness();
  
  /*--- Loops over all the elements ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)      {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {nNodes = 4; EL_KIND = EL_QUAD;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)   {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)       {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)         {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)    {nNodes = 8; EL_KIND = EL_HEXA;}
    
    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      
      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      
      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
        element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
        element_container[FEA_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
      }
    }
    
    /*--- In topology mode determine the penalty to apply to the stiffness ---*/
    su2double simp_penalty = 1.0;
    if (topology_mode) {
      simp_penalty = simp_minstiff+(1.0-simp_minstiff)*pow(element_properties[iElem]->GetPhysicalDensity(),simp_exponent);
    }
    
    /*--- Set the properties of the element ---*/
    element_container[FEA_TERM][EL_KIND]->Set_ElProperties(element_properties[iElem]);
    
    /*--- Compute the components of the jacobian and the stress term ---*/
    if (element_based){
      numerics[element_properties[iElem]->GetMat_Mod()]->Compute_Tangent_Matrix(element_container[FEA_TERM][EL_KIND], config);
    }
    else{
      numerics[FEA_TERM]->Compute_Tangent_Matrix(element_container[FEA_TERM][EL_KIND], config);
    }
    
    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();
    
    for (iNode = 0; iNode < NelNodes; iNode++) {
      
      Ta = element_container[FEA_TERM][EL_KIND]->Get_Kt_a(iNode);
      for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = simp_penalty*Ta[iVar];
      
      LinSysRes.SubtractBlock(indexNode[iNode], Res_Stress_i);
      
      for (jNode = 0; jNode < NelNodes; jNode++) {
        
        Kab = element_container[FEA_TERM][EL_KIND]->Get_Kab(iNode, jNode);
        
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) {
            Jacobian_ij[iVar][jVar] = simp_penalty*Kab[iVar*nVar+jVar];
          }
        }
    
        Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_ij);
      }
      
    }
    
  }
  
  
}

void CFEASolver::Compute_StiffMatrix_NodalStressRes(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {
  
  unsigned long iElem, iVar, jVar;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord, val_Sol, val_Ref = 0.0;
  int EL_KIND = 0;
  
  bool prestretch_fem = config->GetPrestretch();
  
  su2double Ks_ab;
  su2double *Kab = NULL;
  su2double *Kk_ab = NULL;
  su2double *Ta = NULL;
  
  su2double *Ta_DE = NULL;
  su2double Ks_ab_DE = 0.0;
  
  unsigned short NelNodes, jNode;
  
  bool incompressible = (config->GetMaterialCompressibility() == INCOMPRESSIBLE_MAT);
  bool de_effects = config->GetDE_Effects();

  bool topology_mode = config->GetTopology_Optimization();
  su2double simp_exponent = config->GetSIMP_Exponent();
  su2double simp_minstiff = config->GetSIMP_MinStiffness();
  
  /*--- Loops over all the elements ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)      {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {nNodes = 4; EL_KIND = EL_QUAD;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)   {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)       {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)         {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)    {nNodes = 8; EL_KIND = EL_HEXA;}
    
    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;

        /*--- Set current coordinate ---*/
        element_container[FEA_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
        if (de_effects) element_container[DE_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
        if (incompressible) element_container[INC_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);

        /*--- Set reference coordinate ---*/
        if (prestretch_fem) {
          val_Ref = node[indexNode[iNode]]->GetPrestretch(iDim);
          element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Ref, iNode, iDim);
          if (de_effects) element_container[DE_TERM][EL_KIND]->SetRef_Coord(val_Ref, iNode, iDim);
          if (incompressible) element_container[INC_TERM][EL_KIND]->SetRef_Coord(val_Ref, iNode, iDim);
        }
        else {
          element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
          if (de_effects) element_container[DE_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
          if (incompressible) element_container[INC_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
        }
      }
    }
    
    /*--- In topology mode determine the penalty to apply to the stiffness ---*/
    su2double simp_penalty = 1.0;
    if (topology_mode) {
      simp_penalty = simp_minstiff+(1.0-simp_minstiff)*pow(element_properties[iElem]->GetPhysicalDensity(),simp_exponent);
    }
    
    /*--- Set the properties of the element ---*/
    element_container[FEA_TERM][EL_KIND]->Set_ElProperties(element_properties[iElem]);
    if (de_effects) element_container[DE_TERM][EL_KIND]->Set_ElProperties(element_properties[iElem]);
    if (incompressible) element_container[INC_TERM][EL_KIND]->Set_ElProperties(element_properties[iElem]);
    
    /*--- If incompressible, we compute the Mean Dilatation term first so the volume is already computed ---*/
    if (incompressible) numerics[FEA_TERM]->Compute_MeanDilatation_Term(element_container[INC_TERM][EL_KIND], config);
    
    /*--- Compute the components of the Jacobian and the stress term for the material ---*/
    if (element_based){
      numerics[element_properties[iElem]->GetMat_Mod()]->Compute_Tangent_Matrix(element_container[FEA_TERM][EL_KIND], config);
    }
    else{
      numerics[FEA_TERM]->Compute_Tangent_Matrix(element_container[FEA_TERM][EL_KIND], config);
    }
    
    /*--- Compute the electric component of the Jacobian and the stress term ---*/
    if (de_effects) numerics[DE_TERM]->Compute_Tangent_Matrix(element_container[DE_TERM][EL_KIND], config);
    
    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();
    
    for (iNode = 0; iNode < NelNodes; iNode++) {
      
      Ta = element_container[FEA_TERM][EL_KIND]->Get_Kt_a(iNode);
      for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = simp_penalty*Ta[iVar];
      
      /*--- Check if this is my node or not ---*/
      LinSysRes.SubtractBlock(indexNode[iNode], Res_Stress_i);
      
      /*--- Retrieve the electric contribution to the Residual ---*/
      if (de_effects){
        Ta_DE = element_container[DE_TERM][EL_KIND]->Get_Kt_a(iNode);
        for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = simp_penalty*Ta_DE[iVar];
        LinSysRes.SubtractBlock(indexNode[iNode], Res_Stress_i);
        
      }
      
      for (jNode = 0; jNode < NelNodes; jNode++) {
        
        /*--- Retrieve the values of the FEA term ---*/
        Kab = element_container[FEA_TERM][EL_KIND]->Get_Kab(iNode, jNode);
        Ks_ab = element_container[FEA_TERM][EL_KIND]->Get_Ks_ab(iNode,jNode);
        if (incompressible) Kk_ab = element_container[INC_TERM][EL_KIND]->Get_Kk_ab(iNode,jNode);
        
        for (iVar = 0; iVar < nVar; iVar++) {
          Jacobian_s_ij[iVar][iVar] = simp_penalty*Ks_ab;
          for (jVar = 0; jVar < nVar; jVar++) {
            Jacobian_c_ij[iVar][jVar] = simp_penalty*Kab[iVar*nVar+jVar];
            if (incompressible) Jacobian_k_ij[iVar][jVar] = simp_penalty*Kk_ab[iVar*nVar+jVar];
          }
        }
        
        Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_c_ij);
        Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_s_ij);
        if (incompressible) Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_k_ij);
        
        /*--- Retrieve the electric contribution to the Jacobian ---*/
        if (de_effects){
          //          Kab_DE = element_container[DE_TERM][EL_KIND]->Get_Kab(iNode, jNode);
          Ks_ab_DE = element_container[DE_TERM][EL_KIND]->Get_Ks_ab(iNode,jNode);
          
          for (iVar = 0; iVar < nVar; iVar++){
            Jacobian_s_ij[iVar][iVar] = simp_penalty*Ks_ab_DE;
          }
          
          Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_s_ij);
        }
        
      }
      
    }
    
  }
  
}

void CFEASolver::Compute_MassMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {
  
  unsigned long iElem, iVar;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord;
  int EL_KIND = 0;
  
  su2double Mab;
  unsigned short NelNodes, jNode;
  
  bool topology_mode = config->GetTopology_Optimization();
  su2double simp_minstiff = config->GetSIMP_MinStiffness();
  
  /*--- Loops over all the elements ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL)    {nNodes = 4; EL_KIND = EL_QUAD;}
    
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}
    
    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
      }
    }
    
    /*--- In topology mode determine the penalty to apply to the mass, linear function of the physical density ---*/
    su2double simp_penalty = 1.0;
    if (topology_mode) {
      simp_penalty = simp_minstiff+(1.0-simp_minstiff)*element_properties[iElem]->GetPhysicalDensity();
    }
    
    /*--- Set the properties of the element ---*/
    element_container[FEA_TERM][EL_KIND]->Set_ElProperties(element_properties[iElem]);
    
    numerics[FEA_TERM]->Compute_Mass_Matrix(element_container[FEA_TERM][EL_KIND], config);
    
    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();

    for (iNode = 0; iNode < NelNodes; iNode++) {
      
      for (jNode = 0; jNode < NelNodes; jNode++) {
        
        Mab = element_container[FEA_TERM][EL_KIND]->Get_Mab(iNode, jNode);
        
        for (iVar = 0; iVar < nVar; iVar++) {
          MassMatrix_ij[iVar][iVar] = simp_penalty*Mab;
        }
        
        MassMatrix.AddBlock(indexNode[iNode], indexNode[jNode], MassMatrix_ij);
        
      }
      
    }
    
  }
  
}

void CFEASolver::Compute_MassRes(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {

  unsigned long iElem, iVar, iPoint;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord;
  int EL_KIND = 0;

  su2double Mab;
  unsigned short NelNodes, jNode;
  
  bool topology_mode = config->GetTopology_Optimization();
  su2double simp_minstiff = config->GetSIMP_MinStiffness();

  /*--- Set vector entries to zero ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
    TimeRes.SetBlock_Zero(iPoint);
  }

  /*--- Loops over all the elements ---*/

  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL)    {nNodes = 4; EL_KIND = EL_QUAD;}

    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}

    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/

    for (iNode = 0; iNode < nNodes; iNode++) {
      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
      }
    }

    /*--- In topology mode determine the penalty to apply to the mass, linear function of the physical density ---*/
    su2double simp_penalty = 1.0;
    if (topology_mode) {
      simp_penalty = simp_minstiff+(1.0-simp_minstiff)*element_properties[iElem]->GetPhysicalDensity();
    }

    /*--- Set the properties of the element ---*/
    element_container[FEA_TERM][EL_KIND]->Set_ElProperties(element_properties[iElem]);

    numerics[FEA_TERM]->Compute_Mass_Matrix(element_container[FEA_TERM][EL_KIND], config);

    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();

    for (iNode = 0; iNode < NelNodes; iNode++) {

      for (jNode = 0; jNode < NelNodes; jNode++) {

        Mab = element_container[FEA_TERM][EL_KIND]->Get_Mab(iNode, jNode);

        for (iVar = 0; iVar < nVar; iVar++) {
          Residual_i[iVar] = simp_penalty * Mab * TimeRes_Aux.GetBlock(indexNode[iNode],iVar);
          Residual_j[iVar] = simp_penalty * Mab * TimeRes_Aux.GetBlock(indexNode[jNode],iVar);
        }

        TimeRes.AddBlock(indexNode[iNode],Residual_i);
        TimeRes.AddBlock(indexNode[jNode],Residual_j);

      }

    }

  }

}

void CFEASolver::Compute_NodalStressRes(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {
  
  
  unsigned long iElem, iVar;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord, val_Sol, val_Ref = 0.0;
  int EL_KIND = 0;
  
  bool prestretch_fem = config->GetPrestretch();
  
  su2double *Ta = NULL;
  unsigned short NelNodes;
  
  bool topology_mode = config->GetTopology_Optimization();
  su2double simp_exponent = config->GetSIMP_Exponent();
  su2double simp_minstiff = config->GetSIMP_MinStiffness();
  
  /*--- Loops over all the elements ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {nNodes = 4; EL_KIND = EL_QUAD;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}
    
    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
        element_container[FEA_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
        if (prestretch_fem) {
          val_Ref = node[indexNode[iNode]]->GetPrestretch(iDim);
          element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Ref, iNode, iDim);
        }
        else {
          element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
        }
      }
    }
    
    /*--- In topology mode determine the penalty to apply to the stiffness ---*/
    su2double simp_penalty = 1.0;
    if (topology_mode) {
      simp_penalty = simp_minstiff+(1.0-simp_minstiff)*pow(element_properties[iElem]->GetPhysicalDensity(),simp_exponent);
    }
    
    /*--- Set the properties of the element ---*/
    element_container[FEA_TERM][EL_KIND]->Set_ElProperties(element_properties[iElem]);
    
    /*--- Compute the components of the jacobian and the stress term ---*/
    if (element_based){
      numerics[element_properties[iElem]->GetMat_Mod()]->Compute_NodalStress_Term(element_container[FEA_TERM][EL_KIND], config);
    }
    else{
      numerics[FEA_TERM]->Compute_NodalStress_Term(element_container[FEA_TERM][EL_KIND], config);
    }
    
    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();
    
    for (iNode = 0; iNode < NelNodes; iNode++) {
      
      Ta = element_container[FEA_TERM][EL_KIND]->Get_Kt_a(iNode);
      for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = simp_penalty*Ta[iVar];
      
      LinSysRes.SubtractBlock(indexNode[iNode], Res_Stress_i);
      
    }
    
  }

}

void CFEASolver::Compute_NodalStress(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {
  
  unsigned long iPoint, iElem, iVar;
  unsigned short iNode, iDim, iStress;
  unsigned short nNodes = 0, nStress;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord, val_Sol, val_Ref = 0.0;
  int EL_KIND = 0;
  
  bool prestretch_fem = config->GetPrestretch();
  
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);
  
  bool topology_mode = config->GetTopology_Optimization();
  su2double simp_exponent = config->GetSIMP_Exponent();
  
  if (nDim == 2) nStress = 3;
  else nStress = 6;
  
  su2double *Ta = NULL;
  
  unsigned short NelNodes;
  
  /*--- Restart stress to avoid adding results from previous time steps ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iStress = 0; iStress < nStress; iStress++) {
      node[iPoint]->SetStress_FEM(iStress, 0.0);
    }
  }
  
  /*--- Loops over all the elements ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {nNodes = 4; EL_KIND = EL_QUAD;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}
    
    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      //      for (iDim = 0; iDim < nDim; iDim++) {
      //        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
      //        val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
      //        element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
      //        element_container[FEA_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
      //      }
      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
        element_container[FEA_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
        if (prestretch_fem) {
          val_Ref = node[indexNode[iNode]]->GetPrestretch(iDim);
          element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Ref, iNode, iDim);
        }
        else {
          element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
        }
      }
    }
    
    /*--- Correct the stresses and reactions for topology optimization densities. ---*/
    su2double simp_penalty = 1.0;
    if (topology_mode) {
      simp_penalty = pow(element_properties[iElem]->GetPhysicalDensity(),simp_exponent);
    }
    
    /*--- Set the properties of the element ---*/
    element_container[FEA_TERM][EL_KIND]->Set_ElProperties(element_properties[iElem]);
    
    /*--- Compute the components of the jacobian and the stress term ---*/
    if (element_based){
      numerics[element_properties[iElem]->GetMat_Mod()]->Compute_Averaged_NodalStress(element_container[FEA_TERM][EL_KIND], config);
    }
    else{
      numerics[FEA_TERM]->Compute_Averaged_NodalStress(element_container[FEA_TERM][EL_KIND], config);
    }
    
    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();
    
    for (iNode = 0; iNode < NelNodes; iNode++) {
      
      /*--- This only works if the problem is nonlinear ---*/
      Ta = element_container[FEA_TERM][EL_KIND]->Get_Kt_a(iNode);
      for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = simp_penalty*Ta[iVar];
      
      LinSysReact.AddBlock(indexNode[iNode], Res_Stress_i);
      
      for (iStress = 0; iStress < nStress; iStress++) {
        node[indexNode[iNode]]->AddStress_FEM(iStress, simp_penalty *
                                              (element_container[FEA_TERM][EL_KIND]->Get_NodalStress(iNode, iStress) /
                                               geometry->node[indexNode[iNode]]->GetnElem()) );
      }
      
    }
    
  }
  
  su2double *Stress;
  su2double VonMises_Stress, MaxVonMises_Stress = 0.0;
  su2double Sxx,Syy,Szz,Sxy,Sxz,Syz,S1,S2;
  
  /*--- For the number of nodes in the mesh ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Get the stresses, added up from all the elements that connect to the node ---*/
    
    Stress  = node[iPoint]->GetStress_FEM();
    
    /*--- Compute the stress averaged from all the elements connecting to the node and the Von Mises stress ---*/
    
    if (nDim == 2) {
      
      Sxx=Stress[0];
      Syy=Stress[1];
      Sxy=Stress[2];
      
      S1=(Sxx+Syy)/2+sqrt(((Sxx-Syy)/2)*((Sxx-Syy)/2)+Sxy*Sxy);
      S2=(Sxx+Syy)/2-sqrt(((Sxx-Syy)/2)*((Sxx-Syy)/2)+Sxy*Sxy);
      
      VonMises_Stress = sqrt(S1*S1+S2*S2-2*S1*S2);
      
    }
    else {
      
      Sxx = Stress[0];
      Syy = Stress[1];
      Szz = Stress[3];
      
      Sxy = Stress[2];
      Sxz = Stress[4];
      Syz = Stress[5];
      
      VonMises_Stress = sqrt(0.5*(   pow(Sxx - Syy, 2.0)
                                  + pow(Syy - Szz, 2.0)
                                  + pow(Szz - Sxx, 2.0)
                                  + 6.0*(Sxy*Sxy+Sxz*Sxz+Syz*Syz)
                                  ));
      
    }
    
    node[iPoint]->SetVonMises_Stress(VonMises_Stress);
    
    /*--- Compute the maximum value of the Von Mises Stress ---*/
    
    MaxVonMises_Stress = max(MaxVonMises_Stress, VonMises_Stress);
    
  }
  
#ifdef HAVE_MPI
  
  /*--- Compute MaxVonMises_Stress using all the nodes ---*/
  
  su2double MyMaxVonMises_Stress = MaxVonMises_Stress; MaxVonMises_Stress = 0.0;
  SU2_MPI::Allreduce(&MyMaxVonMises_Stress, &MaxVonMises_Stress, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  
#endif
  
  /*--- Set the value of the MaxVonMises_Stress as the CFEA coeffient ---*/
  
  Total_CFEA = MaxVonMises_Stress;
  
  
  bool outputReactions = false;
  
  if (outputReactions) {
    
    ofstream myfile;
    myfile.open ("Reactions.txt");
    
    unsigned short iMarker;
    unsigned long iVertex;
    su2double val_Reaction;
    
    if (!dynamic) {
      /*--- Loop over all the markers  ---*/
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
        switch (config->GetMarker_All_KindBC(iMarker)) {
            
            /*--- If it corresponds to a clamped boundary  ---*/
            
          case CLAMPED_BOUNDARY:
            
            myfile << "MARKER " << iMarker << ":" << endl;
            
            /*--- Loop over all the vertices  ---*/
            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
              
              /*--- Get node index ---*/
              iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
              
              myfile << "Node " << iPoint << "." << " \t ";
              
              for (iDim = 0; iDim < nDim; iDim++) {
                /*--- Retrieve coordinate ---*/
                val_Coord = geometry->node[iPoint]->GetCoord(iDim);
                myfile << "X" << iDim + 1 << ": " << val_Coord << " \t " ;
              }
              
              for (iVar = 0; iVar < nVar; iVar++) {
                /*--- Retrieve reaction ---*/
                val_Reaction = LinSysReact.GetBlock(iPoint, iVar);
                myfile << "F" << iVar + 1 << ": " << val_Reaction << " \t " ;
              }
              
              myfile << endl;
            }
            myfile << endl;
            break;
        }
    }
    else if (dynamic) {
      
      switch (config->GetKind_TimeIntScheme_FEA()) {
        case (CD_EXPLICIT):
          cout << "NOT IMPLEMENTED YET" << endl;
          break;
        case (NEWMARK_IMPLICIT):
          
          /*--- Loop over all points, and set aux vector TimeRes_Aux = a0*U+a2*U'+a3*U'' ---*/
          for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
            for (iVar = 0; iVar < nVar; iVar++) {
              Residual[iVar] =   a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)        //a0*U(t)
              - a_dt[0]*node[iPoint]->GetSolution(iVar)           //a0*U(t+dt)(k-1)
              + a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)    //a2*U'(t)
              + a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);  //a3*U''(t)
            }
            TimeRes_Aux.SetBlock(iPoint, Residual);
          }

          /*--- Once computed, compute M*TimeRes_Aux ---*/
          MassMatrix.MatrixVectorProduct(TimeRes_Aux,TimeRes,geometry,config);
          
          /*--- Loop over all the markers  ---*/
          for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
            switch (config->GetMarker_All_KindBC(iMarker)) {
                
                /*--- If it corresponds to a clamped boundary  ---*/
                
              case CLAMPED_BOUNDARY:
                
                myfile << "MARKER " << iMarker << ":" << endl;
                
                /*--- Loop over all the vertices  ---*/
                for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                  
                  /*--- Get node index ---*/
                  iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                  
                  myfile << "Node " << iPoint << "." << " \t ";
                  
                  for (iDim = 0; iDim < nDim; iDim++) {
                    /*--- Retrieve coordinate ---*/
                    val_Coord = geometry->node[iPoint]->GetCoord(iDim);
                    myfile << "X" << iDim + 1 << ": " << val_Coord << " \t " ;
                  }
                  
                  for (iVar = 0; iVar < nVar; iVar++) {
                    /*--- Retrieve the time contribution ---*/
                    Res_Time_Cont[iVar] = TimeRes.GetBlock(iPoint, iVar);
                    /*--- Retrieve reaction ---*/
                    val_Reaction = LinSysReact.GetBlock(iPoint, iVar) + Res_Time_Cont[iVar];
                    myfile << "F" << iVar + 1 << ": " << val_Reaction << " \t " ;
                  }
                  
                  myfile << endl;
                }
                myfile << endl;
                break;
            }
          
          
          break;
        case (GENERALIZED_ALPHA):
          cout << "NOT IMPLEMENTED YET" << endl;
          break;
      }
      
    }
    
    
    
    myfile.close();
    
  }
  
}

void CFEASolver::Compute_DeadLoad(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {
  
  unsigned long iElem, iVar;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord;
  int EL_KIND = 0;
  
  su2double *Dead_Load = NULL;
  unsigned short NelNodes;
  
  /*--- Loops over all the elements ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {nNodes = 4; EL_KIND = EL_QUAD;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}
    
    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
      }
    }
    
    /*--- Penalize the dead load, do it by default to avoid unecessary "ifs", since it
          goes to the RHS there is no need to have a minimum value for stability ---*/
    su2double simp_penalty = element_properties[iElem]->GetPhysicalDensity();
    
    /*--- Set the properties of the element ---*/
    element_container[FEA_TERM][EL_KIND]->Set_ElProperties(element_properties[iElem]);
    
    numerics[FEA_TERM]->Compute_Dead_Load(element_container[FEA_TERM][EL_KIND], config);
    
    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();
    
    for (iNode = 0; iNode < NelNodes; iNode++) {
      
      Dead_Load = element_container[FEA_TERM][EL_KIND]->Get_FDL_a(iNode);
      for (iVar = 0; iVar < nVar; iVar++) Res_Dead_Load[iVar] = simp_penalty*Dead_Load[iVar];
      
      node[indexNode[iNode]]->Add_BodyForces_Res(Res_Dead_Load);
      
    }
    
  }
  
  
}

void CFEASolver::Initialize_SystemMatrix(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
}

void CFEASolver::Compute_IntegrationConstants(CConfig *config) {
  
  su2double Delta_t= config->GetDelta_DynTime();
  
  su2double gamma = config->GetNewmark_gamma(), beta = config->GetNewmark_beta();
  
  switch (config->GetKind_TimeIntScheme_FEA()) {
    case (CD_EXPLICIT):
      cout << "NOT IMPLEMENTED YET" << endl;
      break;
    case (NEWMARK_IMPLICIT):
      
      /*--- Integration constants for Newmark scheme ---*/
      
      a_dt[0]= 1 / (beta*pow(Delta_t,2.0));
      a_dt[1]= gamma / (beta*Delta_t);
      a_dt[2]= 1 / (beta*Delta_t);
      a_dt[3]= 1 /(2*beta) - 1;
      a_dt[4]= gamma/beta - 1;
      a_dt[5]= (Delta_t/2) * (gamma/beta - 2);
      a_dt[6]= Delta_t * (1-gamma);
      a_dt[7]= gamma * Delta_t;
      a_dt[8]= 0.0;
      
      break;
      
    case (GENERALIZED_ALPHA):
      
      /*--- Integration constants for Generalized Alpha ---*/
      /*--- Needs to be updated if accounting for structural damping ---*/
      
      //      su2double beta = config->Get_Int_Coeffs(0);
      //      //  su2double gamma =  config->Get_Int_Coeffs(1);
      //      su2double alpha_f = config->Get_Int_Coeffs(2), alpha_m =  config->Get_Int_Coeffs(3);
      //
      //      a_dt[0]= (1 / (beta*pow(Delta_t,2.0))) * ((1 - alpha_m) / (1 - alpha_f)) ;
      //      a_dt[1]= 0.0 ;
      //      a_dt[2]= (1 - alpha_m) / (beta*Delta_t);
      //      a_dt[3]= ((1 - 2*beta)*(1-alpha_m) / (2*beta)) - alpha_m;
      //      a_dt[4]= 0.0;
      //      a_dt[5]= 0.0;
      //      a_dt[6]= Delta_t * (1-delta);
      //      a_dt[7]= delta * Delta_t;
      //      a_dt[8]= (1 - alpha_m) / (beta*pow(Delta_t,2.0));
      
      break;
  }
  
  
}


void CFEASolver::BC_Clamped(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                       unsigned short val_marker) {
  
  unsigned long iPoint, iVertex;
  unsigned long iVar, jVar;
  
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Get node index ---*/
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      if (nDim == 2) {
        Solution[0] = 0.0;  Solution[1] = 0.0;
        Residual[0] = 0.0;  Residual[1] = 0.0;
      }
      else {
        Solution[0] = 0.0;  Solution[1] = 0.0;  Solution[2] = 0.0;
        Residual[0] = 0.0;  Residual[1] = 0.0;  Residual[2] = 0.0;
      }
      
      node[iPoint]->SetSolution(Solution);
      
      if (dynamic) {
        node[iPoint]->SetSolution_Vel(Solution);
        node[iPoint]->SetSolution_Accel(Solution);
      }
      
      
      /*--- Initialize the reaction vector ---*/
      LinSysReact.SetBlock(iPoint, Residual);
      
      LinSysRes.SetBlock(iPoint, Residual);
      LinSysSol.SetBlock(iPoint, Solution);
      
      /*--- STRONG ENFORCEMENT OF THE DISPLACEMENT BOUNDARY CONDITION ---*/
      
      /*--- Delete the columns for a particular node ---*/
      
      for (iVar = 0; iVar < nPoint; iVar++) {
        if (iVar==iPoint) {
          Jacobian.SetBlock(iVar,iPoint,mId_Aux);
        }
        else {
          Jacobian.SetBlock(iVar,iPoint,mZeros_Aux);
        }
      }
      
      /*--- Delete the rows for a particular node ---*/
      for (jVar = 0; jVar < nPoint; jVar++) {
        if (iPoint!=jVar) {
          Jacobian.SetBlock(iPoint,jVar,mZeros_Aux);
        }
      }
      
      /*--- If the problem is dynamic ---*/
      /*--- Enforce that in the previous time step all nodes had 0 U, U', U'' ---*/
      
      if(dynamic) {
        
        node[iPoint]->SetSolution_time_n(Solution);
        node[iPoint]->SetSolution_Vel_time_n(Solution);
        node[iPoint]->SetSolution_Accel_time_n(Solution);
        
      }
      
    }
    
  }
  
}

void CFEASolver::BC_Clamped_Post(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                            unsigned short val_marker) {
  
  unsigned long iPoint, iVertex;
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Get node index ---*/
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    if (nDim == 2) {
      Solution[0] = 0.0;  Solution[1] = 0.0;
    }
    else {
      Solution[0] = 0.0;  Solution[1] = 0.0;  Solution[2] = 0.0;
    }
    
    node[iPoint]->SetSolution(Solution);
    
    if (dynamic) {
      node[iPoint]->SetSolution_Vel(Solution);
      node[iPoint]->SetSolution_Accel(Solution);
    }
    
  }
  
}

void CFEASolver::BC_DispDir(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                            unsigned short val_marker) {

  unsigned short iDim, jDim;

  su2double DispDirVal = config->GetDisp_Dir_Value(config->GetMarker_All_TagBound(val_marker));
  su2double DispDirMult = config->GetDisp_Dir_Multiplier(config->GetMarker_All_TagBound(val_marker));
  su2double *Disp_Dir_Local= config->GetDisp_Dir(config->GetMarker_All_TagBound(val_marker));
  su2double Disp_Dir_Unit[3] = {0.0, 0.0, 0.0}, Disp_Dir[3] = {0.0, 0.0, 0.0};
  su2double Disp_Dir_Mod = 0.0;

  for (iDim = 0; iDim < nDim; iDim++)
    Disp_Dir_Mod += Disp_Dir_Local[iDim]*Disp_Dir_Local[iDim];

  Disp_Dir_Mod = sqrt(Disp_Dir_Mod);

  for (iDim = 0; iDim < nDim; iDim++)
    Disp_Dir_Unit[iDim] = Disp_Dir_Local[iDim] / Disp_Dir_Mod;

  su2double TotalDisp;

  su2double CurrentTime=config->GetCurrent_DynTime();
  su2double ModAmpl = 1.0;

  bool Ramp_Load = config->GetRamp_Load();
  su2double Ramp_Time = config->GetRamp_Time();
  su2double Transfer_Time = 0.0;

  if (Ramp_Load){
    if (Ramp_Time == 0.0)
      ModAmpl = 1.0;
    else
      Transfer_Time = CurrentTime / Ramp_Time;

    switch (config->GetDynamic_LoadTransfer()) {
    case INSTANTANEOUS:
      ModAmpl = 1.0;
      break;
    case POL_ORDER_1:
      ModAmpl = Transfer_Time;
      break;
    case POL_ORDER_3:
      ModAmpl = -2.0 * pow(Transfer_Time,3.0) + 3.0 * pow(Transfer_Time,2.0);
      break;
    case POL_ORDER_5:
      ModAmpl = 6.0 * pow(Transfer_Time, 5.0) - 15.0 * pow(Transfer_Time, 4.0) + 10 * pow(Transfer_Time, 3.0);
      break;
    case SIGMOID_10:
      ModAmpl = (1 / (1+exp(-1.0 * 10.0 * (Transfer_Time - 0.5)) ) );
      break;
    case SIGMOID_20:
      ModAmpl = (1 / (1+exp(-1.0 * 20.0 * (Transfer_Time - 0.5)) ) );
      break;
    }
    ModAmpl = max(ModAmpl,0.0);
    ModAmpl = min(ModAmpl,1.0);
  }
  else{
    ModAmpl = 1.0;
  }

  TotalDisp = ModAmpl * DispDirVal * DispDirMult;

  for (iDim = 0; iDim < nDim; iDim++)
    Disp_Dir[iDim] = TotalDisp * Disp_Dir_Unit[iDim];

  unsigned long iNode, iVertex;
  unsigned long iPoint, jPoint;

  su2double valJacobian_ij_00 = 0.0;
  su2double auxJacobian_ij[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index ---*/

    iNode = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->node[iNode]->GetDomain()) {

      if (nDim == 2) {
        Solution[0] = Disp_Dir[0];  Solution[1] = Disp_Dir[1];
        Residual[0] = Disp_Dir[0];  Residual[1] = Disp_Dir[1];
      }
      else {
        Solution[0] = Disp_Dir[0];  Solution[1] = Disp_Dir[1];  Solution[2] = Disp_Dir[2];
        Residual[0] = Disp_Dir[0];  Residual[1] = Disp_Dir[1];  Residual[2] = Disp_Dir[2];
      }

      /*--- Initialize the reaction vector ---*/

      LinSysRes.SetBlock(iNode, Residual);
      LinSysSol.SetBlock(iNode, Solution);

      /*--- STRONG ENFORCEMENT OF THE DISPLACEMENT BOUNDARY CONDITION ---*/

      /*--- Delete the full row for node iNode ---*/
      for (jPoint = 0; jPoint < nPoint; jPoint++){

        /*--- Check whether the block is non-zero ---*/
        valJacobian_ij_00 = Jacobian.GetBlock(iNode, jPoint,0,0);

        if (valJacobian_ij_00 != 0.0 ){
          if (iNode != jPoint) {
            Jacobian.SetBlock(iNode,jPoint,mZeros_Aux);
          }
          else{
            Jacobian.SetBlock(iNode,jPoint,mId_Aux);
          }
        }
      }

      /*--- Delete the columns for a particular node ---*/

      for (iPoint = 0; iPoint < nPoint; iPoint++){

        /*--- Check if the term K(iPoint, iNode) is 0 ---*/
        valJacobian_ij_00 = Jacobian.GetBlock(iPoint,iNode,0,0);

        /*--- If the node iNode has a crossed dependency with the point iPoint ---*/
        if (valJacobian_ij_00 != 0.0 ){

          /*--- Retrieve the Jacobian term ---*/
          for (iDim = 0; iDim < nDim; iDim++){
            for (jDim = 0; jDim < nDim; jDim++){
              auxJacobian_ij[iDim][jDim] = Jacobian.GetBlock(iPoint,iNode,iDim,jDim);
            }
          }

          /*--- Multiply by the imposed displacement ---*/
          for (iDim = 0; iDim < nDim; iDim++){
            Residual[iDim] = 0.0;
            for (jDim = 0; jDim < nDim; jDim++){
              Residual[iDim] += auxJacobian_ij[iDim][jDim] * Disp_Dir[jDim];
            }
          }

          if (iNode != iPoint) {
            /*--- The term is substracted from the residual (right hand side) ---*/
            LinSysRes.SubtractBlock(iPoint, Residual);
            /*--- The Jacobian term is now set to 0 ---*/
            Jacobian.SetBlock(iPoint,iNode,mZeros_Aux);
          }

        }

      }

    }

  }


}

void CFEASolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,  CNumerics **numerics,unsigned short iMesh) {
  
  unsigned short iVar;
  unsigned long iPoint, total_index;
  
  bool first_iter = (config->GetIntIter() == 0);
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);    // Nonlinear analysis.
  bool disc_adj_fem = (config->GetKind_Solver() == DISC_ADJ_FEM);
  
  su2double solNorm = 0.0, solNorm_recv = 0.0;
  
  if (disc_adj_fem) {

    if (nonlinear_analysis) {

      /*--- For nonlinear discrete adjoint, we have 3 convergence criteria ---*/

      /*--- UTOL = norm(Delta_U(k)): ABSOLUTE, norm of the incremental displacements ------------*/
      /*--- RTOL = norm(Residual(k): ABSOLUTE, norm of the residual (T-F) -----------------------*/
      /*--- ETOL = Delta_U(k) * Residual(k): ABSOLUTE, norm of the product disp * res -----------*/

      Conv_Check[0] = LinSysSol.norm();               // Norm of the delta-solution vector
      Conv_Check[1] = LinSysRes.norm();               // Norm of the residual
      Conv_Check[2] = dotProd(LinSysSol, LinSysRes);  // Position for the energy tolerance

      /*--- MPI solution ---*/

      Set_MPI_Solution(geometry, config);
    }
    else {
      /*--- If the problem is linear, the only check we do is the RMS of the displacements ---*/
      /*---  Compute the residual Ax-f ---*/

      Jacobian.ComputeResidual(LinSysSol, LinSysRes, LinSysAux);

      /*--- Set maximum residual to zero ---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        SetRes_RMS(iVar, 0.0);
        SetRes_Max(iVar, 0.0, 0);
      }

      /*--- Compute the residual ---*/

      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          total_index = iPoint*nVar+iVar;
          AddRes_RMS(iVar, LinSysAux[total_index]*LinSysAux[total_index]);
          AddRes_Max(iVar, fabs(LinSysAux[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
        }
      }

      /*--- MPI solution ---*/

      Set_MPI_Solution(geometry, config);

      /*--- Compute the root mean square residual ---*/

      SetResidual_RMS(geometry, config);

    }

  }
  else {
    if (nonlinear_analysis){

      /*--- If the problem is nonlinear, we have 3 convergence criteria ---*/

      /*--- UTOL = norm(Delta_U(k)) / norm(U(k)) --------------------------*/
      /*--- RTOL = norm(Residual(k)) / norm(Residual(0)) ------------------*/
      /*--- ETOL = Delta_U(k) * Residual(k) / Delta_U(0) * Residual(0) ----*/

      if (first_iter){
        Conv_Ref[0] = 1.0;                                        // Position for the norm of the solution
        Conv_Ref[1] = max(LinSysRes.norm(), EPS);                 // Position for the norm of the residual
        Conv_Ref[2] = max(dotProd(LinSysSol, LinSysRes), EPS);    // Position for the energy tolerance

        /*--- Make sure the computation runs at least 2 iterations ---*/
        Conv_Check[0] = 1.0;
        Conv_Check[1] = 1.0;
        Conv_Check[2] = 1.0;

        /*--- If absolute, we check the norms ---*/
        switch (config->GetResidual_Criteria_FEM()) {
          case RESFEM_ABSOLUTE:
            Conv_Check[0] = LinSysSol.norm();         // Norm of the delta-solution vector
            Conv_Check[1] = LinSysRes.norm();         // Norm of the residual
            Conv_Check[2] = dotProd(LinSysSol, LinSysRes);  // Position for the energy tolerance
            break;
        }
      }
      else {
        /*--- Compute the norm of the solution vector Uk ---*/
        for (iPoint = 0; iPoint < nPointDomain; iPoint++){
          for (iVar = 0; iVar < nVar; iVar++){
            solNorm += node[iPoint]->GetSolution(iVar) * node[iPoint]->GetSolution(iVar);
          }
        }

        // We need to communicate the norm of the solution and compute the RMS throughout the different processors

#ifdef HAVE_MPI
        /*--- We sum the squares of the norms across the different processors ---*/
        SU2_MPI::Allreduce(&solNorm, &solNorm_recv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        solNorm_recv         = solNorm;
#endif

        Conv_Ref[0] = max(sqrt(solNorm_recv), EPS);           // Norm of the solution vector

        switch (config->GetResidual_Criteria_FEM()) {
          case RESFEM_RELATIVE:
            Conv_Check[0] = LinSysSol.norm() / Conv_Ref[0];         // Norm of the delta-solution vector
            Conv_Check[1] = LinSysRes.norm() / Conv_Ref[1];         // Norm of the residual
            Conv_Check[2] = dotProd(LinSysSol, LinSysRes) / Conv_Ref[2];  // Position for the energy tolerance
            break;
          case RESFEM_ABSOLUTE:
            Conv_Check[0] = LinSysSol.norm();         // Norm of the delta-solution vector
            Conv_Check[1] = LinSysRes.norm();         // Norm of the residual
            Conv_Check[2] = dotProd(LinSysSol, LinSysRes);  // Position for the energy tolerance
            break;
          default:
            Conv_Check[0] = LinSysSol.norm() / Conv_Ref[0];         // Norm of the delta-solution vector
            Conv_Check[1] = LinSysRes.norm() / Conv_Ref[1];         // Norm of the residual
            Conv_Check[2] = dotProd(LinSysSol, LinSysRes) / Conv_Ref[2];  // Position for the energy tolerance
            break;
        }

      }

      /*--- MPI solution ---*/

      Set_MPI_Solution(geometry, config);

    } else {

      /*--- If the problem is linear, the only check we do is the RMS of the displacements ---*/

      /*---  Compute the residual Ax-f ---*/

      Jacobian.ComputeResidual(LinSysSol, LinSysRes, LinSysAux);

      /*--- Set maximum residual to zero ---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        SetRes_RMS(iVar, 0.0);
        SetRes_Max(iVar, 0.0, 0);
      }

      /*--- Compute the residual ---*/

      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          total_index = iPoint*nVar+iVar;
          AddRes_RMS(iVar, LinSysAux[total_index]*LinSysAux[total_index]);
          AddRes_Max(iVar, fabs(LinSysAux[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
        }
      }

      /*--- MPI solution ---*/

      Set_MPI_Solution(geometry, config);

      /*--- Compute the root mean square residual ---*/

      SetResidual_RMS(geometry, config);
    }
    
  }
  
}

void CFEASolver::BC_Normal_Displacement(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                                   unsigned short val_marker) { }

void CFEASolver::BC_Normal_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                           unsigned short val_marker) {
  
  /*--- Retrieve the normal pressure and the application conditions for the considered boundary ---*/
  
  su2double NormalLoad = config->GetLoad_Value(config->GetMarker_All_TagBound(val_marker));
  su2double TotalLoad = 0.0;
  
  su2double CurrentTime=config->GetCurrent_DynTime();
  su2double ModAmpl = 1.0;
  
  su2double Ramp_Time = config->GetRamp_Time();

  ModAmpl = Compute_LoadCoefficient(CurrentTime, Ramp_Time, config);

  TotalLoad = ModAmpl * NormalLoad;
  
  /*--- Do only if there is a load applied.
   *--- This reduces the computational cost for cases in which we want boundaries with no load.
   */
  
  if (TotalLoad != 0.0) {
    
    unsigned long iElem;
    unsigned short nNodes = 0;
    su2double Length_Elem_ref = 0.0,  Area_Elem_ref = 0.0;
    su2double Length_Elem_curr = 0.0, Area_Elem_curr = 0.0;
    unsigned long indexNode[4]   = {0,0,0,0};
    unsigned long indexVertex[4] = {0,0,0,0};
    su2double nodeCoord_ref[4][3], nodeCoord_curr[4][3];
    su2double *nodal_normal, nodal_normal_unit[3];
    su2double normal_ref_unit[3], normal_curr_unit[3];
    su2double Norm, dot_Prod;
    su2double val_Coord, val_Sol;
    unsigned short iNode, iDim;
    unsigned long iVertex, iPoint;
    su2double a_ref[3], b_ref[3], AC_ref[3], BD_ref[3];
    su2double a_curr[3], b_curr[3], AC_curr[3], BD_curr[3];
    
    /*--- Determine whether the load conditions are applied in the reference or in the current configuration ---*/
    
    bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
    bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS); // Nonlinear analysis.
    
    for (iNode = 0; iNode < 4; iNode++) {
      for (iDim = 0; iDim < 3; iDim++) {
        nodeCoord_ref[iNode][iDim]  = 0.0;
        nodeCoord_curr[iNode][iDim] = 0.0;
      }
    }
    
    for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {
      
      /*--- Identify the kind of boundary element ---*/
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == LINE)           nNodes = 2;
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE)       nNodes = 3;
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL)  nNodes = 4;
      
      /*--- Retrieve the boundary reference and current coordinates ---*/
      for (iNode = 0; iNode < nNodes; iNode++) {
        indexNode[iNode] = geometry->bound[val_marker][iElem]->GetNode(iNode);
        for (iDim = 0; iDim < nDim; iDim++) {
          val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
          val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
          /*--- Assign values to the container ---*/
          nodeCoord_ref[iNode][iDim]  = val_Coord;
          nodeCoord_curr[iNode][iDim] = val_Sol;
        }
      }
      
      /*--- We need the indices of the vertices, which are "Dual Grid Info" ---*/
      for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
        iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
        for (iNode = 0; iNode < nNodes; iNode++) {
          if (iPoint == indexNode[iNode]) indexVertex[iNode] = iVertex;
        }
      }
      
      /*--- Retrieve the reference normal for one of the points. They go INSIDE the structural domain. ---*/
      nodal_normal = geometry->vertex[val_marker][indexVertex[0]]->GetNormal();
      Norm = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Norm += nodal_normal[iDim]*nodal_normal[iDim];
      }
      Norm = sqrt(Norm);
      for (iDim = 0; iDim < nDim; iDim++) {
        nodal_normal_unit[iDim] = nodal_normal[iDim] / Norm;
      }
      
      /*--- Compute area (3D), and length of the surfaces (2D), and the unitary normal vector in current configuration ---*/
      
      if (nDim == 2) {
        
        /*-- Compute the vector a in reference and current configurations ---*/
        for (iDim = 0; iDim < nDim; iDim++) {
          a_ref[iDim]  = nodeCoord_ref[0][iDim] -nodeCoord_ref[1][iDim];
          a_curr[iDim] = nodeCoord_curr[0][iDim]-nodeCoord_curr[1][iDim];
        }
        
        /*-- Compute the length of the boundary element in reference and current configurations ---*/
        Length_Elem_curr = sqrt(a_curr[0]*a_curr[0]+a_curr[1]*a_curr[1]);
        Length_Elem_ref  = sqrt(a_ref[0]*a_ref[0]+a_ref[1]*a_ref[1]);
        
        /*-- Compute the length of the boundary element in reference and current configurations ---*/
        normal_ref_unit[0] =   a_ref[1] /Length_Elem_ref;
        normal_ref_unit[1] = -(a_ref[0])/Length_Elem_ref;
        
        normal_curr_unit[0] =   a_curr[1] /Length_Elem_curr;
        normal_curr_unit[1] = -(a_curr[0])/Length_Elem_curr;
        
        /*-- Dot product to check the element orientation in the reference configuration ---*/
        dot_Prod = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          dot_Prod += normal_ref_unit[iDim] * nodal_normal_unit[iDim];
        }
        
        /*--- If dot_Prod > 0, the normal goes inside the structural domain. ---*/
        /*--- If dot_Prod < 0, the normal goes outside the structural domain. ---*/
        /*--- We adopt the criteria of the normal going inside the domain, so if dot_Prod < 1, we change the orientation. ---*/
        if (dot_Prod < 0) {
          for (iDim = 0; iDim < nDim; iDim++) {
            normal_ref_unit[iDim]  = -1.0*normal_ref_unit[iDim];
            normal_curr_unit[iDim] = -1.0*normal_curr_unit[iDim];
          }
        }
        
        if (linear_analysis) {
          Residual[0] = (1.0/2.0) * TotalLoad * Length_Elem_ref * normal_ref_unit[0];
          Residual[1] = (1.0/2.0) * TotalLoad * Length_Elem_ref * normal_ref_unit[1];
          
          node[indexNode[0]]->Add_SurfaceLoad_Res(Residual);
          node[indexNode[1]]->Add_SurfaceLoad_Res(Residual);
        }
        else if (nonlinear_analysis) {
          Residual[0] = (1.0/2.0) * TotalLoad * Length_Elem_curr * normal_curr_unit[0];
          Residual[1] = (1.0/2.0) * TotalLoad * Length_Elem_curr * normal_curr_unit[1];
          
          node[indexNode[0]]->Add_SurfaceLoad_Res(Residual);
          node[indexNode[1]]->Add_SurfaceLoad_Res(Residual);
        }
        
      }
      
      if (nDim == 3) {
        
        if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE) {
          
          for (iDim = 0; iDim < nDim; iDim++) {
            a_ref[iDim] = nodeCoord_ref[1][iDim]-nodeCoord_ref[0][iDim];
            b_ref[iDim] = nodeCoord_ref[2][iDim]-nodeCoord_ref[0][iDim];
            
            a_curr[iDim] = nodeCoord_curr[1][iDim]-nodeCoord_curr[0][iDim];
            b_curr[iDim] = nodeCoord_curr[2][iDim]-nodeCoord_curr[0][iDim];
          }
          
          su2double Ni=0, Nj=0, Nk=0;
          
          /*--- Reference configuration ---*/
          Ni = a_ref[1]*b_ref[2] - a_ref[2]*b_ref[1];
          Nj = a_ref[2]*b_ref[0] - a_ref[0]*b_ref[2];
          Nk = a_ref[0]*b_ref[1] - a_ref[1]*b_ref[0];
          
          Area_Elem_ref = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
          
          normal_ref_unit[0] = Ni / Area_Elem_ref;
          normal_ref_unit[1] = Nj / Area_Elem_ref;
          normal_ref_unit[2] = Nk / Area_Elem_ref;
          
          /*--- Current configuration ---*/
          Ni = a_curr[1]*b_curr[2] - a_curr[2]*b_curr[1];
          Nj = a_curr[2]*b_curr[0] - a_curr[0]*b_curr[2];
          Nk = a_curr[0]*b_curr[1] - a_curr[1]*b_curr[0];
          
          Area_Elem_curr = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
          
          normal_curr_unit[0] = Ni / Area_Elem_curr;
          normal_curr_unit[1] = Nj / Area_Elem_curr;
          normal_curr_unit[2] = Nk / Area_Elem_curr;
          
          /*-- Dot product to check the element orientation in the reference configuration ---*/
          dot_Prod = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            dot_Prod += normal_ref_unit[iDim] * nodal_normal_unit[iDim];
          }
          
          /*--- If dot_Prod > 0, the normal goes inside the structural domain. ---*/
          /*--- If dot_Prod < 0, the normal goes outside the structural domain. ---*/
          /*--- We adopt the criteria of the normal going inside the domain, so if dot_Prod < 1, we change the orientation. ---*/
          if (dot_Prod < 0) {
            for (iDim = 0; iDim < nDim; iDim++) {
              normal_ref_unit[iDim]  = -1.0*normal_ref_unit[iDim];
              normal_curr_unit[iDim] = -1.0*normal_curr_unit[iDim];
            }
          }
          
          if (linear_analysis) {
            Residual[0] = (1.0/3.0) * TotalLoad * Area_Elem_ref * normal_ref_unit[0];
            Residual[1] = (1.0/3.0) * TotalLoad * Area_Elem_ref * normal_ref_unit[1];
            Residual[2] = (1.0/3.0) * TotalLoad * Area_Elem_ref * normal_ref_unit[2];
            
            node[indexNode[0]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[1]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[2]]->Add_SurfaceLoad_Res(Residual);
          }
          else if (nonlinear_analysis) {
            Residual[0] = (1.0/3.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[0];
            Residual[1] = (1.0/3.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[1];
            Residual[2] = (1.0/3.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[2];
            
            node[indexNode[0]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[1]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[2]]->Add_SurfaceLoad_Res(Residual);
          }
          
        }
        
        else if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL) {
          
          for (iDim = 0; iDim < nDim; iDim++) {
            AC_ref[iDim] = nodeCoord_ref[2][iDim]-nodeCoord_ref[0][iDim];
            BD_ref[iDim] = nodeCoord_ref[3][iDim]-nodeCoord_ref[1][iDim];
            
            AC_curr[iDim] = nodeCoord_curr[2][iDim]-nodeCoord_curr[0][iDim];
            BD_curr[iDim] = nodeCoord_curr[3][iDim]-nodeCoord_curr[1][iDim];
          }
          
          su2double Ni=0, Nj=0, Nk=0;
          
          /*--- Reference configuration ---*/
          Ni=AC_ref[1]*BD_ref[2]-AC_ref[2]*BD_ref[1];
          Nj=-AC_ref[0]*BD_ref[2]+AC_ref[2]*BD_ref[0];
          Nk=AC_ref[0]*BD_ref[1]-AC_ref[1]*BD_ref[0];
          
          Area_Elem_ref = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
          
          normal_ref_unit[0] = Ni / Area_Elem_ref;
          normal_ref_unit[1] = Nj / Area_Elem_ref;
          normal_ref_unit[2] = Nk / Area_Elem_ref;
          
          /*--- Current configuration ---*/
          Ni=AC_curr[1]*BD_curr[2]-AC_curr[2]*BD_curr[1];
          Nj=-AC_curr[0]*BD_curr[2]+AC_curr[2]*BD_curr[0];
          Nk=AC_curr[0]*BD_curr[1]-AC_curr[1]*BD_curr[0];
          
          Area_Elem_curr = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
          
          normal_curr_unit[0] = Ni / Area_Elem_curr;
          normal_curr_unit[1] = Nj / Area_Elem_curr;
          normal_curr_unit[2] = Nk / Area_Elem_curr;
          
          /*-- Dot product to check the element orientation in the reference configuration ---*/
          dot_Prod = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            dot_Prod += normal_ref_unit[iDim] * nodal_normal_unit[iDim];
          }
          
          /*--- If dot_Prod > 0, the normal goes inside the structural domain. ---*/
          /*--- If dot_Prod < 0, the normal goes outside the structural domain. ---*/
          /*--- We adopt the criteria of the normal going inside the domain, so if dot_Prod < 1, we change the orientation. ---*/
          if (dot_Prod < 0) {
            for (iDim = 0; iDim < nDim; iDim++) {
              normal_ref_unit[iDim]  = -1.0*normal_ref_unit[iDim];
              normal_curr_unit[iDim] = -1.0*normal_curr_unit[iDim];
            }
          }
          
          if (linear_analysis) {
            Residual[0] = (1.0/4.0) * TotalLoad * Area_Elem_ref * normal_ref_unit[0];
            Residual[1] = (1.0/4.0) * TotalLoad * Area_Elem_ref * normal_ref_unit[1];
            Residual[2] = (1.0/4.0) * TotalLoad * Area_Elem_ref * normal_ref_unit[2];
            
            node[indexNode[0]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[1]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[2]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[3]]->Add_SurfaceLoad_Res(Residual);
          }
          else if (nonlinear_analysis) {
            Residual[0] = (1.0/4.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[0];
            Residual[1] = (1.0/4.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[1];
            Residual[2] = (1.0/4.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[2];
            
            node[indexNode[0]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[1]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[2]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[3]]->Add_SurfaceLoad_Res(Residual);
          }
          
        }
        
      }
      
      
    }
    
  }
  
}

void CFEASolver::BC_Dir_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                        unsigned short val_marker) {
  
  su2double a[3], b[3], AC[3], BD[3];
  unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3=0;
  su2double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3= NULL;
  su2double Length_Elem = 0.0, Area_Elem = 0.0;
  unsigned short iDim;
  
  su2double LoadDirVal = config->GetLoad_Dir_Value(config->GetMarker_All_TagBound(val_marker));
  su2double LoadDirMult = config->GetLoad_Dir_Multiplier(config->GetMarker_All_TagBound(val_marker));
  su2double *Load_Dir_Local= config->GetLoad_Dir(config->GetMarker_All_TagBound(val_marker));
  
  su2double TotalLoad;
  
  su2double CurrentTime=config->GetCurrent_DynTime();
  su2double Ramp_Time = config->GetRamp_Time();

  su2double ModAmpl = 1.0;
  
  ModAmpl = Compute_LoadCoefficient(CurrentTime, Ramp_Time, config);

  TotalLoad = ModAmpl * LoadDirVal * LoadDirMult;
  
  /*--- Compute the norm of the vector that was passed in the config file ---*/
  su2double Norm = 1.0;
  if (nDim==2) Norm=sqrt(Load_Dir_Local[0]*Load_Dir_Local[0]+Load_Dir_Local[1]*Load_Dir_Local[1]);
  if (nDim==3) Norm=sqrt(Load_Dir_Local[0]*Load_Dir_Local[0]+Load_Dir_Local[1]*Load_Dir_Local[1]+Load_Dir_Local[2]*Load_Dir_Local[2]);
  
  for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {
    
    Point_0 = geometry->bound[val_marker][iElem]->GetNode(0);     Coord_0 = geometry->node[Point_0]->GetCoord();
    Point_1 = geometry->bound[val_marker][iElem]->GetNode(1);     Coord_1 = geometry->node[Point_1]->GetCoord();
    if (nDim == 3) {
      
      Point_2 = geometry->bound[val_marker][iElem]->GetNode(2);  Coord_2 = geometry->node[Point_2]->GetCoord();
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL) {
        Point_3 = geometry->bound[val_marker][iElem]->GetNode(3);  Coord_3 = geometry->node[Point_3]->GetCoord();
      }
      
    }
    
    /*--- Compute area (3D), and length of the surfaces (2D) ---*/
    
    if (nDim == 2) {
      
      for (iDim = 0; iDim < nDim; iDim++) a[iDim] = Coord_0[iDim]-Coord_1[iDim];
      
      Length_Elem = sqrt(a[0]*a[0]+a[1]*a[1]);
      
    }
    
    if (nDim == 3) {
      
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE) {
        
        for (iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = Coord_1[iDim]-Coord_0[iDim];
          b[iDim] = Coord_2[iDim]-Coord_0[iDim];
        }
        
        su2double Ni=0 , Nj=0, Nk=0;
        
        Ni=a[1]*b[2]-a[2]*b[1];
        Nj=-a[0]*b[2]+a[2]*b[0];
        Nk=a[0]*b[1]-a[1]*b[0];
        
        Area_Elem = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
        
      }
      
      else if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL) {
        
        for (iDim = 0; iDim < nDim; iDim++) {
          AC[iDim] = Coord_2[iDim]-Coord_0[iDim];
          BD[iDim] = Coord_3[iDim]-Coord_1[iDim];
        }
        
        su2double Ni=0 , Nj=0, Nk=0;
        
        Ni=AC[1]*BD[2]-AC[2]*BD[1];
        Nj=-AC[0]*BD[2]+AC[2]*BD[0];
        Nk=AC[0]*BD[1]-AC[1]*BD[0];
        
        Area_Elem = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
        
      }
    }
    
    if (nDim == 2) {
      
      Residual[0] = (1.0/2.0)*Length_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
      Residual[1] = (1.0/2.0)*Length_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
      
      node[Point_0]->Add_SurfaceLoad_Res(Residual);
      node[Point_1]->Add_SurfaceLoad_Res(Residual);
      
    }
    
    else {
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE) {
        
        Residual[0] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
        Residual[1] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
        Residual[2] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[2]/Norm;
        
        node[Point_0]->Add_SurfaceLoad_Res(Residual);
        node[Point_1]->Add_SurfaceLoad_Res(Residual);
        node[Point_2]->Add_SurfaceLoad_Res(Residual);
        
      }
      else if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL) {
        
        Residual[0] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
        Residual[1] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
        Residual[2] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[2]/Norm;
        
        node[Point_0]->Add_SurfaceLoad_Res(Residual);
        node[Point_1]->Add_SurfaceLoad_Res(Residual);
        node[Point_2]->Add_SurfaceLoad_Res(Residual);
        node[Point_3]->Add_SurfaceLoad_Res(Residual);
        
      }
      
    }
    
  }
  
}

void CFEASolver::BC_Sine_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                         unsigned short val_marker) { }

void CFEASolver::BC_Damper(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                      unsigned short val_marker) {

  unsigned short iVar;
  su2double dampValue, dampC;
  su2double dampConstant = config->GetDamper_Constant(config->GetMarker_All_TagBound(val_marker));
  unsigned long nVertex = geometry->GetnVertex(val_marker);

  /*--- The damping is distributed evenly over all the nodes in the marker ---*/
  dampC = dampConstant / (nVertex + EPS);

  unsigned long Point_0, Point_1, Point_2, Point_3;
  unsigned long iElem;

  for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {

    Point_0 = geometry->bound[val_marker][iElem]->GetNode(0);
    Point_1 = geometry->bound[val_marker][iElem]->GetNode(1);

    for (iVar = 0; iVar < nVar; iVar++){

        dampValue = - 1.0 * dampC * node[Point_0]->GetSolution_Vel(iVar);
        node[Point_0]->Set_SurfaceLoad_Res(iVar, dampValue);

        dampValue = - 1.0 * dampC * node[Point_1]->GetSolution_Vel(iVar);
        node[Point_1]->Set_SurfaceLoad_Res(iVar, dampValue);
    }

    if (nDim == 3) {

      Point_2 = geometry->bound[val_marker][iElem]->GetNode(2);

      for (iVar = 0; iVar < nVar; iVar++){
          dampValue = - 1.0 * dampC * node[Point_2]->GetSolution_Vel(iVar);
          node[Point_2]->Set_SurfaceLoad_Res(iVar, dampValue);
      }

      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL) {
        Point_3 = geometry->bound[val_marker][iElem]->GetNode(3);
        for (iVar = 0; iVar < nVar; iVar++){
            dampValue = - 1.0 * dampC * node[Point_3]->GetSolution_Vel(iVar);
            node[Point_3]->Set_SurfaceLoad_Res(iVar, dampValue);
        }
      }

    }

  }


}

void CFEASolver::Integrate_FSI_Loads(CGeometry *geometry, CConfig *config) {

  unsigned short iDim, iNode, nNode;
  unsigned long iPoint, iElem, nElem;

  unsigned short iMarkerInt, nMarkerInt = config->GetMarker_n_ZoneInterface()/2,
                 iMarker, nMarker = config->GetnMarker_All();

  /*--- Temporary storage to store the forces on the element faces ---*/
  vector<su2double> forces;

  /*--- Loop through the FSI interface pairs ---*/
  /*--- 1st pass to compute forces ---*/
  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; ++iMarkerInt) {
    /*--- Find the marker index associated with the pair ---*/
    for (iMarker = 0; iMarker < nMarker; ++iMarker)
      if (config->GetMarker_All_ZoneInterface(iMarker) == iMarkerInt)
        break;
    /*--- The current mpi rank may not have this marker ---*/
    if (iMarker == nMarker) continue;

    nElem = geometry->GetnElem_Bound(iMarker);

    for (iElem = 0; iElem < nElem; ++iElem) {
      /*--- Define the boundary element ---*/
      unsigned long nodes[4];
      su2double coords[4][3];
      bool quad = geometry->bound[iMarker][iElem]->GetVTK_Type() == QUADRILATERAL;
      nNode = quad? 4 : nDim;
      
      for (iNode = 0; iNode < nNode; ++iNode) {
        nodes[iNode] = geometry->bound[iMarker][iElem]->GetNode(iNode);
        for (iDim = 0; iDim < nDim; ++iDim)
          coords[iNode][iDim] = geometry->node[nodes[iNode]]->GetCoord(iDim)+
                                node[nodes[iNode]]->GetSolution(iDim);
      }

      /*--- Compute the area ---*/
      su2double area = 0.0;

      if (nDim == 2)
        area = (coords[0][0]-coords[1][0])*(coords[0][0]-coords[1][0])+
               (coords[0][1]-coords[1][1])*(coords[0][1]-coords[1][1]);

      if (nDim == 3) {
        su2double a[3], b[3], Ni, Nj, Nk;

        if (!quad) { // sides of the triangle
          for (iDim = 0; iDim < 3; iDim++) {
            a[iDim] = coords[1][iDim]-coords[0][iDim];
            b[iDim] = coords[2][iDim]-coords[0][iDim];
          }
        }
        else { // diagonals of the quadrilateral
          for (iDim = 0; iDim < 3; iDim++) {
            a[iDim] = coords[2][iDim]-coords[0][iDim];
            b[iDim] = coords[3][iDim]-coords[1][iDim];
          }
        }
        /*--- Area = 0.5*||a x b|| ---*/
        Ni = a[1]*b[2]-a[2]*b[1];
        Nj =-a[0]*b[2]+a[2]*b[0];
        Nk = a[0]*b[1]-a[1]*b[0];

        area = 0.25*(Ni*Ni+Nj*Nj+Nk*Nk);
      }
      area = sqrt(area);

      /*--- Integrate ---*/
      passivedouble weight = 1.0/nNode;
      su2double force[3] = {0.0, 0.0, 0.0};

      for (iNode = 0; iNode < nNode; ++iNode)
        for (iDim = 0; iDim < nDim; ++iDim)
          force[iDim] += weight*area*node[nodes[iNode]]->Get_FlowTraction(iDim);

      for (iDim = 0; iDim < nDim; ++iDim) forces.push_back(force[iDim]);
    }
  }

  /*--- 2nd pass to set values. This is to account for overlap in the markers. ---*/
  /*--- By putting the integrated values back into the nodes no changes have to be made elsewhere. ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); ++iPoint)
    node[iPoint]->Clear_FlowTraction();
  
  vector<su2double>::iterator force_it = forces.begin();
  
  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; ++iMarkerInt) {
    /*--- Find the marker index associated with the pair ---*/
    for (iMarker = 0; iMarker < nMarker; ++iMarker)
      if (config->GetMarker_All_ZoneInterface(iMarker) == iMarkerInt)
        break;
    /*--- The current mpi rank may not have this marker ---*/
    if (iMarker == nMarker) continue;

    nElem = geometry->GetnElem_Bound(iMarker);

    for (iElem = 0; iElem < nElem; ++iElem) {
      bool quad = geometry->bound[iMarker][iElem]->GetVTK_Type() == QUADRILATERAL;
      nNode = quad? 4 : nDim;
      passivedouble weight = 1.0/nNode;

      su2double force[3];
      for (iDim = 0; iDim < nDim; ++iDim) force[iDim] = *(force_it++)*weight;

      for (iNode = 0; iNode < nNode; ++iNode) {
        iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
        node[iPoint]->Add_FlowTraction(force);
      }
    }
  }

#ifdef HAVE_MPI
  /*--- Perform a global reduction, every rank will get the nodal values of all halo elements ---*/
  /*--- This should be cheaper than the "normal" way, since very few points are both halo and interface ---*/
  vector<unsigned long> halo_point_loc, halo_point_glb;
  vector<su2double> halo_force;

  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; ++iMarkerInt) {
    /*--- Find the marker index associated with the pair ---*/
    for (iMarker = 0; iMarker < nMarker; ++iMarker)
      if (config->GetMarker_All_ZoneInterface(iMarker) == iMarkerInt)
        break;
    /*--- The current mpi rank may not have this marker ---*/
    if (iMarker == nMarker) continue;

    nElem = geometry->GetnElem_Bound(iMarker);

    for (iElem = 0; iElem < nElem; ++iElem) {
      bool quad = geometry->bound[iMarker][iElem]->GetVTK_Type() == QUADRILATERAL;
      nNode = quad? 4 : nDim;

      /*--- If this is an halo element we share the nodal forces ---*/
      for (iNode = 0; iNode < nNode; ++iNode)
        if (!geometry->node[geometry->bound[iMarker][iElem]->GetNode(iNode)]->GetDomain())
          break;

      if (iNode < nNode) {
        for (iNode = 0; iNode < nNode; ++iNode) {
          iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
          /*--- local is for when later we update the values in this rank ---*/
          halo_point_loc.push_back(iPoint);
          halo_point_glb.push_back(geometry->node[iPoint]->GetGlobalIndex());
          for (iDim = 0; iDim < nDim; ++iDim)
            halo_force.push_back(node[iPoint]->Get_FlowTraction(iDim));
        }
      }
    }
  }
  /*--- Determine the size of the arrays we need ---*/
  unsigned long nHaloLoc = halo_point_loc.size();
  unsigned long nHaloMax;
  MPI_Allreduce(&nHaloLoc,&nHaloMax,1,MPI_UNSIGNED_LONG,MPI_MAX,MPI_COMM_WORLD);

  /*--- Shared arrays, all the: number of halo points; halo point global indices; respective forces ---*/
  unsigned long *halo_point_num = new unsigned long[size];
  unsigned long *halo_point_all = new unsigned long[size*nHaloMax];
  su2double *halo_force_all = new su2double[size*nHaloMax*nDim];
  
  /*--- If necessary put dummy values in halo_point_glb to get a valid pointer ---*/
  if (halo_point_glb.empty()) halo_point_glb.resize(1);
  /*--- Pad halo_force to avoid (observed) issues in the adjoint when nHaloLoc!=nHaloMax ---*/
  while (halo_force.size() < nHaloMax*nDim) halo_force.push_back(0.0);
  
  MPI_Allgather(&nHaloLoc,1,MPI_UNSIGNED_LONG,halo_point_num,1,MPI_UNSIGNED_LONG,MPI_COMM_WORLD);
  MPI_Allgather(&halo_point_glb[0],nHaloLoc,MPI_UNSIGNED_LONG,halo_point_all,nHaloMax,MPI_UNSIGNED_LONG,MPI_COMM_WORLD);
  SU2_MPI::Allgather(&halo_force[0],nHaloMax*nDim,MPI_DOUBLE,halo_force_all,nHaloMax*nDim,MPI_DOUBLE,MPI_COMM_WORLD);

  /*--- Find shared points with other ranks and update our values ---*/
  for (int proc = 0; proc < size; ++proc)
  if (proc != rank) {
    unsigned long offset = proc*nHaloMax;
    for (iPoint = 0; iPoint < halo_point_num[proc]; ++iPoint) {
      unsigned long iPoint_glb = halo_point_all[offset+iPoint];
      ptrdiff_t pos = find(halo_point_glb.begin(),halo_point_glb.end(),iPoint_glb)-halo_point_glb.begin();
      if (pos < long(halo_point_glb.size())) {
        unsigned long iPoint_loc = halo_point_loc[pos];
        node[iPoint_loc]->Add_FlowTraction(&halo_force_all[(offset+iPoint)*nDim]);
      }
    }
  }

  delete [] halo_point_num;
  delete [] halo_point_all;
  delete [] halo_force_all;
#endif
}

su2double CFEASolver::Compute_LoadCoefficient(su2double CurrentTime, su2double RampTime, CConfig *config){

  su2double LoadCoeff = 1.0;

  bool Ramp_Load = config->GetRamp_Load();
  bool Sine_Load = config->GetSine_Load();
  bool Ramp_And_Release = config->GetRampAndRelease_Load();

  su2double SineAmp = 0.0, SineFreq = 0.0, SinePhase = 0.0;

  su2double TransferTime = 1.0;

  bool restart = config->GetRestart(); // Restart analysis
  bool fsi = config->GetFSI_Simulation();  // FSI simulation.
  bool stat_fsi = (config->GetDynamic_Analysis() == STATIC);

  /*--- This offset introduces the ramp load in dynamic cases starting from the restart point. ---*/
  bool offset = (restart && fsi && (!stat_fsi));
  su2double DeltaT = config->GetDelta_DynTime();
  su2double OffsetTime = 0.0;
  OffsetTime = DeltaT * (config->GetDyn_RestartIter()-1);

  /*--- Polynomial functions from https://en.wikipedia.org/wiki/Smoothstep ---*/

  if ((Ramp_Load) && (RampTime > 0.0)){

    TransferTime = CurrentTime / RampTime;

    if (offset) TransferTime = (CurrentTime - OffsetTime) / RampTime;

    switch (config->GetDynamic_LoadTransfer()) {
    case INSTANTANEOUS:
      LoadCoeff = 1.0;
      break;
    case POL_ORDER_1:
      LoadCoeff = TransferTime;
      break;
    case POL_ORDER_3:
      LoadCoeff = -2.0 * pow(TransferTime,3.0) + 3.0 * pow(TransferTime,2.0);
      break;
    case POL_ORDER_5:
      LoadCoeff = 6.0 * pow(TransferTime, 5.0) - 15.0 * pow(TransferTime, 4.0) + 10 * pow(TransferTime, 3.0);
      break;
    case SIGMOID_10:
      LoadCoeff = (1 / (1+exp(-1.0 * 10.0 * (TransferTime - 0.5)) ) );
      break;
    case SIGMOID_20:
      LoadCoeff = (1 / (1+exp(-1.0 * 20.0 * (TransferTime - 0.5)) ) );
      break;
    }

    if (TransferTime > 1.0) LoadCoeff = 1.0;

    LoadCoeff = max(LoadCoeff,0.0);
    LoadCoeff = min(LoadCoeff,1.0);

  }
  else if (Sine_Load){

      /*--- Retrieve amplitude, frequency (Hz) and phase (rad) ---*/
      SineAmp   = config->GetLoad_Sine()[0];
      SineFreq  = config->GetLoad_Sine()[1];
      SinePhase = config->GetLoad_Sine()[2];

      LoadCoeff = SineAmp * sin(2*PI_NUMBER*SineFreq*CurrentTime + SinePhase);
  }

  /*--- Add possibility to release the load after the ramp---*/
  if ((Ramp_And_Release) && (CurrentTime >  RampTime)){
    LoadCoeff = 0.0;
  }

  /*--- Store the force coefficient ---*/

  SetForceCoeff(LoadCoeff);

  return LoadCoeff;


}

void CFEASolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

void CFEASolver::ImplicitNewmark_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned long iPoint, jPoint;
  unsigned short iVar, jVar;
  
  bool first_iter = (config->GetIntIter() == 0);
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.
  bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
  bool newton_raphson = (config->GetKind_SpaceIteScheme_FEA() == NEWTON_RAPHSON);    // Newton-Raphson method
  bool fsi = config->GetFSI_Simulation();                        // FSI simulation.
  
  bool body_forces = config->GetDeadLoad();                      // Body forces (dead loads).
  
  bool incremental_load = config->GetIncrementalLoad();
  
  if (!dynamic) {
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      /*--- Add the external contribution to the residual    ---*/
      /*--- (the terms that are constant over the time step) ---*/
      if (incremental_load) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = loadIncrement * node[iPoint]->Get_SurfaceLoad_Res(iVar);
        }
      }
      else {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = node[iPoint]->Get_SurfaceLoad_Res(iVar);
        }
        //Res_Ext_Surf = node[iPoint]->Get_SurfaceLoad_Res();
      }
      
      LinSysRes.AddBlock(iPoint, Res_Ext_Surf);
      
      /*--- Add the contribution to the residual due to body forces ---*/
      
      if (body_forces) {
        if (incremental_load) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = loadIncrement * node[iPoint]->Get_BodyForces_Res(iVar);
          }
        }
        else {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = node[iPoint]->Get_BodyForces_Res(iVar);
          }
          //Res_Dead_Load = node[iPoint]->Get_BodyForces_Res();
        }
        
        LinSysRes.AddBlock(iPoint, Res_Dead_Load);
      }

      /*---  Add the contribution to the residual due to flow loads (FSI contribution) ---*/
      if (fsi) {
        if (incremental_load){
          for (iVar = 0; iVar < nVar; iVar++){
            Res_FSI_Cont[iVar] = loadIncrement * node[iPoint]->Get_FlowTraction(iVar);
          }
        }
        else {
            for (iVar = 0; iVar < nVar; iVar++){
              Res_FSI_Cont[iVar] = node[iPoint]->Get_FlowTraction(iVar);
            }
        }
        LinSysRes.AddBlock(iPoint, Res_FSI_Cont);
      }
    }
    
  }
  
  if (dynamic) {
    
    /*--- Add the mass matrix contribution to the Jacobian ---*/
    
    /*
     * If the problem is nonlinear, we need to add the Mass Matrix contribution to the Jacobian at the beginning
     * of each time step. If the solution method is Newton Rapshon, we repeat this step at the beginning of each
     * iteration, as the Jacobian is recomputed
     *
     * If the problem is linear, we add the Mass Matrix contribution to the Jacobian everytime because for
     * correct differentiation the Jacobian is recomputed every time step.
     *
     */
    if ((nonlinear_analysis && (newton_raphson || first_iter)) || linear_analysis) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (jPoint = 0; jPoint < nPoint; jPoint++) {
          for(iVar = 0; iVar < nVar; iVar++) {
            for (jVar = 0; jVar < nVar; jVar++) {
              Jacobian_ij[iVar][jVar] = a_dt[0] * MassMatrix.GetBlock(iPoint, jPoint, iVar, jVar);
            }
          }
          Jacobian.AddBlock(iPoint, jPoint, Jacobian_ij);
        }
      }
    }
    
    
    /*--- Loop over all points, and set aux vector TimeRes_Aux = a0*U+a2*U'+a3*U'' ---*/

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Residual[iVar] =   a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)        //a0*U(t)
        - a_dt[0]*node[iPoint]->GetSolution(iVar)           //a0*U(t+dt)(k-1)
        + a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)    //a2*U'(t)
        + a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);  //a3*U''(t)
      }
      TimeRes_Aux.SetBlock(iPoint, Residual);
    }
    
    /*--- Once computed, compute M*TimeRes_Aux ---*/
    MassMatrix.MatrixVectorProduct(TimeRes_Aux,TimeRes,geometry,config);
    /*--- Add the components of M*TimeRes_Aux to the residual R(t+dt) ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      
      /*--- Dynamic contribution ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Time_Cont[iVar] = TimeRes.GetBlock(iPoint, iVar);
      }
      //Res_Time_Cont = TimeRes.GetBlock(iPoint);
      LinSysRes.AddBlock(iPoint, Res_Time_Cont);
      
      /*--- External surface load contribution ---*/
      if (incremental_load) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = loadIncrement * node[iPoint]->Get_SurfaceLoad_Res(iVar);
        }
      }
      else {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = node[iPoint]->Get_SurfaceLoad_Res(iVar);
        }
        //Res_Ext_Surf = node[iPoint]->Get_SurfaceLoad_Res();
      }
      LinSysRes.AddBlock(iPoint, Res_Ext_Surf);
      
      
      /*--- Body forces contribution (dead load) ---*/
      
      if (body_forces) {
        if (incremental_load) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = loadIncrement * node[iPoint]->Get_BodyForces_Res(iVar);
          }
        }
        else {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = node[iPoint]->Get_BodyForces_Res(iVar);
          }
          //Res_Dead_Load = node[iPoint]->Get_BodyForces_Res();
        }
        
        LinSysRes.AddBlock(iPoint, Res_Dead_Load);
      }
      
      /*--- FSI contribution (flow loads) ---*/
      if (fsi) {
        if (incremental_load) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_FSI_Cont[iVar] = loadIncrement * node[iPoint]->Get_FlowTraction(iVar);
          }
        }
        else {
          Res_FSI_Cont = node[iPoint]->Get_FlowTraction();
        }
        LinSysRes.AddBlock(iPoint, Res_FSI_Cont);
      }
    }
  }
  
  
}

void CFEASolver::ImplicitNewmark_Update(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint;
  
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);          // Dynamic simulations.
  
  /*--- Update solution ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    for (iVar = 0; iVar < nVar; iVar++) {
      
      /*--- Displacements component of the solution ---*/
       
      node[iPoint]->Add_DeltaSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
      
    }
    
  }
  
  if (dynamic) {
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      for (iVar = 0; iVar < nVar; iVar++) {
        
        /*--- Acceleration component of the solution ---*/
        /*--- U''(t+dt) = a0*(U(t+dt)-U(t))+a2*(U'(t))+a3*(U''(t)) ---*/
        
        Solution[iVar]=a_dt[0]*(node[iPoint]->GetSolution(iVar) -
                                node[iPoint]->GetSolution_time_n(iVar)) -
        a_dt[2]* node[iPoint]->GetSolution_Vel_time_n(iVar) -
        a_dt[3]* node[iPoint]->GetSolution_Accel_time_n(iVar);
      }
      
      /*--- Set the acceleration in the node structure ---*/
      
      node[iPoint]->SetSolution_Accel(Solution);
      
      for (iVar = 0; iVar < nVar; iVar++) {
        
        /*--- Velocity component of the solution ---*/
        /*--- U'(t+dt) = U'(t)+ a6*(U''(t)) + a7*(U''(t+dt)) ---*/
        
        Solution[iVar]=node[iPoint]->GetSolution_Vel_time_n(iVar)+
        a_dt[6]* node[iPoint]->GetSolution_Accel_time_n(iVar) +
        a_dt[7]* node[iPoint]->GetSolution_Accel(iVar);
        
      }
      
      /*--- Set the velocity in the node structure ---*/
      
      node[iPoint]->SetSolution_Vel(Solution);
      
    }
    
  }
  
  /*--- Perform the MPI communication of the solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  
}

void CFEASolver::ImplicitNewmark_Relaxation(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint;
  su2double *valSolutionPred;
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);
  
  /*--- Update solution and set it to be the solution after applying relaxation---*/
  
  for (iPoint=0; iPoint < nPointDomain; iPoint++) {
    
    valSolutionPred = node[iPoint]->GetSolution_Pred();
    
    node[iPoint]->SetSolution(valSolutionPred);
  }
  
  if (dynamic){
    
    /*--- Compute velocities and accelerations ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      for (iVar = 0; iVar < nVar; iVar++) {

        /*--- Acceleration component of the solution ---*/
        /*--- U''(t+dt) = a0*(U(t+dt)-U(t))+a2*(U'(t))+a3*(U''(t)) ---*/

        Solution[iVar]=a_dt[0]*(node[iPoint]->GetSolution(iVar) -
            node[iPoint]->GetSolution_time_n(iVar)) -
            a_dt[2]* node[iPoint]->GetSolution_Vel_time_n(iVar) -
            a_dt[3]* node[iPoint]->GetSolution_Accel_time_n(iVar);
      }

      /*--- Set the acceleration in the node structure ---*/

      node[iPoint]->SetSolution_Accel(Solution);

      for (iVar = 0; iVar < nVar; iVar++) {

        /*--- Velocity component of the solution ---*/
        /*--- U'(t+dt) = U'(t)+ a6*(U''(t)) + a7*(U''(t+dt)) ---*/

        Solution[iVar]=node[iPoint]->GetSolution_Vel_time_n(iVar)+
            a_dt[6]* node[iPoint]->GetSolution_Accel_time_n(iVar) +
            a_dt[7]* node[iPoint]->GetSolution_Accel(iVar);

      }

      /*--- Set the velocity in the node structure ---*/

      node[iPoint]->SetSolution_Vel(Solution);

    }
  
  }
  
  /*--- Perform the MPI communication of the solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*--- After the solution has been communicated, set the 'old' predicted solution as the solution ---*/
  /*--- Loop over n points (as we have already communicated everything ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      node[iPoint]->SetSolution_Pred_Old(iVar,node[iPoint]->GetSolution(iVar));
    }
  }
  
  
}


void CFEASolver::GeneralizedAlpha_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned long iPoint, jPoint;
  unsigned short iVar, jVar;
  
  bool first_iter = (config->GetIntIter() == 0);
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.
  bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
  bool newton_raphson = (config->GetKind_SpaceIteScheme_FEA() == NEWTON_RAPHSON);    // Newton-Raphson method
  bool fsi = config->GetFSI_Simulation();                        // FSI simulation.
  
  bool body_forces = config->GetDeadLoad();                      // Body forces (dead loads).
  
  su2double alpha_f = config->Get_Int_Coeffs(2);
  
  bool incremental_load = config->GetIncrementalLoad();
  
  if (!dynamic) {
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      /*--- Add the external contribution to the residual    ---*/
      /*--- (the terms that are constant over the time step) ---*/
      if (incremental_load) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = loadIncrement * node[iPoint]->Get_SurfaceLoad_Res(iVar);
        }
      }
      else {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = node[iPoint]->Get_SurfaceLoad_Res(iVar);
        }
        //Res_Ext_Surf = node[iPoint]->Get_SurfaceLoad_Res();
      }
      
      LinSysRes.AddBlock(iPoint, Res_Ext_Surf);
      
      /*--- Add the contribution to the residual due to body forces ---*/
      
      if (body_forces) {
        if (incremental_load) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = loadIncrement * node[iPoint]->Get_BodyForces_Res(iVar);
          }
        }
        else {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = node[iPoint]->Get_BodyForces_Res(iVar);
          }
          //Res_Dead_Load = node[iPoint]->Get_BodyForces_Res();
        }
        
        LinSysRes.AddBlock(iPoint, Res_Dead_Load);
      }
      
    }
    
  }
  
  if (dynamic) {
    
    /*--- Add the mass matrix contribution to the Jacobian ---*/
    
    /*
     * If the problem is nonlinear, we need to add the Mass Matrix contribution to the Jacobian at the beginning
     * of each time step. If the solution method is Newton Rapshon, we repeat this step at the beginning of each
     * iteration, as the Jacobian is recomputed
     *
     * If the problem is linear, we add the Mass Matrix contribution to the Jacobian everytime because for
     * correct differentiation the Jacobian is recomputed every time step.
     *
     */
    if ((nonlinear_analysis && (newton_raphson || first_iter)) || linear_analysis) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (jPoint = 0; jPoint < nPoint; jPoint++) {
          for(iVar = 0; iVar < nVar; iVar++) {
            for (jVar = 0; jVar < nVar; jVar++) {
              Jacobian_ij[iVar][jVar] = a_dt[0] * MassMatrix.GetBlock(iPoint, jPoint, iVar, jVar);
            }
          }
          Jacobian.AddBlock(iPoint, jPoint, Jacobian_ij);
        }
      }
    }
    
    
    /*--- Loop over all points, and set aux vector TimeRes_Aux = a0*U+a2*U'+a3*U'' ---*/
    
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Residual[iVar] =   a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)        //a0*U(t)
        - a_dt[0]*node[iPoint]->GetSolution(iVar)           //a0*U(t+dt)(k-1)
        + a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)    //a2*U'(t)
        + a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);  //a3*U''(t)
      }
      TimeRes_Aux.SetBlock(iPoint, Residual);
    }

    /*--- Once computed, compute M*TimeRes_Aux ---*/
    MassMatrix.MatrixVectorProduct(TimeRes_Aux,TimeRes,geometry,config);
    /*--- Add the components of M*TimeRes_Aux to the residual R(t+dt) ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      /*--- Dynamic contribution ---*/
      //Res_Time_Cont = TimeRes.GetBlock(iPoint);
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Time_Cont[iVar] = TimeRes.GetBlock(iPoint, iVar);
      }
      LinSysRes.AddBlock(iPoint, Res_Time_Cont);
      /*--- External surface load contribution ---*/
      if (incremental_load) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = loadIncrement * ( (1 - alpha_f) * node[iPoint]->Get_SurfaceLoad_Res(iVar) +
                                                alpha_f  * node[iPoint]->Get_SurfaceLoad_Res_n(iVar) );
        }
      }
      else {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = (1 - alpha_f) * node[iPoint]->Get_SurfaceLoad_Res(iVar) +
          alpha_f  * node[iPoint]->Get_SurfaceLoad_Res_n(iVar);
        }
      }
      LinSysRes.AddBlock(iPoint, Res_Ext_Surf);
      
      /*--- Add the contribution to the residual due to body forces.
       *--- It is constant over time, so it's not necessary to distribute it. ---*/
      
      if (body_forces) {
        if (incremental_load) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = loadIncrement * node[iPoint]->Get_BodyForces_Res(iVar);
          }
        }
        else {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = node[iPoint]->Get_BodyForces_Res(iVar);
          }
          //Res_Dead_Load = node[iPoint]->Get_BodyForces_Res();
        }
        
        LinSysRes.AddBlock(iPoint, Res_Dead_Load);
      }
      
      /*--- Add FSI contribution ---*/
      if (fsi) {
        if (incremental_load) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_FSI_Cont[iVar] = loadIncrement * ( (1 - alpha_f) * node[iPoint]->Get_FlowTraction(iVar) +
                                                  alpha_f  * node[iPoint]->Get_FlowTraction_n(iVar) );
          }
        }
        else {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_FSI_Cont[iVar] = (1 - alpha_f) * node[iPoint]->Get_FlowTraction(iVar) +
            alpha_f  * node[iPoint]->Get_FlowTraction_n(iVar);
          }
        }
        LinSysRes.AddBlock(iPoint, Res_FSI_Cont);
      }
    }
  }
  
}

void CFEASolver::GeneralizedAlpha_UpdateDisp(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint;
  
  /*--- Update solution ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    for (iVar = 0; iVar < nVar; iVar++) {
      
      /*--- Displacements component of the solution ---*/
      
      node[iPoint]->Add_DeltaSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
      
    }
    
  }
  
  /*--- Perform the MPI communication of the solution, displacements only ---*/
  
  Set_MPI_Solution_DispOnly(geometry, config);
  
}

void CFEASolver::GeneralizedAlpha_UpdateSolution(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint;
  
  su2double alpha_f = config->Get_Int_Coeffs(2), alpha_m =  config->Get_Int_Coeffs(3);
  
  /*--- Compute solution at t_n+1, and update velocities and accelerations ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    for (iVar = 0; iVar < nVar; iVar++) {
      
      /*--- Compute the solution from the previous time step and the solution computed at t+1-alpha_f ---*/
      /*--- U(t+dt) = 1/alpha_f*(U(t+1-alpha_f)-alpha_f*U(t)) ---*/
      
      Solution[iVar]=(1 / (1 - alpha_f))*(node[iPoint]->GetSolution(iVar) -
                                          alpha_f * node[iPoint]->GetSolution_time_n(iVar));
      

    }
    
    /*--- Set the solution in the node structure ---*/
    
    node[iPoint]->SetSolution(Solution);
    
    for (iVar = 0; iVar < nVar; iVar++) {
      
      /*--- Acceleration component of the solution ---*/
      /*--- U''(t+dt-alpha_m) = a8*(U(t+dt)-U(t))+a2*(U'(t))+a3*(U''(t)) ---*/
      
      Solution_Interm[iVar]=a_dt[8]*( node[iPoint]->GetSolution(iVar) -
                                     node[iPoint]->GetSolution_time_n(iVar)) -
      a_dt[2]* node[iPoint]->GetSolution_Vel_time_n(iVar) -
      a_dt[3]* node[iPoint]->GetSolution_Accel_time_n(iVar);
      
      /*--- Compute the solution from the previous time step and the solution computed at t+1-alpha_f ---*/
      /*--- U''(t+dt) = 1/alpha_m*(U''(t+1-alpha_m)-alpha_m*U''(t)) ---*/
      
      Solution[iVar]=(1 / (1 - alpha_m))*(Solution_Interm[iVar] - alpha_m * node[iPoint]->GetSolution_Accel_time_n(iVar));
    }
    
    /*--- Set the acceleration in the node structure ---*/
    
    node[iPoint]->SetSolution_Accel(Solution);
    
    for (iVar = 0; iVar < nVar; iVar++) {
      
      /*--- Velocity component of the solution ---*/
      /*--- U'(t+dt) = U'(t)+ a6*(U''(t)) + a7*(U''(t+dt)) ---*/
      
      Solution[iVar]=node[iPoint]->GetSolution_Vel_time_n(iVar)+
      a_dt[6]* node[iPoint]->GetSolution_Accel_time_n(iVar) +
      a_dt[7]* node[iPoint]->GetSolution_Accel(iVar);
      
    }
    
    /*--- Set the velocity in the node structure ---*/
    
    node[iPoint]->SetSolution_Vel(Solution);
    
  }
  
  /*--- Perform the MPI communication of the solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
}

void CFEASolver::GeneralizedAlpha_UpdateLoads(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned long iPoint;
  bool fsi = config->GetFSI_Simulation();
  
  /*--- Set the load conditions of the time step n+1 as the load conditions for time step n ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    node[iPoint]->Set_SurfaceLoad_Res_n();
    if (fsi) node[iPoint]->Set_FlowTraction_n();
  }
  
}

void CFEASolver::Solve_System(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned long IterLinSol = 0, iPoint, total_index;
  unsigned short iVar;
  
  /*--- Initialize residual and solution at the ghost points ---*/
  
  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }
    
  }
  
  CSysSolve femSystem;
  IterLinSol = femSystem.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
  
  /*--- The the number of iterations of the linear solver ---*/
  
  SetIterLinSolver(IterLinSol);
  
}


void CFEASolver::PredictStruct_Displacement(CGeometry **fea_geometry,
                                                       CConfig *fea_config, CSolver ***fea_solution) {
  
  unsigned short predOrder = fea_config->GetPredictorOrder();
  su2double Delta_t = fea_config->GetDelta_DynTime();
  unsigned long iPoint, iDim;
  su2double *solDisp, *solVel, *solVel_tn, *valPred;

  //To nPointDomain: we need to communicate the predicted solution after setting it
  for (iPoint=0; iPoint < nPointDomain; iPoint++) {
    if (predOrder==0) fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Pred();
    else if (predOrder==1) {
      
      solDisp = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
      solVel = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel();
      valPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
      
      for (iDim=0; iDim < nDim; iDim++) {
        valPred[iDim] = solDisp[iDim] + Delta_t*solVel[iDim];
      }
      
    }
    else if (predOrder==2) {
      
      solDisp = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
      solVel = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel();
      solVel_tn = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel_time_n();
      valPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
      
      for (iDim=0; iDim < nDim; iDim++) {
        valPred[iDim] = solDisp[iDim] + 0.5*Delta_t*(3*solVel[iDim]-solVel_tn[iDim]);
      }
      
    }
    else {
      cout<< "Higher order predictor not implemented. Solving with order 0." << endl;
      fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Pred();
    }
  }
  
}

void CFEASolver::ComputeAitken_Coefficient(CGeometry **fea_geometry, CConfig *fea_config,
                                                      CSolver ***fea_solution, unsigned long iOuterIter) {
  
  unsigned long iPoint, iDim;
  su2double rbuf_numAitk = 0, sbuf_numAitk = 0;
  su2double rbuf_denAitk = 0, sbuf_denAitk = 0;
  
  su2double *dispPred, *dispCalc, *dispPred_Old, *dispCalc_Old;
  su2double deltaU[3] = {0.0, 0.0, 0.0}, deltaU_p1[3] = {0.0, 0.0, 0.0};
  su2double delta_deltaU[3] = {0.0, 0.0, 0.0};
  su2double CurrentTime=fea_config->GetCurrent_DynTime();
  su2double WAitkDyn_tn1, WAitkDyn_Max, WAitkDyn_Min, WAitkDyn;
  
  unsigned short RelaxMethod_FSI = fea_config->GetRelaxation_Method_FSI();
  
  ofstream historyFile_FSI;
  bool writeHistFSI = fea_config->GetWrite_Conv_FSI();
  if (writeHistFSI && (rank == MASTER_NODE)) {
    char cstrFSI[200];
    string filenameHistFSI = fea_config->GetConv_FileName_FSI();
    strcpy (cstrFSI, filenameHistFSI.data());
    historyFile_FSI.open (cstrFSI, std::ios_base::app);
  }
  
  
  /*--- Only when there is movement, and a dynamic coefficient is requested, it makes sense to compute the Aitken's coefficient ---*/
  
    
    if (RelaxMethod_FSI == NO_RELAXATION) {
      
      if (writeHistFSI && (rank == MASTER_NODE)) {
        
        SetWAitken_Dyn(1.0);
        
        if (iOuterIter == 0) historyFile_FSI << " " << endl ;
        historyFile_FSI << setiosflags(ios::fixed) << setprecision(4) << CurrentTime << "," ;
        historyFile_FSI << setiosflags(ios::fixed) << setprecision(1) << iOuterIter << "," ;
        if (iOuterIter == 0) historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << 1.0 ;
        else historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << 1.0 << "," ;
      }
      
    }
    else if (RelaxMethod_FSI == FIXED_PARAMETER) {
      
      if (writeHistFSI && (rank == MASTER_NODE)) {
        
        SetWAitken_Dyn(fea_config->GetAitkenStatRelax());
        
        if (iOuterIter == 0) historyFile_FSI << " " << endl ;
        historyFile_FSI << setiosflags(ios::fixed) << setprecision(4) << CurrentTime << "," ;
        historyFile_FSI << setiosflags(ios::fixed) << setprecision(1) << iOuterIter << "," ;
        if (iOuterIter == 0) historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << fea_config->GetAitkenStatRelax() ;
        else historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << fea_config->GetAitkenStatRelax() << "," ;
      }
      
    }
    else if (RelaxMethod_FSI == AITKEN_DYNAMIC) {
      
      if (iOuterIter == 0) {
        
        WAitkDyn_tn1 = GetWAitken_Dyn_tn1();
        WAitkDyn_Max = fea_config->GetAitkenDynMaxInit();
        WAitkDyn_Min = fea_config->GetAitkenDynMinInit();
        
        WAitkDyn = min(WAitkDyn_tn1, WAitkDyn_Max);
        WAitkDyn = max(WAitkDyn, WAitkDyn_Min);
        
        SetWAitken_Dyn(WAitkDyn);
        if (writeHistFSI && (rank == MASTER_NODE)) {
          if (iOuterIter == 0) historyFile_FSI << " " << endl ;
          historyFile_FSI << setiosflags(ios::fixed) << setprecision(4) << CurrentTime << "," ;
          historyFile_FSI << setiosflags(ios::fixed) << setprecision(1) << iOuterIter << "," ;
          historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << WAitkDyn ;
        }
        
      }
      else {
        // To nPointDomain; we need to communicate the values
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
          
          dispPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
          dispPred_Old = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred_Old();
          dispCalc = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
          dispCalc_Old = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Old();
          
          for (iDim = 0; iDim < nDim; iDim++) {
            
            /*--- Compute the deltaU and deltaU_n+1 ---*/
            deltaU[iDim] = dispCalc_Old[iDim] - dispPred_Old[iDim];
            deltaU_p1[iDim] = dispCalc[iDim] - dispPred[iDim];
            
            /*--- Compute the difference ---*/
            delta_deltaU[iDim] = deltaU_p1[iDim] - deltaU[iDim];
            
            /*--- Add numerator and denominator ---*/
            sbuf_numAitk += deltaU[iDim] * delta_deltaU[iDim];
            sbuf_denAitk += delta_deltaU[iDim] * delta_deltaU[iDim];
            
          }
          
        }
        
#ifdef HAVE_MPI
        SU2_MPI::Allreduce(&sbuf_numAitk, &rbuf_numAitk, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(&sbuf_denAitk, &rbuf_denAitk, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        rbuf_numAitk = sbuf_numAitk;
        rbuf_denAitk = sbuf_denAitk;
#endif
        
        WAitkDyn = GetWAitken_Dyn();
        
        if (rbuf_denAitk > 1E-15) {
          WAitkDyn = - 1.0 * WAitkDyn * rbuf_numAitk / rbuf_denAitk ;
        }
        
        WAitkDyn = max(WAitkDyn, 0.1);
        WAitkDyn = min(WAitkDyn, 1.0);
        
        SetWAitken_Dyn(WAitkDyn);
        
        if (writeHistFSI && (rank == MASTER_NODE)) {
          historyFile_FSI << setiosflags(ios::fixed) << setprecision(4) << CurrentTime << "," ;
          historyFile_FSI << setiosflags(ios::fixed) << setprecision(1) << iOuterIter << "," ;
          historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << WAitkDyn << "," ;
        }
        
      }
      
    }
    else {
      if (rank == MASTER_NODE) cout << "No relaxation method used. " << endl;
    }
  
  if (writeHistFSI && (rank == MASTER_NODE)) {historyFile_FSI.close();}
  
}

void CFEASolver::SetAitken_Relaxation(CGeometry **fea_geometry,
                                                 CConfig *fea_config, CSolver ***fea_solution) {
  
  unsigned long iPoint, iDim;
  unsigned short RelaxMethod_FSI;
  su2double *dispPred, *dispCalc;
  su2double WAitken;
  
  RelaxMethod_FSI = fea_config->GetRelaxation_Method_FSI();
  
  /*--- Only when there is movement it makes sense to update the solutions... ---*/
    
    if (RelaxMethod_FSI == NO_RELAXATION) {
      WAitken = 1.0;
    }
    else if (RelaxMethod_FSI == FIXED_PARAMETER) {
      WAitken = fea_config->GetAitkenStatRelax();
    }
    else if (RelaxMethod_FSI == AITKEN_DYNAMIC) {
      WAitken = GetWAitken_Dyn();
    }
    else {
      WAitken = 1.0;
    }
    
    // To nPointDomain; we need to communicate the solutions (predicted, old and old predicted) after this routine
    for (iPoint=0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Retrieve pointers to the predicted and calculated solutions ---*/
      dispPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
      dispCalc = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
      
      /*--- Set predicted solution as the old predicted solution ---*/
      fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Pred_Old();
      
      /*--- Set calculated solution as the old solution (needed for dynamic Aitken relaxation) ---*/
      fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Old(dispCalc);
      
      /*--- Apply the Aitken relaxation ---*/
      for (iDim=0; iDim < nDim; iDim++) {
        dispPred[iDim] = (1.0 - WAitken)*dispPred[iDim] + WAitken*dispCalc[iDim];
      }
      
    }
    
  
}

void CFEASolver::Update_StructSolution(CGeometry **fea_geometry,
                                                  CConfig *fea_config, CSolver ***fea_solution) {
  
  unsigned long iPoint;
  su2double *valSolutionPred;
  
  for (iPoint=0; iPoint < nPointDomain; iPoint++) {
    
    valSolutionPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
    
    fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution(valSolutionPred);
    
  }
  
  /*--- Perform the MPI communication of the solution, displacements only ---*/
  
  Set_MPI_Solution_DispOnly(fea_geometry[MESH_0], fea_config);
  
}


void CFEASolver::Compute_OFRefGeom(CGeometry *geometry, CSolver **solver_container, CConfig *config){

  unsigned short iVar;
  unsigned long iPoint;
  unsigned long nTotalPoint = 1;

  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);

  unsigned long ExtIter = config->GetExtIter();

  su2double reference_geometry = 0.0, current_solution = 0.0;

  bool fsi = config->GetFSI_Simulation();

  su2double objective_function = 0.0, objective_function_reduce = 0.0;
  su2double weight_OF = 1.0;

  su2double objective_function_averaged = 0.0;

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&nPointDomain,  &nTotalPoint,  1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
    nTotalPoint        = nPointDomain;
#endif

  weight_OF = config->GetRefGeom_Penalty() / nTotalPoint;

  for (iPoint = 0; iPoint < nPointDomain; iPoint++){

    for (iVar = 0; iVar < nVar; iVar++){

      /*--- Retrieve the value of the reference geometry ---*/
      reference_geometry = node[iPoint]->GetReference_Geometry(iVar);

      /*--- Retrieve the value of the current solution ---*/
      current_solution = node[iPoint]->GetSolution(iVar);

      /*--- The objective function is the sum of the difference between solution and difference, squared ---*/
      objective_function += weight_OF * (current_solution - reference_geometry)*(current_solution - reference_geometry);
    }

  }

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&objective_function,  &objective_function_reduce,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    objective_function_reduce        = objective_function;
#endif

  Total_OFRefGeom = objective_function_reduce + PenaltyValue;

  Global_OFRefGeom += Total_OFRefGeom;
  objective_function_averaged = Global_OFRefGeom / (ExtIter + 1.0 + EPS);

  bool direct_diff = ((config->GetDirectDiff() == D_YOUNG) ||
                      (config->GetDirectDiff() == D_POISSON) ||
                      (config->GetDirectDiff() == D_RHO) ||
                      (config->GetDirectDiff() == D_RHO_DL) ||
                      (config->GetDirectDiff() == D_EFIELD) ||
                      (config->GetDirectDiff() == D_MACH) ||
                      (config->GetDirectDiff() == D_PRESSURE));

  if ((direct_diff) && (rank == MASTER_NODE)){

    ofstream myfile_res;

    if (config->GetDirectDiff() == D_YOUNG) myfile_res.open ("Output_Direct_Diff_E.txt", ios::app);
    if (config->GetDirectDiff() == D_POISSON) myfile_res.open ("Output_Direct_Diff_Nu.txt", ios::app);
    if (config->GetDirectDiff() == D_RHO) myfile_res.open ("Output_Direct_Diff_Rho.txt", ios::app);
    if (config->GetDirectDiff() == D_RHO_DL) myfile_res.open ("Output_Direct_Diff_Rho_DL.txt", ios::app);
    if (config->GetDirectDiff() == D_EFIELD) myfile_res.open ("Output_Direct_Diff_EField.txt", ios::app);

    if (config->GetDirectDiff() == D_MACH) myfile_res.open ("Output_Direct_Diff_Mach.txt");
    if (config->GetDirectDiff() == D_PRESSURE) myfile_res.open ("Output_Direct_Diff_Pressure.txt");

    myfile_res.precision(15);

    myfile_res << scientific << Total_OFRefGeom << "\t";

    if (dynamic) myfile_res << scientific << objective_function_averaged << "\t";

    su2double local_forward_gradient = 0.0;
    su2double averaged_gradient = 0.0;

    local_forward_gradient = SU2_TYPE::GetDerivative(Total_OFRefGeom);

    if (fsi) {
      Total_ForwardGradient = local_forward_gradient;
      averaged_gradient     = Total_ForwardGradient / (ExtIter + 1.0);
    }
    else {
      Total_ForwardGradient += local_forward_gradient;
      averaged_gradient      = Total_ForwardGradient / (ExtIter + 1.0);
    }

    myfile_res << scientific << local_forward_gradient << "\t";

    myfile_res << scientific << averaged_gradient << "\t";

    myfile_res << endl;

    myfile_res.close();

    if (config->GetDirectDiff() == D_YOUNG)   cout << "Objective function: " << Total_OFRefGeom << ". Global derivative of the Young Modulus: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_POISSON) cout << "Objective function: " << Total_OFRefGeom << ". Global derivative of the Poisson's ratio: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_RHO)     cout << "Objective function: " << Total_OFRefGeom << ". Global derivative of the structural density: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_RHO_DL)  cout << "Objective function: " << Total_OFRefGeom << ". Global derivative of the dead weight: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_EFIELD)  cout << "Objective function: " << Total_OFRefGeom << ". Global derivative of the electric field: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_MACH)    cout << "Objective function: " << Total_OFRefGeom << ". Global derivative of the Mach number: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_PRESSURE) cout << "Objective function: " << Total_OFRefGeom << ". Global derivative of the freestream Pressure: " << Total_ForwardGradient << "." << endl;
  }
  else
  {

      // TODO: Temporary output file for the objective function. Will be integrated in the output once is refurbished.
   if (rank == MASTER_NODE){
      cout << "Objective function: " << Total_OFRefGeom << "." << endl;
      ofstream myfile_res;
      myfile_res.open ("of_refgeom.dat");
      myfile_res.precision(15);
      if (dynamic) myfile_res << scientific << objective_function_averaged << endl;
      else myfile_res << scientific << Total_OFRefGeom << endl;
      myfile_res.close();
      if (fsi){
          ofstream myfile_his;
          myfile_his.open ("history_refgeom.dat",ios::app);
          myfile_his.precision(15);
          myfile_his << scientific << Total_OFRefGeom << endl;
          myfile_his.close();
      }
    }

  }

}

void CFEASolver::Compute_OFRefNode(CGeometry *geometry, CSolver **solver_container, CConfig *config){

  unsigned short iVar;
  unsigned long iPoint;

  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);

  unsigned long ExtIter = config->GetExtIter();

  su2double reference_geometry = 0.0, current_solution = 0.0;

  bool fsi = config->GetFSI_Simulation();

  su2double objective_function = 0.0, objective_function_reduce = 0.0;
  su2double distance_sq = 0.0 ;
  su2double weight_OF = 1.0;

  su2double objective_function_averaged = 0.0;

  /*--- TEMPORARY, for application in dynamic TestCase ---*/
  su2double difX = 0.0, difX_reduce = 0.0;
  su2double difY = 0.0, difY_reduce = 0.0;

  weight_OF = config->GetRefNode_Penalty();

  for (iPoint = 0; iPoint < nPointDomain; iPoint++){

    if (geometry->node[iPoint]->GetGlobalIndex() == config->GetRefNode_ID() ){

      for (iVar = 0; iVar < nVar; iVar++){

        /*--- Retrieve the value of the reference geometry ---*/
        reference_geometry = config->GetRefNode_Displacement(iVar);

        /*--- Retrieve the value of the current solution ---*/
        current_solution = node[iPoint]->GetSolution(iVar);

        /*--- The objective function is the sum of the difference between solution and difference, squared ---*/
        distance_sq +=  (current_solution - reference_geometry)*(current_solution - reference_geometry);
      }

      objective_function = weight_OF * sqrt(distance_sq);

      difX = node[iPoint]->GetSolution(0) - config->GetRefNode_Displacement(0);
      difY = node[iPoint]->GetSolution(1) - config->GetRefNode_Displacement(1);

    }

  }

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&objective_function,  &objective_function_reduce,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&difX,  &difX_reduce,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&difY,  &difY_reduce,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    objective_function_reduce        = objective_function;
    difX_reduce                      = difX;
    difY_reduce                      = difY;
#endif

  Total_OFRefNode = objective_function_reduce + PenaltyValue;

  Global_OFRefNode += Total_OFRefNode;
  objective_function_averaged = Global_OFRefNode / (ExtIter + 1.0 + EPS);

  bool direct_diff = ((config->GetDirectDiff() == D_YOUNG) ||
                      (config->GetDirectDiff() == D_POISSON) ||
                      (config->GetDirectDiff() == D_RHO) ||
                      (config->GetDirectDiff() == D_RHO_DL) ||
                      (config->GetDirectDiff() == D_EFIELD) ||
                      (config->GetDirectDiff() == D_MACH) ||
                      (config->GetDirectDiff() == D_PRESSURE));

  if ((direct_diff) && (rank == MASTER_NODE)){

    ofstream myfile_res;

    if (config->GetDirectDiff() == D_YOUNG) myfile_res.open ("Output_Direct_Diff_E.txt", ios::app);
    if (config->GetDirectDiff() == D_POISSON) myfile_res.open ("Output_Direct_Diff_Nu.txt", ios::app);
    if (config->GetDirectDiff() == D_RHO) myfile_res.open ("Output_Direct_Diff_Rho.txt", ios::app);
    if (config->GetDirectDiff() == D_RHO_DL) myfile_res.open ("Output_Direct_Diff_Rho_DL.txt", ios::app);
    if (config->GetDirectDiff() == D_EFIELD) myfile_res.open ("Output_Direct_Diff_EField.txt", ios::app);

    if (config->GetDirectDiff() == D_MACH) myfile_res.open ("Output_Direct_Diff_Mach.txt");
    if (config->GetDirectDiff() == D_PRESSURE) myfile_res.open ("Output_Direct_Diff_Pressure.txt");

    myfile_res.precision(15);

    myfile_res << scientific << Total_OFRefNode << "\t";

    if (dynamic) myfile_res << scientific << objective_function_averaged << "\t";

    su2double local_forward_gradient = 0.0;
    su2double averaged_gradient = 0.0;

    local_forward_gradient = SU2_TYPE::GetDerivative(Total_OFRefNode);

    if (fsi) {
      Total_ForwardGradient = local_forward_gradient;
      averaged_gradient     = Total_ForwardGradient / (ExtIter + 1.0);
    }
    else {
      Total_ForwardGradient += local_forward_gradient;
      averaged_gradient      = Total_ForwardGradient / (ExtIter + 1.0);
    }

    myfile_res << scientific << local_forward_gradient << "\t";

    myfile_res << scientific << averaged_gradient << "\t";

    myfile_res << scientific << difX_reduce << "\t";
    myfile_res << scientific << difY_reduce << "\t";

    myfile_res << endl;

    myfile_res.close();

    if (config->GetDirectDiff() == D_YOUNG)   cout << "Objective function: " << Total_OFRefNode << ". Global derivative of the Young Modulus: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_POISSON) cout << "Objective function: " << Total_OFRefNode << ". Global derivative of the Poisson's ratio: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_RHO)     cout << "Objective function: " << Total_OFRefNode << ". Global derivative of the structural density: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_RHO_DL)  cout << "Objective function: " << Total_OFRefNode << ". Global derivative of the dead weight: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_EFIELD)  cout << "Objective function: " << Total_OFRefNode << ". Global derivative of the electric field: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_MACH)    cout << "Objective function: " << Total_OFRefNode << ". Global derivative of the Mach number: " << Total_ForwardGradient << "." << endl;
    if (config->GetDirectDiff() == D_PRESSURE) cout << "Objective function: " << Total_OFRefNode << ". Global derivative of the freestream Pressure: " << Total_ForwardGradient << "." << endl;
  }
  else
  {

    // TODO: Temporary output file for the objective function. Will be integrated in 
    // the output once is refurbished.
    if (rank == MASTER_NODE){
      cout << "Objective function: " << Total_OFRefNode << "." << endl;
      ofstream myfile_res;
      myfile_res.open ("of_refnode.dat");
      myfile_res.precision(15);
      if (dynamic) myfile_res << scientific << objective_function_averaged << endl;
      else myfile_res << scientific << Total_OFRefNode << endl;
      myfile_res.close();

      ofstream myfile_his;
      myfile_his.open ("history_refnode.dat",ios::app);
      myfile_his.precision(15);
      myfile_his << ExtIter << "\t";
      myfile_his << scientific << Total_OFRefNode << "\t";
      myfile_his << scientific << objective_function_averaged << "\t";
      myfile_his << scientific << difX_reduce << "\t";
      myfile_his << scientific << difY_reduce << endl;
      myfile_his.close();

    }

  }


}

void CFEASolver::Compute_OFVolFrac(CGeometry *geometry, CSolver **solver_container, CConfig *config)
{
  /*--- Perform a volume average of the physical density of the elements for topology optimization ---*/

  unsigned long iElem, nElem = geometry->GetnElem();
  su2double total_volume = 0.0, integral = 0.0;

  for (iElem=0; iElem<nElem; ++iElem) {
    /*--- count only elements that belong to the partition ---*/
    if ( geometry->node[geometry->elem[iElem]->GetNode(0)]->GetDomain() ){
      su2double volume = geometry->elem[iElem]->GetVolume();
      total_volume += volume;
      integral += volume*element_properties[iElem]->GetPhysicalDensity();
    }
  }
  
#ifdef HAVE_MPI
  {
    su2double tmp;
    SU2_MPI::Allreduce(&total_volume,&tmp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    total_volume = tmp;
    SU2_MPI::Allreduce(&integral,&tmp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    integral = tmp;
  }
#endif

  Total_OFVolFrac = integral/total_volume;

  // TODO: Temporary output file for the objective function. Will be integrated in the output once is refurbished.
  if (rank == MASTER_NODE){
    cout << "Objective function: " << Total_OFVolFrac << "." << endl;
    ofstream myfile_res;
    myfile_res.open ("of_volfrac.dat");
    myfile_res.precision(15);
    myfile_res << scientific << Total_OFVolFrac << endl;
    myfile_res.close();
  }
}

void CFEASolver::Stiffness_Penalty(CGeometry *geometry, CSolver **solver, CNumerics **numerics, CConfig *config){

  unsigned long iElem;
  unsigned short iNode, iDim, nNodes = 0;

  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord, val_Sol;
  int EL_KIND = 0;

  su2double elementVolume, dvValue, ratio;
  su2double weightedValue = 0.0;
  su2double weightedValue_reduce = 0.0;
  su2double totalVolume = 0.0;
  su2double totalVolume_reduce = 0.0;

  /*--- Loops over the elements in the domain ---*/

  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL)    {nNodes = 4; EL_KIND = EL_QUAD;}

    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}

    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    for (iNode = 0; iNode < nNodes; iNode++) {
        indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
        for (iDim = 0; iDim < nDim; iDim++) {
            val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
            val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
            element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
            element_container[FEA_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
        }
    }

    // Avoid double-counting elements:
    // Only add the value if the first node is in the domain
    if (geometry->node[indexNode[0]]->GetDomain()){

        // Compute the area/volume of the element
        if (nDim == 2)
        	elementVolume = element_container[FEA_TERM][EL_KIND]->ComputeArea();
        else
        	elementVolume = element_container[FEA_TERM][EL_KIND]->ComputeVolume();

        // Compute the total volume
        totalVolume += elementVolume;

        // Retrieve the value of the design variable
        dvValue = numerics[FEA_TERM]->Get_DV_Val(element_properties[iElem]->GetDV());

        // Add the weighted sum of the value of the design variable
        weightedValue += dvValue * elementVolume;

    }

  }

// Reduce value across processors for parallelization

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&weightedValue,  &weightedValue_reduce,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&totalVolume,  &totalVolume_reduce,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    weightedValue_reduce        = weightedValue;
    totalVolume_reduce          = totalVolume;
#endif

    ratio = 1.0 - weightedValue_reduce/totalVolume_reduce;

    PenaltyValue = config->GetTotalDV_Penalty() * ratio * ratio;



}

void CFEASolver::ComputeResidual_Multizone(CGeometry *geometry, CConfig *config){

  unsigned short iVar;
  unsigned long iPoint;
  su2double residual;

  /*--- Set Residuals to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++){
      SetRes_BGS(iVar,0.0);
      SetRes_Max_BGS(iVar,0.0,0);
  }

  /*--- Set the residuals ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
          residual = node[iPoint]->GetSolution(iVar) - node[iPoint]->Get_BGSSolution_k(iVar);
          AddRes_BGS(iVar,residual*residual);
          AddRes_Max_BGS(iVar,fabs(residual),geometry->node[iPoint]->GetGlobalIndex(),geometry->node[iPoint]->GetCoord());
      }
  }

  SetResidual_BGS(geometry, config);

}


void CFEASolver::UpdateSolution_BGS(CGeometry *geometry, CConfig *config){

  unsigned long iPoint;

  /*--- To nPoint: The solution must be communicated beforehand ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++){

    node[iPoint]->Set_BGSSolution_k();

  }

}

void CFEASolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  unsigned short iVar, nSolVar;
  unsigned long index;

  ifstream restart_file;
  string restart_filename, filename, text_line;

  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry[MESH_0]->GetnZone();

  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);
  bool fluid_structure = config->GetFSI_Simulation();
  bool discrete_adjoint = config->GetDiscrete_Adjoint();

  if (dynamic) nSolVar = 3 * nVar;
  else nSolVar = nVar;

  su2double *Sol = new su2double[nSolVar];

  /*--- Skip coordinates ---*/

  unsigned short skipVars = geometry[MESH_0]->GetnDim();

  /*--- Restart the solution from file information ---*/

  filename = config->GetSolution_FEMFileName();

  /*--- If multizone, append zone name ---*/

  if (nZone > 1)
    filename = config->GetMultizone_FileName(filename, iZone);

  if (dynamic) {
    filename = config->GetUnsteady_FileName(filename, val_iter);
  }

  /*--- Read all lines in the restart file ---*/

  int counter = 0;
  long iPoint_Local; unsigned long iPoint_Global = 0; unsigned long iPoint_Global_Local = 0;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, filename);
  }

  /*--- Load data from the restart into correct containers. ---*/

  counter = 0;
  for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1] + skipVars;
      for (iVar = 0; iVar < nSolVar; iVar++) Sol[iVar] = Restart_Data[index+iVar];

      for (iVar = 0; iVar < nVar; iVar++) {
        node[iPoint_Local]->SetSolution(iVar, Sol[iVar]);
        if (dynamic) {
          node[iPoint_Local]->SetSolution_time_n(iVar, Sol[iVar]);
          node[iPoint_Local]->SetSolution_Vel(iVar, Sol[iVar+nVar]);
          node[iPoint_Local]->SetSolution_Vel_time_n(iVar, Sol[iVar+nVar]);
          node[iPoint_Local]->SetSolution_Accel(iVar, Sol[iVar+2*nVar]);
          node[iPoint_Local]->SetSolution_Accel_time_n(iVar, Sol[iVar+2*nVar]);
        }
        if (fluid_structure && !dynamic) {
          node[iPoint_Local]->SetSolution_Pred(iVar, Sol[iVar]);
          node[iPoint_Local]->SetSolution_Pred_Old(iVar, Sol[iVar]);
        }
        if (fluid_structure && discrete_adjoint){
          node[iPoint_Local]->SetSolution_Old(iVar, Sol[iVar]);
        }
      }
      iPoint_Global_Local++;

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }
    
  }

  /*--- Detect a wrong solution file ---*/

  if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }
#ifndef HAVE_MPI
  rbuf_NotMatching = sbuf_NotMatching;
#else
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (rbuf_NotMatching != 0) {
    SU2_MPI::Error(string("The solution file ") + filename + string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- MPI. If dynamic, we also need to communicate the old solution ---*/

  solver[MESH_0][FEA_SOL]->Set_MPI_Solution(geometry[MESH_0], config);
  if (dynamic) solver[MESH_0][FEA_SOL]->Set_MPI_Solution_Old(geometry[MESH_0], config);
  if (fluid_structure && !dynamic){
      solver[MESH_0][FEA_SOL]->Set_MPI_Solution_Pred(geometry[MESH_0], config);
      solver[MESH_0][FEA_SOL]->Set_MPI_Solution_Pred_Old(geometry[MESH_0], config);
  }

  delete [] Sol;

  /*--- Delete the class memory that is used to load the restart. ---*/

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;
  
}

void CFEASolver::RegisterVariables(CGeometry *geometry, CConfig *config, bool reset)
{
  /*--- Register the element density to get the derivatives required for
  material-based topology optimization, this is done here because element_properties
  is a member of CFEASolver only. ---*/
  if (!config->GetTopology_Optimization()) return;

  for (unsigned long iElem = 0; iElem < geometry->GetnElem(); iElem++)
    element_properties[iElem]->RegisterDensity();
}

void CFEASolver::ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config)
{
  /*--- Extract and output derivatives for topology optimization, this is done
  here because element_properties is a member of CFEASolver only and the output
  structure only supports nodal values (these are elemental). ---*/
  if (!config->GetTopology_Optimization()) return;
    
  unsigned long iElem,
                nElem = geometry->GetnElem(),
                nElemDomain = geometry->GetGlobal_nElemDomain();

  /*--- Allocate and initialize an array onto which the derivatives of every partition
  will be reduced, this is to output results in the correct order, it is not a very
  memory efficient solution... single precision is enough for output. ---*/
  float *send_buf = new float[nElemDomain], *rec_buf = NULL;
  for(iElem=0; iElem<nElemDomain; ++iElem) send_buf[iElem] = 0.0;
    
  for(iElem=0; iElem<nElem; ++iElem) {
    unsigned long iElem_global = geometry->elem[iElem]->GetGlobalIndex();
    send_buf[iElem_global] = SU2_TYPE::GetValue(element_properties[iElem]->GetAdjointDensity());
  }

#ifdef HAVE_MPI
  if (rank == MASTER_NODE) rec_buf = new float[nElemDomain];
  /*--- Need to use this version of Reduce instead of the wrapped one because we use float ---*/
  MPI_Reduce(send_buf,rec_buf,nElemDomain,MPI_FLOAT,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD);
#else
  rec_buf = send_buf;
#endif

  /*--- The master writes the file ---*/
  if (rank == MASTER_NODE) {
    string filename = config->GetTopology_Optim_FileName();
    ofstream file;
    file.open(filename.c_str());
    for(iElem=0; iElem<nElemDomain; ++iElem) file << rec_buf[iElem] << endl;
    file.close();
  }
  
  delete [] send_buf;
#ifdef HAVE_MPI
  if (rank == MASTER_NODE) delete [] rec_buf;
#endif

}

void CFEASolver::FilterElementDensities(CGeometry *geometry, CConfig *config)
{
  /*--- Apply a filter to the design densities of the elements to generate the
  physical densities which are the ones used to penalize their stiffness. ---*/
  
  unsigned short type;
  su2double param, radius;
  
  vector<pair<unsigned short,su2double> > kernels;
  vector<su2double> filter_radius;
  for (unsigned short iKernel=0; iKernel<config->GetTopology_Optim_Num_Kernels(); ++iKernel)
  {
    config->GetTopology_Optim_Kernel(iKernel,type,param,radius);
    kernels.push_back(make_pair(type,param));
    filter_radius.push_back(radius);
  }

  unsigned long iElem, nElem = geometry->GetnElem();
  
  su2double *design_rho = new su2double [nElem],
            *physical_rho = new su2double [nElem];
  
  for (iElem=0; iElem<nElem; ++iElem)
    design_rho[iElem] = element_properties[iElem]->GetDesignDensity();
  
  geometry->FilterValuesAtElementCG(filter_radius, kernels, design_rho, physical_rho);
  
  /*--- Apply projection ---*/
  config->GetTopology_Optim_Projection(type,param);
  switch (type) {
    case NO_PROJECTION: break;
    case HEAVISIDE_UP:
      for (iElem=0; iElem<nElem; ++iElem)
        physical_rho[iElem] = 1.0-exp(-param*physical_rho[iElem])+physical_rho[iElem]*exp(-param);
      break;
    case HEAVISIDE_DOWN:
      for (iElem=0; iElem<nElem; ++iElem)
        physical_rho[iElem] = exp(-param*(1.0-physical_rho[iElem]))-(1.0-physical_rho[iElem])*exp(-param);
      break;
    default:
      SU2_MPI::Error("Unknown type of projection function",CURRENT_FUNCTION);
  }
  
  for (iElem=0; iElem<nElem; ++iElem)
    element_properties[iElem]->SetPhysicalDensity(physical_rho[iElem]);
  
  delete [] design_rho;
  delete [] physical_rho;
}


/* ----------------------- CSD ---------------------------- */

CModalSolver::CModalSolver(void) : CSolver() {
  
    nElement                = 0;
    nDim                    = 0;
    nVar                    = 0;
    nMarker                 = 0;
    nPoint                  = 0;
	nModes	                = 0;
	omega 	                = NULL;
	damping	                = NULL;
    modalForceLast          = NULL;
    modalForce              = NULL;
    SolRest                 = NULL;
    generalizedDisplacement = NULL;
    generalizedVelocity     = NULL;
}

CModalSolver::CModalSolver(CGeometry *geometry, CConfig *config) : CSolver() {

    unsigned long iPoint;
    unsigned short iMode, iVar, iDim;
    su2double Uinf = 1.0;

    
//     bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.
    //   bool fsi = config->GetFSI_Simulation();
    
//     element_based = false;          // A priori we don't have an element-based input file (most of the applications will be like this)
//     element_based = false;          // A priori we don't have an element-based input file (most of the applications will be like this)

	omega 	                = NULL;
	damping	                = NULL;
    modalForceLast          = NULL;
    modalForce              = NULL;
    SolRest                 = NULL;
    generalizedDisplacement = NULL;
    generalizedVelocity     = NULL;

    blasFunctions   = NULL;
    Ass             = NULL;
    Bss             = NULL;
    EMatrixInv      = NULL;

    DMatrix                 = NULL;
    EMatrix                 = NULL;
    AMatrix                 = NULL;
    AsMatrix                = NULL;
    EMatrixInv              = NULL;


    Conv_Check[0] = 0; Conv_Check[1] = 0; Conv_Check[2] = 0;
    
    nElement    = geometry->GetnElem();
    nDim        = geometry->GetnDim();
    nEqn        = 2;
    nModes      = config->GetNumberOfModes();
    nVar        = nEqn*nModes;                 // n. of vars for state-space calculation
    nMarker     = geometry->GetnMarker();
    nPoint      = geometry->GetnPoint();
    omega       = new su2double[nModes];
    omega[0]    = 106.69842;
    cout<< "nModes here is "<< nModes<<endl;
    cout << "modal solver initialized:" << nDim << "\t" << nVar << "\t" << nPoint << endl;
    
    // CModalVariable contains displacements in cartesian coordinates (solution variable)
    // TODO: check if these need to be absolute positions or are added.
    node    = new CVariable*[nPoint];
    // restart solution, i.e. generalized displacements for now???
    SolRest = new su2double[nVar];

// <<<<<<< HEAD
//     if (config->GetDynamic_Analysis()){
//         cout<<"Dynamic Analysis \n";
//         if (config->GetDynamic_Method() == MODAL_HARMONIC_BALANCE)
//         {
//             unsigned short nInst = 3;
//             su2double HB_Period =  0.05891103435003335;
//             HB_Period /= config->GetTime_Ref();
//             su2double deltaT = HB_Period/(su2double)(nInst);
// //            su2double HB_Omega[3]={0,106.69842,-106.69842};
//             su2double HB_Const = 2.0*PI_NUMBER/(su2double)(nInst);
//             su2double HB_t[3]={HB_Const*0,HB_Const*1,HB_Const*2};
//             // EInv Matrix
//             Initialize_Transformation_Matrix(nModes, nInst);
//             // E Matrix
//             Initialize_Transformation_Matrix_Inv(nModes, nInst);
//             // A Matrix
// //            Initialize_A_Matrix(nModes, nInst);
//             // D Matrix
//             Initialize_HB_Operator(nModes, nInst);
//             // As Matrix
//             Initialize_As_Matrix(nModes,nInst);
//         }
//     }
// =======
    if (config->GetDynamic_Analysis()){
        cout<<"Dynamic Analysis \n";
        if (config->GetDynamic_Method() == MODAL_HARMONIC_BALANCE)
        {
            unsigned short nInst = 3, i, j, k;
            nModes = 1;
            nVar *= nInst;
            su2double HB_Period =  0.05891103435003335;
            HB_Period /= config->GetTime_Ref();
            su2double HB_Const = 2.0*PI_NUMBER/(su2double)(nInst);
            su2double deltaT = HB_Period/(su2double)(nInst);
            su2double omega[3] = {0.0, 106.69842, -106.69842};
            HB_Omega = new su2double[nInst*nEqn];

            HVector = new su2double[2*nModes*nInst];
            for (k =0; k < (nModes); k++)
                for (i =0; i < (nModes); i++)
                    for (j =0; j < (nInst); j++)
                        HVector[(k*nEqn)+(i*nModes)+j] = sin(HB_Const*j);



            HB_Omega[0] = 0; HB_Omega[1] = 106.69842; HB_Omega[2] = -106.69842;
            HB_Omega[3] = 0; HB_Omega[4] = 106.69842; HB_Omega[5] = -106.69842;
            su2double HB_t[3]={HB_Const*0,HB_Const*1,HB_Const*2};
            // EInv Matrix
            InitializeHBMatrices(nModes, nInst);

            // EqnOfMotion
            Initialize_EqnOfMotion(nModes, nInst, QSolVector);
        }
    }
// >>>>>>> acac050e3cb0912cefc312ff94a8c162d39c5666

    /*--- Initialize from zero everywhere. ---*/
    for (iVar = 0; iVar < nVar; iVar++) SolRest[iVar] = 0.0;
    
    /*--- Initialize the BGS residuals in FSI problems. ---*/
    if (config->GetMultizone_Residual()){

        FSI_Residual      = 0.0;
        RelaxCoeff        = 1.0;
        ForceCoeff        = 1.0;

        Residual_BGS      = new su2double[nVar];
        Residual_Max_BGS  = new su2double[nVar];
        for (iVar = 0; iVar < nVar; iVar++) Residual_BGS[iVar]  = 0.0;
        for (iVar = 0; iVar < nVar; iVar++) Residual_Max_BGS[iVar]  = 0.0;

        /*--- Define some structures for locating max residuals ---*/
        Point_Max_BGS       = new unsigned long[nVar];  
        for (iVar = 0; iVar < nVar; iVar++) Point_Max_BGS[iVar]  = 0;
        Point_Max_Coord_BGS = new su2double*[nVar];
        for (iVar = 0; iVar < nVar; iVar++) {
            Point_Max_Coord_BGS[iVar] = new su2double[nDim];
            for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord_BGS[iVar][iDim] = 0.0;
        }
    }
    else{
        ForceCoeff        = 1.0;
    }
  
    // read in mode shapes (TODO: rename function)
    ReadCSD_Mesh_Nastran(config,geometry);
    
    for (iMode=0; iMode < nModes; ++iMode) {
        cout << "2-Mode " << iMode+1 << " Frequency:\t" << " = " << omega[iMode] << endl;
    }
}

CModalSolver::~CModalSolver(void) {
  
    unsigned short iMode;    //, jVar, iMode;
//   unsigned long iElem;
  

    if (generalizedDisplacement != NULL) {
        for(iMode = 0; iMode < nModes; ++iMode) delete generalizedDisplacement[iMode];
        delete [] generalizedDisplacement;
    }
    if (generalizedVelocity != NULL) {
        for(iMode = 0; iMode < nModes; ++iMode) delete generalizedVelocity[iMode];
        delete [] generalizedVelocity;
    }
    
    if (modalForceLast != NULL) delete [] modalForceLast;
    if (modalForce != NULL) delete [] modalForce;
    if (SolRest != NULL) delete [] SolRest;
    if (damping != NULL) delete [] damping;
    if (omega != NULL) delete [] omega;

    if( blasFunctions ) delete blasFunctions;
}

void CModalSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics, unsigned short iMesh, unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output) {
  
    /*
     * TODO:    adjoint (differentiation)
     *          nonlinear structure  
     *          body and other loads
     */
  
    unsigned long iPoint;
    cout << "modal solver Preprocessing completed\n";
//   unsigned short iMode;
//   bool initial_calc = (config->GetExtIter() == 0);  // Checks if it is the first calculation.
//   bool first_iter = (config->GetIntIter() == 0);    // Checks if it is the first iteration
// //   bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.
// //   bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
//   bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
//   bool newton_raphson = (config->GetKind_SpaceIteScheme_FEA() == NEWTON_RAPHSON);    // Newton-Raphson method
//   bool restart = config->GetRestart();                        // Restart analysis
//   bool initial_calc_restart = (SU2_TYPE::Int(config->GetExtIter()) == config->GetDyn_RestartIter()); // Initial calculation for restart
//   
// //   bool disc_adj_fem = (config->GetKind_Solver() == DISC_ADJ_FEM);   // Discrete adjoint FEM solver
//   bool incremental_load = config->GetIncrementalLoad();             // If an incremental load is applied
//   bool body_forces = config->GetDeadLoad();                         // Body forces (dead loads).
//   
//   bool fsi = config->GetFSI_Simulation();
//   bool consistent_interpolation = (!config->GetConservativeInterpolation() ||
//                                   (config->GetKindInterpolation() == WEIGHTED_AVERAGE));
//     
//     nModes  = geometry->getnModes();
//     omega   = new su2double[nModes];
//     damping = new su2double[nModes];
    //   nPointDomain  = geometry->GetnPointDomain();

}

void CModalSolver::ReadCSD_Mesh_Ansys(CConfig *config) {
    //TODO: add options to config file (hardwire test case for now?!)
    //
    unsigned long iPoint,iMode,nModesPoints;
    unsigned short nModes,iDim;
    vector<unsigned long>::iterator it;
    su2double frequency, Uinf;
    su2double Coord_3D[3];
    su2double *XV,*YV,*ZV;
    ifstream mode_file, mesh_file;
    string::size_type position;
    string line, dummy, text_line;
    char cstr[200];

    mesh_file.close();
    Uinf = 0.96*pow(config->GetGamma()*config->GetGas_Constant()*config->GetTemperature_FreeStream(),0.5);
    cout<< " Mach Inf       :: "<<config->GetMach_Motion()<<endl;
    cout<< " Gamma          :: "<<config->GetGamma()<<endl;
    cout<< " Temperature    :: "<<config->GetTemperature_FreeStream()<<endl;

    /* --- Read in modes' frequencies and mode shapes vectors --- */
    strcpy(cstr,"Mode_shaoe_node_coord_mode_1.dat");

    mode_file.open(cstr, ios::in);

    // open the input file
    if (mode_file.fail()) {
        cout << "Modes File could not be opened! Terminating!" << endl;
        exit(EXIT_FAILURE);
    }

    getline(mode_file, line);
    nModes = 1; //atoi(line.c_str());

//    getline(mode_file, line);
//    nModesPoints = atoi(line.c_str());
//    cout << "npoints/nModesPoints: " << nPoint << "\t" << nModesPoints << endl;
//    assert(nModesPoints == nPoint);

//    getline(mode_file, line);
//    refLength = atof(line.c_str());
//
//    cout << "N. of modes = " << nModes << "\n";
//    cout << "N. of nodes in CSD mesh = "<< nPoint << "\n";
//    cout << "Reference Length = " << refLength << endl;

    // initialize modes shapes, frequency and coordinate arrays
    modeShapes.reserve(nDim*nPoint*nModes);
    omega           = new su2double [nModes];
    modalForce      = new su2double [nModes];
    modalForceLast  = new su2double [nModes];
    XV              = new su2double [nPoint];
    YV              = new su2double [nPoint];
    ZV              = new su2double [nPoint];

    generalizedDisplacement   = new su2double *[nModes];
    generalizedVelocity       = new su2double *[nModes];

    for (iMode=0; iMode < nModes; ++iMode) generalizedDisplacement[iMode]   = new su2double [3];
    for (iMode=0; iMode < nModes; ++iMode) generalizedVelocity[iMode]       = new su2double [3];

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
        node[iPoint] = new CModalVariable(SolRest, nDim, nVar, nModes, config);
    }

    for (iMode=0; iMode < nModes; ++iMode) {
        getline(mode_file, line);
        istringstream iss(line);
        iss >> frequency;
        omega[iMode] = frequency;//2.0*PI_NUMBER*frequency*refLength/Uinf; //frequency; // TODO: compute Uinf from config;
        cout<< "Omega :: "<<omega[iMode]<<" Freq ::"<<  frequency<<endl;
        iss.clear();
    }

    for (iMode=0; iMode < nModes; ++iMode) {
        for (iPoint = 0 ; iPoint < nPoint; iPoint++) {
            getline(mode_file, line);
            istringstream iss(line);
            iss >> Coord_3D[0]; iss >> Coord_3D[1]; iss >> Coord_3D[2];
            XV[iPoint] = Coord_3D[0]/refLength;
            YV[iPoint] = Coord_3D[1]/refLength;
            ZV[iPoint] = Coord_3D[2]/refLength;
            for(iDim=0; iDim < nDim; ++iDim) node[iPoint]->SetModeVector(iMode,iDim,Coord_3D[iDim]);
            iss.clear();
        }
        for (iPoint = 0 ; iPoint < nPoint; iPoint++) modeShapes.push_back(XV[iPoint]);
        for (iPoint = 0 ; iPoint < nPoint; iPoint++) modeShapes.push_back(YV[iPoint]);
        for (iPoint = 0 ; iPoint < nPoint; iPoint++) modeShapes.push_back(ZV[iPoint]);
    }

    for (iMode=0; iMode < nModes; ++iMode) {
        cout << "Mode " << iMode+1 << " Frequency:\t" << " = " << omega[iMode] << endl;
    }

    mode_file.close();

    delete [] XV;
    delete [] YV;
    delete [] ZV;
}

void CModalSolver::ReadCSD_Mesh_Nastran(CConfig *config,CGeometry *geometry){
    //TODO: add options to config file (hardwire test case for now?!)
    //
    unsigned long iPoint,iMode,nModesPoints;
	unsigned short number_of_modes,iDim;
	vector<unsigned long>::iterator it;
	su2double frequency;
	su2double Coord_3D[3];
    su2double *XV,*YV,*ZV;
	ifstream mode_file,mesh_file;
    string::size_type position;
	string line,dummy,text_line;
	char cstr[200];

    mesh_file.close();
	Uinf = config->GetMach()*pow(config->GetGamma()*config->GetGas_Constant()*config->GetTemperature_FreeStream(),0.5);
    Qinf = ONE2 * config->GetGamma()*config->GetPressure_FreeStream()*config->GetMach()*config->GetMach();
//     massrat = config->GetPressure_FreeStream()/(config->GetGas_Constant()*config->GetTemperature_FreeStream())*pow(refLength,5);
//     su2double flutter_index = config->GetAeroelastic_Flutter_Speed_Index();
	cout<< " Mach Inf       :: "<<config->GetMach()<<endl;
    cout<< " Gamma          :: "<<config->GetGamma()<<endl;
    cout<< " Temperature    :: "<<config->GetTemperature_FreeStream() << endl;
    cout<< " Pressure_inf   :: "<<config->GetPressure_FreeStream() << endl;
    cout<< " Density_inf    :: "<<config->GetDensity_FreeStream() << endl;
    cout<< " U_inf          :: "<< Uinf << endl;
    cout<< " Q              :: "<< Qinf <<endl;
    cout<< " Density        :: "<< config->GetPressure_FreeStream()/(config->GetGas_Constant()*config->GetTemperature_FreeStream())<<endl;
    /* --- Read in modes' frequencies and mode shapes vectors ---2935.059 */
    strcpy(cstr,"modesFile.dat");

	mode_file.open(cstr, ios::in);

	// open the input file
	if (mode_file.fail()) {
		cout << "Modes File could not be opened! Terminating!" << endl;
		exit(EXIT_FAILURE);
	}

	getline(mode_file, line);
	number_of_modes = atoi(line.c_str());
	nModes = number_of_modes;

	getline(mode_file, line);
	nModesPoints = atoi(line.c_str());
    cout << "npoints/nModesPoints: " << nPoint << "\t" << nModesPoints << endl;
    assert(nModesPoints == nPoint);
    
	getline(mode_file, line);
	refLength = atof(line.c_str());

	cout << "N. of modes = " << nModes << "\n";
    cout << "N. of nodes in CSD mesh = "<< nPoint << "\n";
    cout << "Reference Length = " << refLength << endl;

	// initialize modes shapes, frequency and coordinate arrays
	modeShapes.reserve(nDim*nPoint*nModes);
    omega           = new su2double [nModes];
    modalForce      = new su2double [nModes];
    modalForceLast  = new su2double [nModes];
    XV              = new su2double [nPoint];
    YV              = new su2double [nPoint];
    ZV              = new su2double [nPoint];
    
    generalizedDisplacement   = new su2double *[nModes];
    generalizedVelocity       = new su2double *[nModes];
    
    for (iMode=0; iMode < nModes; ++iMode) generalizedDisplacement[iMode]   = new su2double [3];
    for (iMode=0; iMode < nModes; ++iMode) generalizedVelocity[iMode]       = new su2double [3];

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
        node[iPoint] = new CModalVariable(SolRest, nDim, nVar, nModes, config);
    }
    
	for (iMode=0; iMode < nModes; ++iMode) {
		getline(mode_file, line);
		istringstream iss(line); 
		iss >> frequency;
		omega[iMode] = 2.0*PI_NUMBER*frequency/(Uinf); //frequency;
		cout<< "vel_inf: " << Uinf << "; Omega :: "<<omega[iMode]<<" Freq ::"<<  frequency<<endl;
		iss.clear();
	}

	for (iMode=0; iMode < nModes; ++iMode) {
		for (iPoint = 0 ; iPoint < nPoint; iPoint++) {
			getline(mode_file, line);
			istringstream iss(line); 
			iss >> Coord_3D[0]; iss >> Coord_3D[1]; iss >> Coord_3D[2];
            XV[iPoint] = Coord_3D[0]/refLength;
            YV[iPoint] = Coord_3D[1]/refLength;
            ZV[iPoint] = Coord_3D[2]/refLength;
            
            for(iDim=0; iDim < nDim; ++iDim) node[iPoint]->SetModeVector(iMode,iDim,Coord_3D[iDim]);
            iss.clear();    //WARNING: why is this not scaled with refLength?
		}
		for (iPoint = 0 ; iPoint < nPoint; iPoint++) modeShapes.push_back(XV[iPoint]);
        for (iPoint = 0 ; iPoint < nPoint; iPoint++) modeShapes.push_back(YV[iPoint]);
        for (iPoint = 0 ; iPoint < nPoint; iPoint++) modeShapes.push_back(ZV[iPoint]);
	}

	for (iMode=0; iMode < nModes; ++iMode) {
        cout << "Mode " << iMode+1 << " Frequency:\t" << " = " << omega[iMode] << endl;
    }

	mode_file.close();
    
    //--- scale CSD mesh ---//
    for (iPoint = 0 ; iPoint < geometry->GetnPoint(); iPoint++){
        for(iDim = 0; iDim < nDim ; ++iDim){
            Coord_3D[iDim] = geometry->node[iPoint]->GetCoord(iDim);
            geometry->node[iPoint]->SetCoord(iDim,Coord_3D[iDim]/refLength);
        }
    }

    // --- initialize state-space matrices ---//
    Initialize_StateSpace_Matrices(0);
    
//     su2double *Ctmp;
// //     Cmatrix = new su2double[nModes*nModes];
//     Cmatrix = (su2double *) malloc( nModes*nModes*sizeof(su2double) );
//     unsigned long i,j;
//     for ( j = 0; j < nModes; j++ )
//     {
//       for ( i = 0; i < nModes; i++ )
//       {
//         Cmatrix[i+j*nModes] = ( su2double ) ( 10 * ( i + 1 ) + ( j + 1 ) );
//       }
//     }
//     const su2double *Atmp = &(Ass[0]);
//     const su2double *Btmp = &(Bss[0]);
//     Ctmp = (su2double *) malloc( nModes*2*nModes*sizeof(su2double) );
// 
//     cout << "\n\n  Blas matrix multiplication C matrix\n\nCol:    ";    
//     // test
//     dgemm('n','n',2*nModes, nModes, 2*nModes, 1.0L, Ass, 2*nModes, Bss, 2*nModes, 0.L, Ctmp, 2*nModes );
//   
//     //  Write the header.
//     for (int j = 0; j < nModes; j++ )
//     {
//       cout << setw(7) << j << "       ";
//     }
//     cout << "\n";
//     cout << "  Row\n";
//     cout << "\n";
//     for ( int i = 0; i < 2*nModes; i++ )
//     {
//       for ( int j = 0; j < nModes; j++ ) cout << setw(12) << Ctmp[i+j*2*nModes] << "  ";
//       cout << "\n";
//     }
//     su2double *x, *y;
//     x = new su2double[2*nModes];
//     y = new su2double[2*nModes];
//     for (int i = 0; i < 2*nModes; i++ ) x[i] = ( su2double ) ( i + 1 );
//   
//     bool trans = false;
//     dgemv( trans, 2*nModes,2*nModes, 1.0, Ass,2*nModes,x,1,0.,y, 1 );
//     for ( int j = 0; j < 2*nModes; j++ ) cout << setw(5) << x[j] << " ";
//     cout << endl;
//     for ( int j = 0; j < 2*nModes; j++ ) cout << setw(5) << y[j] << " ";
//     cout << endl;
//     exit(0);
    
    for (iMode=0; iMode < nModes; ++iMode) cout << "X0=" << generalizedDisplacement[iMode] << endl;
    for (iMode=0; iMode < nModes; ++iMode) cout << "V0=" << generalizedVelocity[iMode] << endl;
    //--- SPM ---//
    delete [] XV;
    delete [] YV;
    delete [] ZV;
}

void CModalSolver::RK2(CGeometry *geometry, CSolver **solver_container, CConfig *config){

    unsigned short iMode, iDim;
    su2double *qsol;
    su2double dy[4] = {0.0, 0.0, 0.0, 0.0};
    
    su2double dt = 0.015;//config->GetDelta_UnstTime();
    cout<<" Time Step is :: "<< dt <<endl;
    qsol = new su2double[2*nModes];

    cout << "solving structural equations of motion using two-stage RK method "<< endl;
    // solution array includes X, Y, Z displacements, obtained from gen. vars;
    // and must be parsed to node variable
    // generalizedXXX contains solution from state-space modal problem

    // Get the modal forces from the interpolation
    ComputeModalFluidForces(geometry, config);

    su2double DampingRatio = 0.0;

    
    for( iMode = 0; iMode < nModes; ++iMode) {
        dy[0] = generalizedVelocity[iMode][1];
        dy[1] = modalForceLast[iMode] - omega[iMode]*omega[iMode]*generalizedDisplacement[iMode][1] - DampingRatio*omega[iMode]*omega[iMode]*generalizedVelocity[iMode][1];

        dy[2] = generalizedVelocity[iMode][0];
        dy[3] = modalForce[iMode] - omega[iMode]*omega[iMode]*generalizedDisplacement[iMode][0] - DampingRatio*omega[iMode]*omega[iMode]*generalizedVelocity[iMode][0];
        
        qsol[2*iMode]   = generalizedDisplacement[iMode][1] + ONE2*dt*(dy[2] + dy[0]);
        qsol[2*iMode+1] = generalizedVelocity[iMode][1] + ONE2*dt*(dy[3] + dy[1]);
    }

    // Update old and new solution
    for( iMode = 0; iMode < nModes; ++iMode) {
        generalizedDisplacement[iMode][1]   = generalizedDisplacement[iMode][0];
        generalizedVelocity[iMode][1]       = generalizedVelocity[iMode][0];
        generalizedDisplacement[iMode][0]   = qsol[2*iMode];
        generalizedVelocity[iMode][0]       = qsol[2*iMode+1];
    }

    UpdateStructuralNodes();

    delete [] qsol;
}

void CModalSolver::HB_RK4(CGeometry *geometry, CSolver **solver_container, CConfig *config){

    su2double *k1, *k2, *k3, *k4, *QSol_RK;
    su2double coeff1 = 1.0/6.0, coeff2 = 2.0/6.0, coeff3 = 2.0/6.0, coeff4 = 1.0/6.0;
    su2double dt = 0.005;
    unsigned short i, j; // iMode;
    unsigned short nVars = (nModes*nInst*nEqn);
    su2double flutter_index = config->GetAeroelastic_Flutter_Speed_Index();
    
    QSol_RK = new su2double[nVars];

    ComputeModalFluidForces(geometry, config);

    k1 = new su2double[nVars];
    k2 = new su2double[nVars];
    k3 = new su2double[nVars];
    k4 = new su2double[nVars];
    for (i=0; i<nVars; i++){
        // RK Stage 1
        k1[i] = dt * SystemVector[i] ;
        k2[i] = 0.0 ;
        k3[i] = 0.0 ;
        k4[i] = 0.0 ;
    }

    // RK Stage 2
    for(i=0; i<nVars; i++)
        QSol_RK[i] = QSolVector[i]+(k1[i]*0.5);
    Initialize_EqnOfMotion(nModes, nInst, QSol_RK);
    for (i=0; i<nVars; i++)
        k2[i] = dt * dt * 0.5 * SystemVector[i];

    // RK Stage 3
    for(i=0; i<nVars; i++)
        QSol_RK[i] = QSolVector[i]+(k2[i]*0.5);
    Initialize_EqnOfMotion(nModes, nInst, QSol_RK);
    for (i=0; i<nVars; i++)
        k2[i] = dt * dt * 0.5 * SystemVector[i];

    // RK Stage 4
    for(i=0; i<nVars; i++)
        QSol_RK[i] = QSolVector[i]+(k3[i]*0.5);
    Initialize_EqnOfMotion(nModes, nInst, QSol_RK);
    for (i=0; i<nVars; i++)
        k2[i] = dt * dt *  SystemVector[i];

    for (i = 0; i< nVars; i++){
        QSolVector[i] = QSolVector_Old[i] + (coeff1*k1[1]+coeff2*k1[2]+coeff3*k1[3]+coeff4*k1[4]);
    }

//     for( iMode = 0; iMode < nModes; ++iMode) {
//         generalizedDisplacement[iMode][1]   = generalizedDisplacement[iMode][0];
//         generalizedVelocity[iMode][1]       = generalizedVelocity[iMode][0];
//         generalizedDisplacement[iMode][0]   = qsol[2*iMode];
//         generalizedVelocity[iMode][0]       = qsol[2*iMode+1];
//     }

    UpdateStructuralNodes();


}

void CModalSolver::RK4(CGeometry *geometry, CSolver **solver_container, CConfig *config){

// <<<<<<< HEAD
    unsigned short iMode, iDim, irk, nStage,ind;
    su2double *qsol,*dy,*ForceVec;
    su2double *rk1,*rk2,*rk3,*rk4,*yout;
//     su2double dy[4] = {0.0, 0.0, 0.0, 0.0};
    su2double rkcoeff[4] = {0.0, 0.0, 0.0, 0.0};
    bool trans = false;
    su2double dt = config->GetTime_Step();
    su2double flutter_index = config->GetAeroelastic_Flutter_Speed_Index();
    
    qsol    = new su2double[2*nModes];
    dy      = new su2double[2*nModes];
    yout    = new su2double[2*nModes];
    ForceVec= new su2double[2*nModes];

    // R-K number of stages and coefficients
    nStage = 4;
    switch(nStage){
        case 4:
            rkcoeff[0]  = rkcoeff[3] = 0.16666666666666666;
            rkcoeff[1]  = rkcoeff[2] = 0.33333333333333;
            rk1         = new su2double[2*nModes];
            rk2         = new su2double[2*nModes];
            rk3         = new su2double[2*nModes];
            rk4         = new su2double[2*nModes];
        case 2:
            rkcoeff[0] = 0.5;
            rkcoeff[1] = 0.5;
            rk1         = new su2double[2*nModes];
            rk2         = new su2double[2*nModes];
    }

    cout << "solving structural equations of motion using n-stage RK method "<< endl;

    // solution array includes X, Y, Z displacements, obtained from gen. vars;
    // and must be parsed to node variable
    // generalizedXXX contains solution from state-space modal problem
    //
    //   d |y|   |  0      I   ||y|          |   0   |
    // ----| | = |             || | + Q/Uinf*|       |
    //  dt |ŷ|   | -K/M   -D/M ||ŷ|          | PHI*F |

    // Get the modal forces from the interpolation
    ComputeModalFluidForces(geometry, config);

    for(iMode = 0; iMode < nModes; ++iMode) qsol[iMode]         = generalizedDisplacement[iMode][0];
    for(iMode = 0; iMode < nModes; ++iMode) qsol[iMode+nModes]  = generalizedVelocity[iMode][0];
    for(iMode = 0; iMode < 2*nModes; ++iMode) dy[iMode]         = 0;
    for(iMode = 0; iMode < 2*nModes; ++iMode) yout[iMode]       = 0;

    for(iMode = 0; iMode < nModes; ++iMode) {
        ForceVec[iMode]             = 0;
        ForceVec[iMode+nModes]      = flutter_index*modalForce[iMode]*Bss[nModes+iMode + iMode*2*nModes]; //(modalForce[iMode]-modalForceLast[iMode])*Bss[nModes+iMode + iMode*2*nModes];
    }

    //rk1    
    dgemv(trans,2*nModes,2*nModes,1.0,Ass,2*nModes,qsol,1,0.,yout, 1 );
    for(iMode = 0; iMode < 2*nModes; ++iMode) rk1[iMode] = yout[iMode] + ForceVec[iMode];
    for(iMode = 0; iMode < 2*nModes; ++iMode) yout[iMode] = 0;
    
    //rk2
    for(iMode = 0; iMode < 2*nModes; ++iMode) dy[iMode] = qsol[iMode] + ONE2*dt*rk1[iMode];
    dgemv(false,2*nModes,2*nModes,1.0,Ass,2*nModes,dy,1,0.,yout, 1 );
    for(iMode = 0; iMode < 2*nModes; ++iMode) rk2[iMode] = yout[iMode] + ForceVec[iMode];
    for(iMode = 0; iMode < 2*nModes; ++iMode) yout[iMode] = 0;

    //rk3
    for(iMode = 0; iMode < 2*nModes; ++iMode) dy[iMode] = qsol[iMode] + ONE2*dt*rk2[iMode];
    dgemv(false,2*nModes,2*nModes,1.0,Ass,2*nModes,dy,1,0.,yout, 1 );
    for(iMode = 0; iMode < 2*nModes; ++iMode) rk3[iMode] = yout[iMode] + ForceVec[iMode];
    for(iMode = 0; iMode < 2*nModes; ++iMode) yout[iMode] = 0;
    
    //rk4
    for(iMode = 0; iMode < 2*nModes; ++iMode) dy[iMode] = qsol[iMode] + ONE2*dt*rk3[iMode];
    dgemv(false,2*nModes,2*nModes,1.0,Ass,2*nModes,dy,1,0.,yout, 1 );
    for(iMode = 0; iMode < 2*nModes; ++iMode) rk4[iMode] = yout[iMode] + ForceVec[iMode];    
    
    for( iMode = 0; iMode < nModes; ++iMode) {
        generalizedDisplacement[iMode][1]   = generalizedDisplacement[iMode][0];
        generalizedVelocity[iMode][1]       = generalizedVelocity[iMode][0];
//         
        generalizedDisplacement[iMode][0] += dt*(rkcoeff[0]*rk1[iMode] + 
        rkcoeff[1]*rk2[iMode] + rkcoeff[2]*rk3[iMode] + rkcoeff[3]*rk4[iMode]);
        ind = iMode+nModes;
        generalizedVelocity[iMode][0] += dt*(rkcoeff[0]*rk1[ind] + 
        rkcoeff[1]*rk2[ind] + rkcoeff[2]*rk3[ind] + rkcoeff[3]*rk4[ind]);
    }
    UpdateStructuralNodes();

//     for( iMode = 0; iMode < nModes; ++iMode) generalizedDisplacement[iMode][0] = 0;
//     for( iMode = 0; iMode < nModes; ++iMode) generalizedVelocity[iMode][0] = 0;
    
    delete [] qsol;
    delete [] rk1;
    delete [] rk2;
    delete [] rk3;
    delete [] rk4;
    delete [] dy;
    delete [] yout;
}

// void CModalSolver::RK4(CGeometry *geometry, CSolver **solver_container, CConfig *config){
// 
//     unsigned short iMode,iDim;
//     su2double *qsol;
//     su2double dy[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
// 
//     su2double dt = 0.025;
//     qsol = new su2double[2*nModes];
// 
//     cout << "solving structural equations of motion using two-stage RK method "<< endl;
//     // solution array includes X, Y, Z displacements, obtained from gen. vars;
//     // and must be parsed to node variable
//     // generalizedXXX contains solution from state-space modal problem
// 
//     // Get the modal forces from the interpolation
//     ComputeModalFluidForces(geometry, config);
// 
//     for( iMode = 0; iMode < nModes; ++iMode) {
// 
//         // RK Step 1
//         // K1  = dt * f (Tn, Yn)
//         dy[0] = dt*generalizedVelocity[iMode][0];
//         dy[1] = dt*(modalForceLast[iMode] - omega[iMode]*omega[iMode]*generalizedDisplacement[iMode][0]);
// 
//         // RK Step 2
//         // K2  = dt * f (Tn + dt/2, Yn + K1/2)
//         dy[2] = dt * (generalizedVelocity[iMode][0] + ONE2*dt*dy[1]) ;
//         dy[3] = dt * (modalForce[iMode] - omega[iMode]*omega[iMode]*(generalizedDisplacement[iMode][0]+ONE2*dt*dy[1])) ;
// 
//         // RK Step 3
//         // K3  = dt * f (Tn + dt/2, Yn + K2/2)
//         dy[4] = dt *(generalizedVelocity[iMode][0] + ONE2*dt*dy[3] );
//         dy[5] = dt *(modalForce[iMode] - omega[iMode]*omega[iMode]*(generalizedDisplacement[iMode][0]+ONE2*dt*dy[3]));
// 
//         // RK Step 4
//         // K4  = dt * f (Tn + dt, Yn + K3)
//         dy[6] = dt *(generalizedVelocity[iMode][0] + dt*dy[5] );
//         dy[7] = dt *(modalForce[iMode] - omega[iMode]*omega[iMode]*(generalizedDisplacement[iMode][0]+dt*dy[1]));
// 
//         // RK4 Solution
//         // Y N+1 = Y N + 1/6 * (K1 + 2K2 + 2K3 + K4)
//         qsol[2*iMode]   = generalizedDisplacement[iMode][0] + 0.16666*dt*(dy[0] + 2*dy[2] + 2*dy[4] + dy[6]);
//         qsol[2*iMode+1] = generalizedVelocity[iMode][0] + 0.16666*(dy[1] + 2*dy[3] + 2*dy[5] + dy[7]);
//     }
// 
//     // Update old and new solution
//     for( iMode = 0; iMode < nModes; ++iMode) {
//         generalizedDisplacement[iMode][1]   = generalizedDisplacement[iMode][0];
//         generalizedVelocity[iMode][1]       = generalizedVelocity[iMode][0];
//         generalizedDisplacement[iMode][0]   = qsol[2*iMode];
//         generalizedVelocity[iMode][0]       = qsol[2*iMode+1];
//     }
// 
//     UpdateStructuralNodes();
// 
//     delete [] qsol;
// }

void CModalSolver::UpdateStructuralNodes() {

    unsigned long iPoint;
    unsigned short iDim, iMode;
    su2double delta, solutionValue;

    // Initialization of the Solution variables
    for(iPoint = 0; iPoint < nPoint; ++iPoint) node[iPoint]->SetSolution_time_n();
    for(iPoint = 0; iPoint < nPoint; ++iPoint) node[iPoint]->SetSolution_Vel_time_n();
    for(iPoint = 0; iPoint < nPoint; ++iPoint)
        for(iDim = 0; iDim < nDim; ++iDim) node[iPoint]->SetSolution(iDim,0.0);
    for(iPoint = 0; iPoint < nPoint; ++iPoint)
        for(iDim = 0; iDim < nDim; ++iDim) node[iPoint]->SetSolution_Vel(iDim,0.0);

    // Loop over nodes to calculate the Gen. Disp. and Gen. Vel.
    for(iMode = 0; iMode < nModes; ++iMode) {
        delta = generalizedDisplacement[iMode][0] - generalizedDisplacement[iMode][1];
        for(iPoint = 0; iPoint < nPoint; ++iPoint) {
            for(iDim = 0; iDim < nDim; ++iDim) {
                solutionValue = delta*node[iPoint]->GetModeVector(iMode,iDim);
                node[iPoint]->Add_DeltaSolution(iDim,solutionValue);
            }
        }
        // repeat to update velocities
        delta = generalizedVelocity[iMode][0] - generalizedVelocity[iMode][1];
        for(iPoint = 0; iPoint < nPoint; ++iPoint) {
            for(iDim = 0; iDim < nDim; ++iDim) {
                solutionValue = delta*node[iPoint]->GetModeVector(iMode,iDim);
                node[iPoint]->Add_DeltaVelSolution(iDim,solutionValue);
            }
        }
    }
}

void CModalSolver::ComputeModalFluidForces(CGeometry *geometry, CConfig *config) {

  unsigned short iDim, iNode, nNode, iMode;
  unsigned long iPoint, iElem, nElem;

  unsigned short iMarkerInt, nMarkerInt = config->GetMarker_n_ZoneInterface()/2,
                 iMarker, nMarker = config->GetnMarker_All();

  /*--- Temporary storage to store the forces on the element faces ---*/
  vector<su2double> forces;
  vector<unsigned long> npointcheck;

  cout << "getting fluid force\n";
  /*--- Loop through the FSI interface pairs ---*/
  /*--- 1st pass to compute forces ---*/
  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; ++iMarkerInt) {
    /*--- Find the marker index associated with the pair ---*/
    for (iMarker = 0; iMarker < nMarker; ++iMarker)
      if (config->GetMarker_All_ZoneInterface(iMarker) == iMarkerInt)
        break;
    /*--- The current mpi rank may not have this marker ---*/
    if (iMarker == nMarker) continue;

    nElem = geometry->GetnElem_Bound(iMarker);

    for (iElem = 0; iElem < nElem; ++iElem) {
      /*--- Define the boundary element ---*/
      unsigned long nodes[4];
      su2double coords[4][3];
      bool quad = geometry->bound[iMarker][iElem]->GetVTK_Type() == QUADRILATERAL;
      nNode = quad? 4 : nDim;

      for (iNode = 0; iNode < nNode; ++iNode) {
        // getting pointer to stored node variable
        nodes[iNode] = geometry->bound[iMarker][iElem]->GetNode(iNode);
        // current location?
        for (iDim = 0; iDim < nDim; ++iDim) {
            coords[iNode][iDim] = geometry->node[nodes[iNode]]->GetCoord(iDim)+
                                node[nodes[iNode]]->GetSolution(iDim);
        }
      }

      /*--- Compute the area ---*/
      su2double area = 0.0;

      if (nDim == 2)
        area = (coords[0][0]-coords[1][0])*(coords[0][0]-coords[1][0])+
               (coords[0][1]-coords[1][1])*(coords[0][1]-coords[1][1]);

      if (nDim == 3) {
        su2double a[3], b[3], Ni, Nj, Nk;

        if (!quad) { // sides of the triangle
          for (iDim = 0; iDim < 3; iDim++) {
            a[iDim] = coords[1][iDim]-coords[0][iDim];
            b[iDim] = coords[2][iDim]-coords[0][iDim];
          }
        }
        else { // diagonals of the quadrilateral
          for (iDim = 0; iDim < 3; iDim++) {
            a[iDim] = coords[2][iDim]-coords[0][iDim];
            b[iDim] = coords[3][iDim]-coords[1][iDim];
          }
        }
        /*--- Area = 0.5*||a x b|| ---*/
        Ni = a[1]*b[2]-a[2]*b[1];
        Nj =-a[0]*b[2]+a[2]*b[0];
        Nk = a[0]*b[1]-a[1]*b[0];

        area = 0.25*(Ni*Ni+Nj*Nj+Nk*Nk);
      }
      area = sqrt(area);

      /*--- Integrate ---*/
      passivedouble weight = 1.0/nNode;
      su2double force[3] = {0.0, 0.0, 0.0};

    for (iNode = 0; iNode < nNode; ++iNode)
        for (iDim = 0; iDim < nDim; ++iDim) force[iDim] += weight*area*node[nodes[iNode]]->Get_FlowTraction(iDim);
        for (iDim = 0; iDim < nDim; ++iDim) {
          forces.push_back(force[iDim]);
// 	  cout << force[0] << "\t"  << force[1] << "\t" << force[2] << endl;
      }
    }
  }
  /*--- 2nd pass to set values. This is to account for overlap in the markers. ---*/
  /*--- By putting the integrated values back into the nodes no changes have to be made elsewhere. ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); ++iPoint) node[iPoint]->Clear_FlowTraction();
  
  vector<su2double>::iterator force_it = forces.begin();
  cout << "list points from elements\n";
  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; ++iMarkerInt) {
    /*--- Find the marker index associated with the pair ---*/
    for (iMarker = 0; iMarker < nMarker; ++iMarker)
      if (config->GetMarker_All_ZoneInterface(iMarker) == iMarkerInt) break;
    /*--- The current mpi rank may not have this marker ---*/
    if (iMarker == nMarker) continue;

    nElem = geometry->GetnElem_Bound(iMarker);

    for (iElem = 0; iElem < nElem; ++iElem) {
      bool quad = geometry->bound[iMarker][iElem]->GetVTK_Type() == QUADRILATERAL;
      nNode = quad? 4 : nDim;
      passivedouble weight = 1.0/nNode;

      su2double force[3];
      for (iDim = 0; iDim < nDim; ++iDim) force[iDim] = *(force_it++)*weight;

      for (iNode = 0; iNode < nNode; ++iNode) {
        iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
        node[iPoint]->Add_FlowTraction(force);
      }
    }
  }
#ifdef HAVE_MPI
  /*--- Perform a global reduction, every rank will get the nodal values of all halo elements ---*/
  /*--- This should be cheaper than the "normal" way, since very few points are both halo and interface ---*/
  vector<unsigned long> halo_point_loc, halo_point_glb;
  vector<su2double> halo_force;

  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; ++iMarkerInt) {
    /*--- Find the marker index associated with the pair ---*/
    for (iMarker = 0; iMarker < nMarker; ++iMarker)
      if (config->GetMarker_All_ZoneInterface(iMarker) == iMarkerInt)
        break;
    /*--- The current mpi rank may not have this marker ---*/
    if (iMarker == nMarker) continue;

    nElem = geometry->GetnElem_Bound(iMarker);

    for (iElem = 0; iElem < nElem; ++iElem) {
      bool quad = geometry->bound[iMarker][iElem]->GetVTK_Type() == QUADRILATERAL;
      nNode = quad? 4 : nDim;

      /*--- If this is an halo element we share the nodal forces ---*/
      for (iNode = 0; iNode < nNode; ++iNode)
        if (!geometry->node[geometry->bound[iMarker][iElem]->GetNode(iNode)]->GetDomain())
          break;

      if (iNode < nNode) {
        for (iNode = 0; iNode < nNode; ++iNode) {
          iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
          /*--- local is for when later we update the values in this rank ---*/
          halo_point_loc.push_back(iPoint);
          halo_point_glb.push_back(geometry->node[iPoint]->GetGlobalIndex());
          for (iDim = 0; iDim < nDim; ++iDim)
            halo_force.push_back(node[iPoint]->Get_FlowTraction(iDim));
        }
      }
    }
  }
  /*--- Determine the size of the arrays we need ---*/
  unsigned long nHaloLoc = halo_point_loc.size();
  unsigned long nHaloMax;
  MPI_Allreduce(&nHaloLoc,&nHaloMax,1,MPI_UNSIGNED_LONG,MPI_MAX,MPI_COMM_WORLD);

  /*--- Shared arrays, all the: number of halo points; halo point global indices; respective forces ---*/
  unsigned long *halo_point_num = new unsigned long[size];
  unsigned long *halo_point_all = new unsigned long[size*nHaloMax];
  su2double *halo_force_all = new su2double[size*nHaloMax*nDim];
  
  /*--- If necessary put dummy values in halo_point_glb to get a valid pointer ---*/
  if (halo_point_glb.empty()) halo_point_glb.resize(1);
  /*--- Pad halo_force to avoid (observed) issues in the adjoint when nHaloLoc!=nHaloMax ---*/
  while (halo_force.size() < nHaloMax*nDim) halo_force.push_back(0.0);
  
  MPI_Allgather(&nHaloLoc,1,MPI_UNSIGNED_LONG,halo_point_num,1,MPI_UNSIGNED_LONG,MPI_COMM_WORLD);
  MPI_Allgather(&halo_point_glb[0],nHaloLoc,MPI_UNSIGNED_LONG,halo_point_all,nHaloMax,MPI_UNSIGNED_LONG,MPI_COMM_WORLD);
  SU2_MPI::Allgather(&halo_force[0],nHaloMax*nDim,MPI_DOUBLE,halo_force_all,nHaloMax*nDim,MPI_DOUBLE,MPI_COMM_WORLD);

  /*--- Find shared points with other ranks and update our values ---*/
  for (int proc = 0; proc < size; ++proc)
  if (proc != rank) {
    unsigned long offset = proc*nHaloMax;
    for (iPoint = 0; iPoint < halo_point_num[proc]; ++iPoint) {
      unsigned long iPoint_glb = halo_point_all[offset+iPoint];
      ptrdiff_t pos = find(halo_point_glb.begin(),halo_point_glb.end(),iPoint_glb)-halo_point_glb.begin();
      if (pos < long(halo_point_glb.size())) {
        unsigned long iPoint_loc = halo_point_loc[pos];
        node[iPoint_loc]->Add_FlowTraction(&halo_force_all[(offset+iPoint)*nDim]);
      }
    }
  }

  delete [] halo_point_num;
  delete [] halo_point_all;
  delete [] halo_force_all;
#endif
  
  /* --- project force onto modes---*/

    for(iMode = 0; iMode < nModes; ++iMode){
        modalForceLast[iMode] = modalForce[iMode];
        modalForce[iMode] = 0.0;
    }

    su2double ScaleFactor = 1.225 *  pow(0.463,5);
        // compute new value
    for(iMode = 0; iMode < nModes; ++iMode){
        for(iPoint = 0; iPoint < nPoint; ++iPoint){
             for(iDim = 0; iDim < nDim; ++iDim){      
                modalForce[iMode] += node[iPoint]->GetModeVector(iMode,iDim)*node[iPoint]->Get_FlowTraction(iDim);
            }
        }
    }
}

void CModalSolver::ComputeModalFluidDamp(CGeometry *geometry, CConfig *config) {

  unsigned short iDim, iNode, nNode, iMode;
  unsigned long iPoint, iElem, nElem;

  unsigned short iMarkerInt, nMarkerInt = config->GetMarker_n_ZoneInterface()/2,
                 iMarker, nMarker = config->GetnMarker_All();

  /*--- Temporary storage to store the forces on the element faces ---*/
  vector<su2double> forces;
  vector<unsigned long> npointcheck;

  cout << "getting fluid damping\n";
  /*--- Loop through the FSI interface pairs ---*/
  /*--- 1st pass to compute forces ---*/
  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; ++iMarkerInt) {
    /*--- Find the marker index associated with the pair ---*/
    for (iMarker = 0; iMarker < nMarker; ++iMarker)
      if (config->GetMarker_All_ZoneInterface(iMarker) == iMarkerInt)
        break;
    /*--- The current mpi rank may not have this marker ---*/
    if (iMarker == nMarker) continue;

    nElem = geometry->GetnElem_Bound(iMarker);

    for (iElem = 0; iElem < nElem; ++iElem) {
      /*--- Define the boundary element ---*/
      unsigned long nodes[4];
      su2double coords[4][3];
      bool quad = geometry->bound[iMarker][iElem]->GetVTK_Type() == QUADRILATERAL;
      nNode = quad? 4 : nDim;

      for (iNode = 0; iNode < nNode; ++iNode) {
          //getting pointer to stored node variable
        nodes[iNode] = geometry->bound[iMarker][iElem]->GetNode(iNode);
        //current location?
        for (iDim = 0; iDim < nDim; ++iDim){
            coords[iNode][iDim] = geometry->node[nodes[iNode]]->GetCoord(iDim)+
                                node[nodes[iNode]]->GetSolution(iDim);
        }
      }

      /*--- Compute the area ---*/
      su2double area = 0.0;

      if (nDim == 2)
        area = (coords[0][0]-coords[1][0])*(coords[0][0]-coords[1][0])+
               (coords[0][1]-coords[1][1])*(coords[0][1]-coords[1][1]);

      if (nDim == 3) {
        su2double a[3], b[3], Ni, Nj, Nk;

        if (!quad) { // sides of the triangle
          for (iDim = 0; iDim < 3; iDim++) {
            a[iDim] = coords[1][iDim]-coords[0][iDim];
            b[iDim] = coords[2][iDim]-coords[0][iDim];
          }
        }
        else { // diagonals of the quadrilateral
          for (iDim = 0; iDim < 3; iDim++) {
            a[iDim] = coords[2][iDim]-coords[0][iDim];
            b[iDim] = coords[3][iDim]-coords[1][iDim];
          }
        }
        /*--- Area = 0.5*||a x b|| ---*/
        Ni = a[1]*b[2]-a[2]*b[1];
        Nj =-a[0]*b[2]+a[2]*b[0];
        Nk = a[0]*b[1]-a[1]*b[0];

        area = 0.25*(Ni*Ni+Nj*Nj+Nk*Nk);
      }
      area = sqrt(area);

      /*--- Integrate ---*/
      passivedouble weight = 1.0/nNode;
      su2double force[3] = {0.0, 0.0, 0.0};

    for (iNode = 0; iNode < nNode; ++iNode)
        for (iDim = 0; iDim < nDim; ++iDim) force[iDim] += weight*area*node[nodes[iNode]]->Get_FlowTraction(iDim);
        for (iDim = 0; iDim < nDim; ++iDim) {
          forces.push_back(force[iDim]);
      }
    }
  }
  /*--- 2nd pass to set values. This is to account for overlap in the markers. ---*/
  /*--- By putting the integrated values back into the nodes no changes have to be made elsewhere. ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); ++iPoint) node[iPoint]->Clear_FlowTraction();

  vector<su2double>::iterator force_it = forces.begin();
  cout << "list points from elements\n";
  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; ++iMarkerInt) {
    /*--- Find the marker index associated with the pair ---*/
    for (iMarker = 0; iMarker < nMarker; ++iMarker)
      if (config->GetMarker_All_ZoneInterface(iMarker) == iMarkerInt) break;
    /*--- The current mpi rank may not have this marker ---*/
    if (iMarker == nMarker) continue;

    nElem = geometry->GetnElem_Bound(iMarker);

    for (iElem = 0; iElem < nElem; ++iElem) {
      bool quad = geometry->bound[iMarker][iElem]->GetVTK_Type() == QUADRILATERAL;
      nNode = quad? 4 : nDim;
      passivedouble weight = 1.0/nNode;

      su2double force[3];
      for (iDim = 0; iDim < nDim; ++iDim) force[iDim] = *(force_it++)*weight;

      for (iNode = 0; iNode < nNode; ++iNode) {
        iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
        node[iPoint]->Add_FlowTraction(force);
      }
    }
  }
#ifdef HAVE_MPI
  /*--- Perform a global reduction, every rank will get the nodal values of all halo elements ---*/
  /*--- This should be cheaper than the "normal" way, since very few points are both halo and interface ---*/
  vector<unsigned long> halo_point_loc, halo_point_glb;
  vector<su2double> halo_force;

  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; ++iMarkerInt) {
    /*--- Find the marker index associated with the pair ---*/
    for (iMarker = 0; iMarker < nMarker; ++iMarker)
      if (config->GetMarker_All_ZoneInterface(iMarker) == iMarkerInt)
        break;
    /*--- The current mpi rank may not have this marker ---*/
    if (iMarker == nMarker) continue;

    nElem = geometry->GetnElem_Bound(iMarker);

    for (iElem = 0; iElem < nElem; ++iElem) {
      bool quad = geometry->bound[iMarker][iElem]->GetVTK_Type() == QUADRILATERAL;
      nNode = quad? 4 : nDim;

      /*--- If this is an halo element we share the nodal forces ---*/
      for (iNode = 0; iNode < nNode; ++iNode)
        if (!geometry->node[geometry->bound[iMarker][iElem]->GetNode(iNode)]->GetDomain())
          break;

      if (iNode < nNode) {
        for (iNode = 0; iNode < nNode; ++iNode) {
          iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
          /*--- local is for when later we update the values in this rank ---*/
          halo_point_loc.push_back(iPoint);
          halo_point_glb.push_back(geometry->node[iPoint]->GetGlobalIndex());
          for (iDim = 0; iDim < nDim; ++iDim)
            halo_force.push_back(node[iPoint]->Get_FlowTraction(iDim));
        }
      }
    }
  }
  /*--- Determine the size of the arrays we need ---*/
  unsigned long nHaloLoc = halo_point_loc.size();
  unsigned long nHaloMax;
  MPI_Allreduce(&nHaloLoc,&nHaloMax,1,MPI_UNSIGNED_LONG,MPI_MAX,MPI_COMM_WORLD);

  /*--- Shared arrays, all the: number of halo points; halo point global indices; respective forces ---*/
  unsigned long *halo_point_num = new unsigned long[size];
  unsigned long *halo_point_all = new unsigned long[size*nHaloMax];
  su2double *halo_force_all = new su2double[size*nHaloMax*nDim];

  /*--- If necessary put dummy values in halo_point_glb to get a valid pointer ---*/
  if (halo_point_glb.empty()) halo_point_glb.resize(1);
  /*--- Pad halo_force to avoid (observed) issues in the adjoint when nHaloLoc!=nHaloMax ---*/
  while (halo_force.size() < nHaloMax*nDim) halo_force.push_back(0.0);

  MPI_Allgather(&nHaloLoc,1,MPI_UNSIGNED_LONG,halo_point_num,1,MPI_UNSIGNED_LONG,MPI_COMM_WORLD);
  MPI_Allgather(&halo_point_glb[0],nHaloLoc,MPI_UNSIGNED_LONG,halo_point_all,nHaloMax,MPI_UNSIGNED_LONG,MPI_COMM_WORLD);
  SU2_MPI::Allgather(&halo_force[0],nHaloMax*nDim,MPI_DOUBLE,halo_force_all,nHaloMax*nDim,MPI_DOUBLE,MPI_COMM_WORLD);

  /*--- Find shared points with other ranks and update our values ---*/
  for (int proc = 0; proc < size; ++proc)
  if (proc != rank) {
    unsigned long offset = proc*nHaloMax;
    for (iPoint = 0; iPoint < halo_point_num[proc]; ++iPoint) {
      unsigned long iPoint_glb = halo_point_all[offset+iPoint];
      ptrdiff_t pos = find(halo_point_glb.begin(),halo_point_glb.end(),iPoint_glb)-halo_point_glb.begin();
      if (pos < long(halo_point_glb.size())) {
        unsigned long iPoint_loc = halo_point_loc[pos];
        node[iPoint_loc]->Add_FlowTraction(&halo_force_all[(offset+iPoint)*nDim]);
      }
    }
  }

  delete [] halo_point_num;
  delete [] halo_point_all;
  delete [] halo_force_all;
#endif

  /* --- project force onto modes---*/

    for(iMode = 0; iMode < nModes; ++iMode){
        modalForceLast[iMode] = modalForce[iMode];
        modalForce[iMode] = 0.0;
    }
    
    
    // compute new value
    for(iMode = 0; iMode < nModes; ++iMode){
        for( iPoint = 0; iPoint < nPoint; ++iPoint){
             for(iDim = 0; iDim < nDim; ++iDim){
                modalForce[iMode] += node[iPoint]->GetModeVector(iMode,iDim)*node[iPoint]->Get_FlowTraction(iDim);
            }
        }
    }

}

void CModalSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,  CNumerics **numerics,unsigned short iMesh) {
    
    su2double deltaMax,delta,dq,dqdot;
    unsigned short iDim,iMode;
    unsigned long iPoint;
    //sanity check: delta displacement fraction w.r.t. reference length
    
    deltaMax = 0;
    for( iPoint = 0; iPoint < nPoint; ++iPoint){
        delta = 0.;
        for(iDim = 0; iDim < nDim; ++iDim) delta += node[iPoint]->GetSolution(iDim);
        if(deltaMax < delta) deltaMax = delta;
    }
    cout << "max. displacement/refLength = " << deltaMax << "/" << refLength << endl;
    cout << "\tGen. Displacement\tGen. Velocity\n";
    for(iMode = 0; iMode < nModes; ++iMode){
        cout <<"Mode " << iMode+1 << "\t";
        dq      = generalizedDisplacement[iMode][0] - generalizedDisplacement[iMode][1];
        dqdot   = generalizedVelocity[iMode][0] - generalizedVelocity[iMode][1];
        cout << dq << "\t" << dqdot << endl;
    }
    
    return;    
}

void CModalSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {
  
  unsigned short iMode, i, j;
  bool incremental_load = config->GetIncrementalLoad();              // If an incremental load is applied
  
  //set initial velocities
  for(iMode = 0; iMode < nModes; ++iMode) generalizedVelocity[iMode][0] = 1e-2;

  if (config->GetDynamic_Method()==MODAL_HARMONIC_BALANCE){
      QSolVector = new su2double[2*nModes*nInst];
      QSolVector_Old = new su2double[2*nModes*nInst];
      for (i =0; i < (nModes*nInst); i++) {
          QSolVector[i] = 0.0;
          QSolVector[nModes * nInst + i] = 1e-2;
          QSolVector_Old[i] = 0.0;
          QSolVector_Old[nModes * nInst + i] = 1e-2;
      }
  }

}

void CModalSolver::ComputeResidual_Multizone(CGeometry *geometry, CConfig *config){
    unsigned short iVar;
    unsigned long iMode;
    su2double residual;

    /*--- Set Residuals to zero ---*/

    for (iMode = 0; iMode < 2*nModes; iMode++){
        SetRes_BGS(iMode,0.0);
        SetRes_Max_BGS(iMode,0.0,0);
    }

    /*--- Set the residuals ---*/ 
    /*--- compute residual for each modal variable ---*/

    // TODO:: AddRes_Max_BGS is just a hack at the moment more accurate implementation should come in soon.
    for (iMode = 0; iMode < nModes; iMode++){
            residual = generalizedDisplacement[iMode][0] - generalizedDisplacement[iMode][1];
            AddRes_BGS(2*iMode,residual*residual);
            AddRes_Max_BGS(2*iMode,fabs(residual),2*iMode,geometry->node[0]->GetCoord());

            residual = generalizedVelocity[iMode][0] - generalizedVelocity[iMode][1];
            AddRes_BGS(2*iMode+1,residual*residual);
            AddRes_Max_BGS(2*iMode+1,fabs(residual),2*iMode+1,geometry->node[0]->GetCoord());
    }

    SetResidual_BGS(geometry, config);
}

void CModalSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config){
    unsigned short iVar, iMarker, MarkerS, MarkerR;
    unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
    su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;

    bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.

    unsigned short nSolVar;

    nSolVar = nVar;

#ifdef HAVE_MPI
    int send_to, receive_from;
  SU2_MPI::Status status;
#endif

    for (iMarker = 0; iMarker < nMarker; iMarker++) {

        if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
            (config->GetMarker_All_SendRecv(iMarker) > 0)) {

            MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
            send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

            nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
            nBufferS_Vector = nVertexS*nSolVar;     nBufferR_Vector = nVertexR*nSolVar;

            /*--- Allocate Receive and send buffers  ---*/
            Buffer_Receive_U = new su2double [nBufferR_Vector];
            Buffer_Send_U = new su2double[nBufferS_Vector];

            /*--- Copy the solution that should be sent ---*/
            for (iVertex = 0; iVertex < nVertexS; iVertex++) {
                iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
                for (iVar = 0; iVar < nVar; iVar++)
                    Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
                if (dynamic) {
                    for (iVar = 0; iVar < nVar; iVar++) {
                        Buffer_Send_U[(iVar+nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Vel(iVar);
                        Buffer_Send_U[(iVar+2*nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Accel(iVar);
                    }
                }
            }

#ifdef HAVE_MPI

            /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

            /*--- Receive information without MPI ---*/
            for (iVertex = 0; iVertex < nVertexR; iVertex++) {
                for (iVar = 0; iVar < nVar; iVar++)
                    Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
                if (dynamic) {
                    for (iVar = nVar; iVar < 3*nVar; iVar++)
                        Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
                }
            }

#endif

            /*--- Deallocate send buffer ---*/
            delete [] Buffer_Send_U;

            /*--- Do the coordinate transformation ---*/
            for (iVertex = 0; iVertex < nVertexR; iVertex++) {

                /*--- Find point and its type of transformation ---*/
                iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();

                /*--- Copy solution variables. ---*/
                for (iVar = 0; iVar < nSolVar; iVar++)
                    SolRest[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];

                cout<<"SolRest[iVar] 5871 "<<SolRest[iVar]<<endl;
                /*--- Store received values back into the variable. ---*/
                for (iVar = 0; iVar < nVar; iVar++)
                    node[iPoint]->SetSolution(iVar, SolRest[iVar]);

                }

            }

            /*--- Deallocate receive buffer ---*/
            delete [] Buffer_Receive_U;

        }

    }

void CModalSolver::ImplicitNewmark_Update(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

    unsigned short iVar;
    unsigned long iPoint;
    unsigned long nPointDomain = geometry->GetnPointDomain();
    unsigned long nPoint = geometry->GetnPoint();
    su2double *valSolutionPred;
    su2double *valSolution;

    if (1) {
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

            valSolutionPred = node[iPoint]->GetSolution_Pred();
            valSolution = node[iPoint]->GetSolution();
            for (unsigned short iDim = 0; iDim<nDim; iDim++) {
                node[iPoint]->SetSolution_Pred_Old(iDim, valSolutionPred[iDim]);
                node[iPoint]->SetSolution_Pred(iDim, valSolution[iDim]);
            }

        }

        /*--- Perform the MPI communication of the solution ---*/

        Set_MPI_Solution(geometry, config);

        /*--- After the solution has been communicated, set the 'old' predicted solution as the solution ---*/
        /*--- Loop over n points (as we have already communicated everything ---*/

    }


}

void CModalSolver::Initialize_StateSpace_Matrices(unsigned short nInst) {
    //    Input, unsigned short  nModes, number of modes, leading dimension of the matrix = (2*nModes)x(2*nModes).
    //    number of columns implied from the state-space problem.
    //
    //    Output are matrix Ass and Bss.
    //    Ass is stored in (2*nModes)x(2*nModes) entries; Bss => is as an M x N matrix;
    //
    //              |   0       I   |                      | 0 |
    // Ass =        |               |;   BsMatrix = Q/Uinf*|   |
    //              |   K/M     D/M |                      | I |

//     su2double *a;
    unsigned long i, j;
    //m, n: number of rows and columns of the matrix.
    const unsigned short m = 2*nModes;
    const unsigned short n = 2*nModes;
    
    su2double DampingRatio = 0.0;

//     Ass = new su2double[m*n];
//     BsMatrix = new su2double[m*nModes];
    Ass = (su2double *) malloc( m*n*sizeof(su2double) );
    Bss = (su2double *) malloc( m*nModes*sizeof(su2double) );

    // initialize
    for ( i = 0; i < m*n; ++i ) Ass[i] = (su2double) (0);
    for ( i = 0; i < m*nModes; ++i ) Bss[i] = (su2double) (0);
    
    //Ass
    // I
    for ( j = nModes; j < n; j++ )   //loop over cols
        for ( i = 0; i < nModes; i++ ) if (i+nModes == j) Ass[i+j*m] = ( su2double ) (1); //rows
    
    //K/M
    for ( j = 0; j < nModes; j++ )   //loop over cols
        for ( i = nModes; i < m; i++ )
            if (i == j+nModes) Ass[i+j*m] = -omega[j]*omega[j]; //rows

    //D/M
    for ( j = nModes; j < n; j++ )   //loop over cols
        for ( i = nModes; i < m; i++ ) 
            if (i == j) Ass[i+j*m] = -2.*DampingRatio*omega[j-nModes]; // su2double (3); //rows

    //  Write the header.
    cout << "\nAs Matrix\n\n   Col:    ";
    for ( j = 0; j < n; j++ ) cout << setw(7) << j << "       ";
    
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
    
    for ( i = 0; i < n; i++ ){  //rows
      for ( j = 0; j < m; j++ ) cout << setw(12) << Ass[i+j*m] << "  ";     //cols
      cout << "\n";
    }

    
    //Bss
    // I
    for ( j = 0; j < nModes; j++ )   //loop over cols
    {
        for ( i = nModes; i < m; i++ ) {
            if (i == j+nModes) Bss[i+j*m] = Qinf/(Uinf*Uinf); //rows
        }
    }
    //  Write the header.
    cout << " \n Bs Matrix " << Uinf << "\t" << Qinf << "\n\nCol:    ";
    for ( j = 0; j < nModes; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
    for ( i = 0; i < m; i++ )
    {
      for ( j = 0; j < nModes; j++ ) cout << setw(12) << Bss[i+j*m] << "  ";
      cout << "\n";
    }
    
}

void CModalSolver::Initialize_HB_Operator(unsigned short nMode, unsigned short nInst) {

    unsigned short iInst, iMode, Di, Dj, Dk, nDij = nMode*nInst;
    su2double Const = 2.0 / nInst, factor;
    su2double iHB = 1;

    DMatrix = new su2double*[2*nDij] ;
    for (Di = 0; Di < 2*nDij; Di++)
        DMatrix[Di] = new su2double[2*nDij];

    for (Di = 0; Di < 2*nDij; Di++)
        for (Dj = 0; Dj<2*nDij; Dj++)
            DMatrix[Di][Dj] = 0.0;

    // Generalized displacement
    cout<< " +++ DMatrix +++ \n";
    for (Di = 0; Di < nDij; Di++) {
        for (Dj = 0; Dj < nDij; Dj++) {
            for (Dk = 0; Dk < iHB; Dk++) {
                factor = (2.0 * PI_NUMBER * (Dk + 1) * (Dj - Di)) / (nInst);
                DMatrix[Di][Dj] += (Dk + 1) * sin(factor);
            }
            DMatrix[Di][Dj] *= Const;
            cout <<"DMatrix("<<Di<<","<<Dj<<")="<<DMatrix[Di][Dj] << "\t";
        }
        cout << endl;
    }

    // Generalized Velocity
//    cout<< " +++ DMatrix +++ \n";
    for (Di = 0; Di < nDij; Di++) {
        for (Dj = 0; Dj < nDij; Dj++) {
            for (Dk = 0; Dk < iHB; Dk++) {
                factor = (2.0 * PI_NUMBER * (Dk + 1) * (Dj - Di)) / (nInst);
                DMatrix[nDij+Di][nDij+Dj] += (Dk + 1) * sin(factor);
            }
            DMatrix[nDij+Di][nDij+Dj] *= Const;
//            cout <<"DMatrix("<<Di<<","<<Dj<<")="<<DMatrix[Di][Dj] << "\t";
        }
//        cout << endl;
    }

    /*!
     * Square of D Matrix
     */
    // Initialize
    D2Matrix = new su2double* [2*nDij];

    for (Di =0; Di < 2*nDij; Di++)
        D2Matrix[Di] = new su2double[2*nDij];

    for (Di = 0; Di < 2*nDij; Di++)
        for (Dj = 0; Dj < 2*nDij; Dj++)
            D2Matrix[Di][Dj] = 0.0;

     // Matrix Multiplication

    for (Di = 0; Di < 2*nDij; Di++)
        for (Dj = 0; Dj < 2*nDij; Dj++)
            for (Dk = 0; Dk < 2*nDij; Dk++)
                D2Matrix[Di][Dj] += DMatrix[Di][Dk]*DMatrix[Dk][Dj];
}

void CModalSolver::Initialize_A_Matrix(unsigned short nMode, unsigned short nInst) {
    /*!
     *  A =
     */
    unsigned short iInst, iMode, Ai, Aj, nA=nMode*nInst;
    AMatrix = new su2double*[nMode] ;

    for (Ai = 0; Ai < nA; Ai++)
        EMatrix[Ai] = new su2double[nInst];

    for (Ai = 0; Ai < nA; Ai++)
        for (Aj = 0; Aj<nA; Aj++)
            AMatrix[Ai][Aj] = 0.0;
}

void CModalSolver::Initialize_Transformation_Matrix_Inv(unsigned short nMode, unsigned short nInst) {
    /*!         |    I  Icost0  Isint0  |
    *  EInv =   |    I  Icost1  Isint1  |
    *           |    I  Icost2  Isint2  |
    */
    unsigned short iInst, iMode, Ei, Ej, nE=nMode*nInst;

    su2double HB_Const = 2.0*PI_NUMBER/(su2double)(nInst);

    EMatrixInv = new su2double*[2*nE] ;
    for (Ei = 0; Ei < 2*nE; Ei++)
        EMatrixInv[Ei] = new su2double[2*nE];

    for (Ei = 0; Ei < 2*nE; Ei++)
        for (Ej = 0; Ej<2*nE; Ej++)
            EMatrixInv[Ei][Ej] = 0.0;
    // Generalized displacement
    for (Ei = 0; Ei < nInst; Ei++){
        for (Ej = 0; Ej < nMode; Ej++){
            EMatrixInv[Ej+(nMode*Ei)][Ej+(nMode*0)] = 1;
//            cout<<Ej+(nMode*Ei)<<" "<<Ej+(nMode*0)<<"\t";
            EMatrixInv[Ej+(nMode*Ei)][Ej+(nMode*1)] = cos(HB_Const*Ei);
//            cout<<Ej+(nMode*Ei)<<" "<<Ej+(nMode*1)<<"\t";
            EMatrixInv[Ej+(nMode*Ei)][Ej+(nMode*2)] = sin(HB_Const*Ei);
//            cout<<Ej+(nMode*Ei)<<" "<<Ej+(nMode*2)<<"\t";
        }
//        cout<< EMatrix[Ei][Ej] << "\t";
//        cout<<endl;
    }
    // Generalized Velocity
    for (Ei = 0; Ei < nInst; Ei++){
        for (Ej = 0; Ej < nMode; Ej++){
            EMatrixInv[nE+Ej+(nMode*Ei)][nE+Ej+(nMode*0)] = 1;
//            cout<<Ej+(nMode*Ei)<<" "<<Ej+(nMode*0)<<"\t";
            EMatrixInv[nE+Ej+(nMode*Ei)][nE+Ej+(nMode*1)] = cos(HB_Const*Ei);
//            cout<<Ej+(nMode*Ei)<<" "<<Ej+(nMode*1)<<"\t";
            EMatrixInv[nE+Ej+(nMode*Ei)][nE+Ej+(nMode*2)] = sin(HB_Const*Ei);
//            cout<<Ej+(nMode*Ei)<<" "<<Ej+(nMode*2)<<"\t";
        }
//        cout<< EMatrix[Ei][Ej] << "\t";
//        cout<<endl;
    }
//    cout<<endl;
    for (Ei = 0; Ei < 2*nE; Ei++) {
        for (Ej = 0; Ej < 2*nE; Ej++) {
            cout << EMatrixInv[Ei][Ej] << "\t";
        }
        cout<<endl;
    }
}

void CModalSolver::Initialize_Transformation_Matrix(unsigned short nMode, unsigned short nInst) {
    /*!       2      |    I       I       I       |
    *  EInv =  /     |    Icost0  Icost1  Icost2  |
    *           2N+1 |    Isint0  Isint1  Isint2  |
    */
    unsigned short iInst, iMode, Ei, Ej, nE=nMode*nInst;

    su2double HB_Const = 2.0*PI_NUMBER/(su2double)(nInst);

    su2double Const = 2.0 / nInst, factor;

    EMatrix = new su2double*[2*nE] ;

    for (Ei = 0; Ei < 2*nE; Ei++)
        EMatrix[Ei] = new su2double[2*nE];

    for (Ei = 0; Ei < 2*nE; Ei++)
        for (Ej = 0; Ej<2*nE; Ej++)
            EMatrix[Ei][Ej] = 0.0;

    // For Modal Displacements
    for (Ei = 0; Ei < nInst; Ei++){
        for (Ej = 0; Ej < nMode; Ej++){
//            cout<<Ej+(nMode*Ei)<<Ej+(nMode*0)<<"\t";
            EMatrix[Ej+(nMode*Ei)][Ej+(nMode*0)] = 0.5*Const;
//            cout<<Ej+(nMode*Ei)<<Ej+(nMode*1)<<"\t";
            EMatrix[Ej+(nMode*Ei)][Ej+(nMode*1)] = cos(HB_Const*Ei)*Const;
//            cout<<Ej+(nMode*Ei)<<Ej+(nMode*2)<<"\t";
            EMatrix[Ej+(nMode*Ei)][Ej+(nMode*2)] = sin(HB_Const*Ei)*Const;
        }
    }

    // For Modal Velocities
    for (Ei = 0; Ei < nInst; Ei++){
        for (Ej = 0; Ej < nMode; Ej++){
//            cout<<Ej+(nMode*Ei)<<Ej+(nMode*0)<<"\t";
            EMatrix[nE+Ej+(nMode*Ei)][nE+Ej+(nMode*0)] = 0.5*Const;
//            cout<<Ej+(nMode*Ei)<<Ej+(nMode*1)<<"\t";
            EMatrix[nE+Ej+(nMode*Ei)][nE+Ej+(nMode*1)] = cos(HB_Const*Ei)*Const;
//            cout<<Ej+(nMode*Ei)<<Ej+(nMode*2)<<"\t";
            EMatrix[nE+Ej+(nMode*Ei)][nE+Ej+(nMode*2)] = sin(HB_Const*Ei)*Const;
        }
    }

    cout<<"+++ InvMatrix +++\n";
    for (Ei = 0; Ei < 2*nE; Ei++) {
        for (Ej = 0; Ej < 2*nE; Ej++) {
            cout << EMatrix[Ei][Ej] << "\t";
        }
        cout<<endl;
    }
}

// void CModalSolver::Initialize_As_Matrix(unsigned short nMode, unsigned short nInst) {
//     /*!             |   0       I   |
//      * AsMatrix =   |               |
//      *              |   K/M     D/M |
//      */
//     unsigned short Asij = nMode*nInst*2;
//     unsigned short i,j;
// 
//     // Initialize Matrix
//     AsMatrix = new su2double*[Asij];
// 
//     for(i=0; i<Asij; i++)
//         AsMatrix[i] = new su2double[Asij];
// 
//     for(i=0; i<Asij; i++)
//         for(j=0; j<Asij; j++)
//             AsMatrix[i][j] = 0.0;
// 
//     // I Matrix
//     for(i = 0; i < Asij/2; i++)
//         AsMatrix[i][i+(Asij/2)] = -1.0;
// 
//     // K/M Matrix
//     for (i=0; i< Asij/2; i++)
//         AsMatrix[i+(Asij/2)][i] = 2.0; // K/M terms
// 
//     // D/M Matrix
//     for (i=0; i< Asij/2; i++)
//         AsMatrix[i+(Asij/2)][i+(Asij/2)] = 3.0; // D/M terms
// 
//     cout<<" +++ As Matrix +++ \n";
//     cout<<" Asij :: "<<Asij<<endl;
//     for (i=0; i < Asij; i++)
//         for (j=0; j < Asij; j++)
//             cout<<AsMatrix[i][j]<<"\t";
//         cout<<endl;
// }

void CModalSolver::dgemm( char transa, char transb, unsigned long m, unsigned long n, unsigned long k, 
  su2double alpha, su2double a[],  unsigned long lda, su2double b[],  unsigned long ldb, su2double beta, su2double c[],  unsigned long ldc ){

//****************************************************************************80
//
//  Purpose: DGEMM computes C = alpha * A * B and related operations.
//  Discussion: DGEMM performs one of the matrix-matrix operations
//     C := alpha * op ( A ) * op ( B ) + beta * C,
//    where op ( X ) is one of: op ( X ) = X   or   op ( X ) = X',
//
//    ALPHA and BETA are scalars, and A, B and C are matrices, with op ( A )
//    an M by K matrix, op ( B ) a K by N matrix and C an N by N matrix.
//  Licensing: This code is distributed under the GNU LGPL license.
//  Modified: 26 July 2019
//  Author: Original FORTRAN77 version by Jack Dongarra.
//    C++ version by John Burkardt.
//
//  Parameters:
//    Input, char TRANSA, specifies the form of op( A ) to be used in
//    the matrix multiplication as follows:
//    'N' or 'n', op ( A ) = A.
//    'T' or 't', op ( A ) = A'.
//    'C' or 'c', op ( A ) = A'.
//
//    Input, char TRANSB, specifies the form of op ( B ) to be used in
//    the matrix multiplication as follows:
//    'N' or 'n', op ( B ) = B.
//    'T' or 't', op ( B ) = B'.
//    'C' or 'c', op ( B ) = B'.
//
//    Input, int M, the number of rows of the  matrix op ( A ) and of the  
//    matrix C.  0 <= M.
//
//    Input, int N, the number  of columns of the matrix op ( B ) and the 
//    number of columns of the matrix C.  0 <= N.
//
//    Input, int K, the number of columns of the matrix op ( A ) and the 
//    number of rows of the matrix op ( B ).  0 <= K.
//
//    Input, double ALPHA, the scalar multiplier 
//    for op ( A ) * op ( B ).
//
//    Input, double A(LDA,KA), where:
//    if TRANSA is 'N' or 'n', KA is equal to K, and the leading M by K
//    part of the array contains A;
//    if TRANSA is not 'N' or 'n', then KA is equal to M, and the leading
//    K by M part of the array must contain the matrix A.
//
//    Input, int LDA, the first dimension of A as declared in the calling 
//    routine.  When TRANSA = 'N' or 'n' then LDA must be at least max ( 1, M ), 
//    otherwise LDA must be at least max ( 1, K ).
//
//    Input, double B(LDB,KB), where:
//    if TRANSB is 'N' or 'n', kB is N, and the leading K by N 
//    part of the array contains B;
//    if TRANSB is not 'N' or 'n', then KB is equal to K, and the leading
//    N by K part of the array must contain the matrix B.
//
//    Input, int LDB, the first dimension of B as declared in the calling 
//    routine.  When TRANSB = 'N' or 'n' then LDB must be at least max ( 1, K ), 
//    otherwise LDB must be at least max ( 1, N ).
//
//    Input, double BETA, the scalar multiplier for C.
//
//    Input/output, double C[LDC*N].
//    On input, the leading M by N part of this array must contain the 
//    matrix C, except when BETA is 0.0, in which case C need not be set.
//    On output, the array C is overwritten by the M by N matrix
//    alpha * op ( A ) * op ( B ) + beta * C.
//
//    Input, int LDC, the first dimension of C as declared in the calling 
//    routine.  max ( 1, M ) <= LDC.
//
  unsigned long i,j,l;  
  unsigned long ncola, nrowa, nrowb, oneL, m2;
  bool nota;
  bool notb;
  su2double temp;
  
  oneL = 1;
//
//  Set NOTA and NOTB as true if A and B respectively are not
//  transposed and set NROWA, NCOLA and NROWB as the number of rows
//  and columns of A and the number of rows of B respectively.
//
  nota = ( ( transa == 'N' ) || ( transa == 'n' ) );

  cout << "options0: "<< "\t" << transa << "\t" << transb << "\t" << nota << endl;
  if ( nota ){
    nrowa   = m;
    ncola   = k;
    m2      = m;
    cout << "options01" << endl;
  }
  else{
    nrowa   = k;
    ncola   = m;
    m2      = m;
    cout << "options02" << endl;
  }

  notb = ( ( transb == 'N' ) || ( transb == 'n' ) );
  cout << "options0: "<< "\t" << transb << "\t" << transb << "\t" << nota << endl;

  if ( notb ){
    nrowb = k;
  }
  else{
    nrowb = n;
  }
//
//  Test the input parameters.
//
  if ( ! ( transa == 'N' || transa == 'n' ||
           transa == 'C' || transa == 'c' ||
           transa == 'T' || transa == 't' ) )
  {
    cerr << "\n";
    cerr << "DGEMM - Fatal error!\n";
    cerr << "  Input TRANSA has an illegal value.\n";
    exit ( 1 );
  }

  if ( ! ( transb == 'N' || transb == 'n' ||
           transb == 'C' || transb == 'c' ||
           transb == 'T' || transb == 't' ) )
  {
    cerr << "\n";
    cerr << "DGEMM - Fatal error!\n";
    cerr << "  Input TRANSB has an illegal value.\n";
    exit ( 1 );
  }

  if ( m < 0 ){
    cerr << "\n";
    cerr << "DGEMM - Fatal error!\n";
    cerr << "  Input M has an illegal value.\n";
    exit ( 1 );
  }

  if ( n < 0 ){
    cerr << "\n";
    cerr << "DGEMM - Fatal error!\n";
    cerr << "  Input N has an illegal value.\n";
    exit ( 1 );
  }

  if ( k  < 0 ){
    cerr << "\n";
    cerr << "DGEMM - Fatal error!\n";
    cerr << "  Input K has an illegal value.\n";
    exit ( 1 );
  }

  if ( lda < max( oneL, nrowa ) ){
    cerr << "\n";
    cerr << "DGEMM - Fatal error!\n";
    cerr << "  Input LDA has an illegal value.\n";
    cout << lda << "\t" << oneL << "\t" << nrowa << endl;
    exit ( 1 );
  }

  if ( ldb < max( oneL, nrowb ) )
  {
    cerr << "\n";
    cerr << "DGEMM - Fatal error!\n";
    cerr << "  Input LDB has an illegal value.\n";
//     cout << ldb << "\t" << oneL << "\t" << nrowb << "\t" << n << "\t" << endl;
    exit ( 1 );
  }

  if ( ldc < max( oneL, m2 ) )
  {
    cerr << "\n";
    cerr << "DGEMM - Fatal error!\n";
    cerr << "  Input LDC has an illegal value.\n";
    exit ( 1 );
  }

//   cout << "cleared checks\n";
//
//  Quick return if possible.
//
  if ( m == 0 ) return;
  if ( n == 0 ) return;

  if ( ( alpha == 0.0 || k == 0 ) && ( beta == 1.0 ) ) return;
//
//  And if alpha is 0.0.
//
  if ( alpha == 0.0 ){
    for ( j = 0; j < n; j++ ) 
        for ( i = 0; i < m; i++ ) c[i+j*ldc] = 0.0;
    return;
  }
//
//  Start the operations.
//
  if ( notb ){
//
//  Form  C := alpha*A*B
//
    if ( nota ){
      for ( j = 0; j < n; j++ ){
        for ( i = 0; i < m; i++ ) c[i+j*ldc] = 0.0;
        for ( l = 0; l < k; l++ ){
            if ( b[l+j*ldb] != 0.0 ){
                temp = alpha * b[l+j*ldb];
                for ( i = 0; i < m; i++ ) c[i+j*ldc] = c[i+j*ldc] + temp * a[i+l*lda];
            }
        }
      }
    }
//
//  Form  C := alpha*A'*B
//
    else{
      for ( j = 0; j < n; j++ ){
        for ( i = 0; i < m; i++ ){
          temp = 0.0;
          for ( l = 0; l < k; l++ ) temp += a[l+i*lda] * b[l+j*ldb];
          c[i+j*ldc] = alpha * temp;
        }
      }
    }
  }
//
//  Form  C := alpha*A*B'
//
  else{
    if ( nota )
    {
      for ( j = 0; j < n; j++ )
      {
        for ( i = 0; i < m; i++ ) c[i+j*ldc] = 0.0;
        for ( l = 0; l < k; l++ ){
          if ( b[j+l*ldb] != 0.0 ){
            temp = alpha * b[j+l*ldb];
            for ( i = 0; i < m; i++ ) c[i+j*ldc] = c[i+j*ldc] + temp * a[i+l*lda];
          }
        }
      }
    }
//
//  Form  C := alpha*A'*B'
//
    else{
      for ( j = 0; j < n; j++ ){
        for ( i = 0; i < m; i++ ){
          temp = 0.0;
          for ( l = 0; l < k; l++ ) temp = temp + a[l+i*lda] * b[j+l*ldb];
          c[i+j*ldc] = alpha * temp;
        }
      }
    }
  }

  return;
}

//****************************************************************************80
void CModalSolver::dgemv( bool trans, unsigned long m, unsigned long n, double alpha, double a[], 
                          unsigned long lda, su2double x[], unsigned long incx, su2double beta, 
                          su2double y[], int incy )
//****************************************************************************80
//
//  Purpose: DGEMV computes y := alpha * A * x + beta * y for general matrix A.
//
//  Discussion: DGEMV performs one of the matrix-vector operations
//      y := alpha*A *x + beta*y
//    or
//      y := alpha*A'*x + beta*y,
//    where alpha and beta are scalars, x and y are vectors and A is an
//    m by n matrix.
//
//  Licensing: This code is distributed under the GNU LGPL license.
//  Modified: 02 April 2014 / 27 July 2019
//  Author: C++ version by John Burkardt
//
//  Parameters:
//    Input, char TRANS, specifies the operation to be performed:
//    'n' or 'N'   y := alpha*A *x + beta*y.
//    't' or 'T'   y := alpha*A'*x + beta*y.
//    'c' or 'C'   y := alpha*A'*x + beta*y.
//
//    - Input, int M, the number of rows of the matrix A. 0 <= M.
//    - Input, int N, the number of columns of the matrix A. 0 <= N.
//    - Input, double ALPHA, the scalar multiplier for A * x.
//    - Input, double A[LDA*N].  The M x N subarray contains the matrix A.
//    - Input, int LDA, the the first dimension of A as declared
//    in the calling routine.  max ( 1, M ) <= LDA.
//    - Input, double X[*], an array containing the vector to be 
//    multiplied by the matrix A.  
//    - If TRANS = 'N' or 'n', then X must contain N entries, stored in INCX 
//    increments in a space of at least ( 1 + ( N - 1 ) * abs ( INCX ) ) 
//    locations.
//    - Otherwise, X must contain M entries, store in INCX increments
//    in a space of at least ( 1 + ( M - 1 ) * abs ( INCX ) ) locations.
//
//    - Input, int INCX, the increment for the elements of X.  INCX must not be zero.
//    - Input, double BETA, the scalar multiplier for Y.
//
//    - Input/output, double Y[*], an array containing the vector to
//    be scaled and incremented by A*X.
//    - If TRANS = 'N' or 'n', then Y must contain M entries, stored in INCY
//    increments in a space of at least ( 1 + ( M - 1 ) * abs ( INCY ) ) 
//    locations.
//    - Otherwise, Y must contain N entries, store in INCY increments
//    in a space of at least ( 1 + ( N - 1 ) * abs ( INCY ) ) locations.
//
//    Input, int INCY, the increment for the elements of Y.  INCY must not be zero. --> TODO: remove feature
//
{
  unsigned long i,j,lenx,leny;
  unsigned short info;
  unsigned long oneL = 1;
  int ix, iy, jx, jy, kx, ky;
  su2double temp;
//
//  Test the input parameters.
//
  info = 0;

  if ( m < 0 ){
    info = 2;
  }
  else if ( n < 0 ){
    info = 3;
  }
  else if ( lda < max ( oneL, m ) )
  {
    info = 6;
  }
  else if ( incx == 0 )
  {
    info = 8;
  }
  else if ( incy == 0 )
  {
    info = 11;
  }

  if ( info != 0 ){
    cerr << "DGEMV: " << info << endl;
    return;
  }
//
//  Quick return if possible.
//
  if ( ( m == 0 ) ||
       ( n == 0 ) ||
       ( ( alpha == 0.0 ) && ( beta == 1.0 ) ) ) return;
//
//  Set LENX and LENY, the lengths of the vectors x and y, and set
//  up the start points in X and Y.
//
  if ( trans ){
    lenx = m;
    leny = n;
  }
  else{
    lenx = n;
    leny = m;
  }

  if ( 0 < incx ){
    kx = 0;
  }
  else{
    kx = 0 - ( lenx - 1 ) * incx;
  }

  if ( 0 < incy ){
    ky = 0;
  }
  else{
    ky = 0 - ( leny - 1 ) * incy;
  }
//
//  Start the operations. In this version the elements of A are
//  accessed sequentially with one pass through A.
//     for ( i = 0; i < m; i++ ) {cout << "y- " << y[i] << endl;}
  if ( alpha == 0.0 ) return;
//
//  Form y := alpha*A*x.
//
  if ( !(trans) ){
    jx = kx;
    if ( incy == 1 ){
      for ( j = 0; j < n; j++ ){
//         if ( abs(x[jx]) > 1e-14 ){
          temp = alpha * x[jx];
//           cout << alpha << "\t" << x[jx] << "\t" << temp << endl;
          for ( i = 0; i < m; i++ ) {
              y[i] = y[i] + temp * a[i+j*lda];
//               cout << "ind: " << i << "\t" << a[i+j*lda] << "\t" << y[i] << endl;
          }
//         }
        jx = jx + incx;
      }
    }
    else
    {
      for ( j = 0; j < n; j++ ){
        if ( abs(x[jx]) > 1e-14 ){
          temp = alpha * x[jx];
//         cout <<"v2 " << alpha << "\t" << x[jx] << "\t" << temp << endl;
          iy = ky;
          for ( i = 0; i < m; i++ )
          {
            y[iy] = y[iy] + temp * a[i+j*lda];
            iy = iy + incy;
//               cout << "v2: " << a[i+j*lda] << "\t" << y[i] << endl;
          }
        }
        jx = jx + incx;
      }
    }
  }
/*
  Form y := alpha*A'*x.
*/
  else{
    jy = ky;
    if ( incx == 1 )
    {
      for ( j = 0; j < n; j++ )
      {
        temp = 0.0;
        for ( i = 0; i < m; i++ )
        {
          temp = temp + a[i+j*lda] * x[i];
        }
        y[jy] = y[jy] + alpha * temp;
        jy = jy + incy;
      }
    }
    else
    {
      for ( j = 0; j < n; j++ )
      {
        temp = 0.0;
        ix = kx;
        for ( i = 0; i < m; i++ )
        {
          temp = temp + a[i+j*lda] * x[ix];
          ix = ix + incx;
        }
        y[jy] = y[jy] + alpha * temp;
        jy = jy + incy;
      }
    }
  }

  return;
}

void CModalSolver::Initialize_EqnOfMotion(unsigned short nMode, unsigned short nInst, su2double *QSol) {
    /*!
     *  Y = -As Q + B F
     */
    unsigned short iMode, iInst, i, j;
    su2double DampCoeff = 0.001;

    su2double *term1, *term2, **term3;
    term1 = new su2double[nMode*nInst*2];
    term2 = new su2double[nMode*nInst*2];
    term3 = new su2double*[nMode*nInst*2];
    SystemVector = new su2double[nMode*nInst*2];

    for (i=0; i<nMode*nInst*2;i++) {
        term1[i] = 0.0;
        term2[i] = 0.0;
        term3[i] = new su2double[nMode*nInst*2];
        SystemVector[i] = 0.0;
    }

    for (i = 0; i<nMode*nInst*2; i ++ )
        for (j = 0; j<nMode*nInst*2; j ++ ) {
            term3[i][j] = 0.0;
        }


    // -As Q
    for (i = 0; i<nMode*nInst*2; i ++ )
        for (j = 0; j<nMode*nInst*2; j ++ ) {
            term1[i] += -1.0*AsMatrix[i][j]*QSol[j];
        }

    // F H
    for (i=0; i<nMode*nInst*2; i++)
        for (i=0; i<nMode*nInst*2; i++)
            term2[i] += FMatrix[i][j]*HVector[j];

    for (i=0; i<nMode*nInst*2; i++)
            SystemVector[i] = term1[i]+term2[i];

}

void CModalSolver::InitializeHBMatrices(unsigned short nMode, unsigned short nInst){

    unsigned short i, j, nEqn=2;


    Initialize_Transformation_Matrix(nMode, nInst);

    // E Matrix
    Initialize_Transformation_Matrix_Inv(nMode, nInst);

    // A Matrix
    //            Initialize_A_Matrix(nModes, nInst);

    // D Matrix
    Initialize_HB_Operator(nMode, nInst);

    // As Matrix
//     Initialize_As_Matrix(nMode,nInst);

    FMatrix = new su2double*[2*nModes*nInst];

    for (i=0; i< (nEqn*nMode*nInst); i++){
        FMatrix[i] = new su2double[2 * nModes * nInst];
        for (j = 0; j < (nEqn * nMode * nInst); j++)
        FMatrix[i][j] = 0.0;
        FMatrix[i][i] = 1.0;
    }

}
