/*!
 * \file blas_structure.cpp
 * \brief Implementation of the functions that either simulate BLAS functionality
          or interface to an actual BLAS implementation.
 * \author E. van der Weide
 * \version 6.0.1 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#include "../include/blas_structure.hpp"
#include <cstring>

/* MKL or BLAS, if supported. */
#if (defined (HAVE_MKL) || defined(HAVE_BLAS)) && !(defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE))

/* Function prototypes for the BLAS routines used. */
extern "C" void dgemm_(char*, char*, const int*, const int*, const int*,
                       const passivedouble*, const passivedouble*, const int*,
                       const passivedouble*, const int*,
                       const passivedouble*, passivedouble*, const int*);

extern "C" void dgemv_(char*, const int*, const int*, const passivedouble*,
                       const passivedouble*, const int*, const passivedouble*,
                       const int*, const passivedouble*, passivedouble*, const int*);
#endif

/* Constructor. Initialize the const member variables, if needed. */
CBlasStructure::CBlasStructure(void)
#if !(defined(HAVE_LIBXSMM) || defined(HAVE_BLAS) || defined(HAVE_MKL)) || (defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE))
  : mc (256), kc (128), nc (128) 
#endif
{}

/* Destructor. Nothing to be done. */
CBlasStructure::~CBlasStructure(void) {}

/* Dense matrix multiplication, gemm functionality. */
void CBlasStructure::gemm(const int M,        const int N,        const int K,
                          const su2double *A, const su2double *B, su2double *C,
                          CConfig *config) {

  /* Initialize the variable for the timing, if profiling is active. */
#ifdef PROFILE
  double timeGemm;
  if( config ) config->GEMM_Tick(&timeGemm);
#endif

#if (defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE)) || !(defined(HAVE_LIBXSMM) || defined(HAVE_MKL) || defined(HAVE_BLAS))
  /* Native implementation of the matrix product. This optimized implementation
     assumes that the matrices are in column major order. This can be
     accomplished by swapping N and M and A and B. This implementation is based
     on https://github.com/flame/how-to-optimize-gemm. */
  for(int i=0; i < M*K; ++i)  cout << "A[" << i << "]= " << A[i] <<  "\t" << M << " " << N << " " << K << endl;
  for(int i=0; i < M*N; ++i)  cout << "B[" << i << "]= " << B[i] <<  endl;  
  for(int i=0; i < M*N; ++i)  cout << "C[" << i << "]= " << C[i] <<  endl;  
  cout << endl;
  gemm_imp(N, M, K, B, A, C);

#else
#ifdef HAVE_LIBXSMM

  /* The gemm function of libxsmm is used to carry out the multiplication.
     Note that libxsmm_gemm expects the matrices in column major order. That's
     why the in the calling sequence A and B and M and N are reversed. */
  su2double alpha = 1.0;
  su2double beta  = 0.0;
  char trans = 'N';

  libxsmm_dgemm(&trans, &trans, &N, &M, &K, &alpha, B, &N, A, &K, &beta, C, &N);

#else // MKL and BLAS

  /* The standard blas routine dgemm is used for the multiplication.
     Call dgemm without transposing the matrices. In that case dgemm expects
     the matrices in column major order, see the comments for libxsmm. */
  su2double alpha = 1.0;
  su2double beta  = 0.0;
  char trans = 'N';

  dgemm_(&trans, &trans, &N, &M, &K, &alpha, B, &N, A, &K, &beta, C, &N);

#endif
#endif

  /* Store the profiling information, if needed. */
#ifdef PROFILE
  if( config ) config->GEMM_Tock(timeGemm, M, N, K);
#endif
}

/* Dense matrix vector multiplication, gemv functionality. */
void CBlasStructure::gemv(const int M,        const int N,   const su2double *A,
                          const su2double *x, su2double *y) {

#if (defined (HAVE_BLAS) || defined(HAVE_MKL)) && !(defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE))

  /* The standard blas routine dgemv is used for the multiplication.
     Note that dgemv expects the matrices in column major order, while
     A is in row major order. This can be solved by using the transpose
     and switching M and N. */
     su2double alpha = 1.0;
     su2double beta  = 0.0;
     int       inc   = 1;
     char trans = 'T';

     dgemv_(&trans, &N, &M, &alpha, A, &N, x, &inc, &beta, y, &inc);

#else

  /* Native implementation of the matix vector product.
     Initialize the elements of y to zero. */
  memset(y, 0, M*sizeof(su2double));  

  /* Carry out the matrix vector product. */
  for(int k=0; k<M; ++k) {
    const su2double *AA = A + k*N;
    for(int l=0; l<N; ++l)
      y[k] += AA[l]*x[l];
  }

#endif
}

#if !(defined(HAVE_LIBXSMM) || defined(HAVE_BLAS) || defined(HAVE_MKL)) || (defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE))

/* Macros for accessing submatrices of a matmul using the leading dimension. */
#define A(i, j) a[(j)*lda + (i)]
#define B(i, j) b[(j)*ldb + (i)]
#define C(i, j) c[(j)*ldc + (i)]

/* Function, which perform the implementation of the gemm functionality.  */
void CBlasStructure::gemm_imp(const int m,        const int n,        const int k,
                              const su2double *a, const su2double *b, su2double *c) {

    
  /* Initialize the elements of c to zero. */
  memset(c, 0, m*n*sizeof(su2double));
  
  for(int i=0; i < m*n; ++i)  cout << "_imp B[" << i << "]= " << a[i] << endl;
  for(int i=0; i < n*k; ++i)  cout << "_imp A[" << i << "]= " << b[i] << endl;  
  for(int i=0; i < m*n; ++i)  cout << "_imp C[" << i << "]= " << c[i] << endl;  
  
  /* Set the leading dimensions of the three matrices. */
  const int lda = m;
  const int ldb = n;
  const int ldc = k;
  
  cout << "lds: " << lda << "\t" << ldb << "\t" << ldc << endl;

  const int kc = 128;
  const int mc = 256;
  const int nc = 128;

  /* The full matrix multiplication is split in several blocks.
     Loop over these blocks. */
  for(int p=0; p<k; p+=kc) {
    int pb = min(k-p, kc);
    for(int j=0; j<n; j+=nc) {
      int jb = min(n-j, nc);
      for(int i=0; i<m; i+=mc) {
        int ib = min(m-i, mc);
        cout << "indices: " << i << "  " << p << "  " << j << " -> " << A(i,p) << "\t" << B(p,j) << endl;
        /* Carry out the multiplication for this block. */
        gemm_inner(ib, jb, pb, &A(i, p), lda, &B(p, j), ldb, &C(i, j), ldc);
      }
    }
  }
}

/* Compute a portion of the c matrix one block at a time.
   Handle ragged edges with calls to a slow but general function. */
void CBlasStructure::gemm_inner(int m, int n, int k, const su2double *a, int lda,
                                const su2double *b, int ldb, su2double *c, int ldc) {

  /* Carry out the multiplication for this block. At the
     moment simply a call to gemm_arbitrary. */
  gemm_arbitrary(m, n, k, a, lda, b, ldb, c, ldc);
}

/* Naive gemm implementation to handle arbitrary sized matrices. */
void CBlasStructure::gemm_arbitrary(int m, int n, int k, const su2double *a, int lda,
                                    const su2double *b, int ldb, su2double *c, int ldc) {

  /* The order of these loops is tuned for column-major matrices. */
  for (int p = 0; p < k; p++) {
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < m; i++) {
        C(i, j) += A(i, p) * B(p, j);
      }
    }
  }
}

#undef C
#undef B
#undef A

#endif
