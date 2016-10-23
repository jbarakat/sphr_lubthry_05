/* LINEAR ALGEBRA
 *  Subroutines for manipulating matrices and vectors (especially matrices with banded structure).
 *  In some cases, BLAS and LAPACK subroutines are used.
 *
 * REFERENCES
 *  von Rosenberg, Methods for the Numerical Solution of Partial Differential Equations (Elsevier, 1969)
 *  Golub & Van Loan, Matrix Computations (The Johns Hopkins University Press, 1983)
 *  Anderson et al, LAPACK Users' Guide (Society for Industrial and Applied Mathematics, 1999)
 *  
 * PARAMETERS
 *  J		[input]			system size
 */

#ifndef LA_H
#define LA_H


/* HEADER FILES */
//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <iomanip>
//#include <string>
#include <cstdlib>
#include <cstdio>
#include <lapacke.h>
#include <cblas.h>
#include <vector>
#include <math.h>

using namespace std;

/* PROTOTYPES */
void laband  (int, double *, int, int, double *, int, int, double *);	
void laband  (int, double *, int, int, double *, double *);
void lablock (int, double *, double *, double *, double *, double *, double *, double *, double *);
void labipent(int, double *, double *, double *, double *, double *, double *, double *, double *);

/* IMPLEMENTATIONS */

/* banded matrix-matrix multiplication A*B = C
 *  A, B, C = JxJ matrices
 *  ua, la  = upper, lower bandwidth of A
 *  ub, lb  = upper, lower bandwidth of B
 */
void laband(int J, double *A, int ua, int la, 
                   double *B, int ub, int lb,
                   double *C){	
	int i, j, k, m;
	int ui, li, Ji;
	int uj, lj;
	int uc = min(J-1, ua + ub);
	int lc = min(J-1, la + lb);
	double c;

	// multiply A*B row-wise and store in C
	for (i = 0; i < J; i++){
		ui = min(J-1, i + ua);
		li = max(0  , i - la);
		uj = min(J-1, i + uc);
		lj = max(0  , i - lc);
		Ji = ui - li + 1;
		
		double a[Ji], b[Ji];
		for (m = 0; m < Ji; m++) a[m] = 0.0;
		for (m = 0; m < Ji; m++) b[m] = 0.0;
		for (k = li; k < ui+1; k++){
			m = k - li;
			a[m] = A[i*J + k];
		}

		for (j = lj; j < uj+1; j++){
			for (k = li; k < ui+1; k++){	
				m	= k - li;
				b[m] = B[k*J + j];
			}
			
			c = 0.0;	
			for (m = 0; m < Ji; m++){
				c += a[m]*b[m];
			}

			C[i*J + j] = c;
		}
	}
}

/* banded matrix-vector multiplication A*x = y
 *  A       = JxJ matrix
 *  x, y		= J-vectors
 *  ua, la  = upper, lower bandwidth of A
 */
void laband(int J, double *A, int ua, int la, 
                   double *x, double *y){	
	int i, j, k, m;
	int ui, li, Ji;
	double sum;
	int JJ = J*J;

	// multiply A*x row-wise and store in y
	for (i = 0; i < J; i++){
		ui = min(J-1, i + ua);
		li = max(0  , i - la);
		Ji = ui - li + 1;
		
		sum = 0.0;
		for (k = li; k < ui+1; k++){
			sum += A[i*J + k]*x[k];
		}

		y[i] = sum;
	}
}

/* solve block matrix-vector system by Gaussian elimination of the whole matrix
 * (employs LAPACK's DGETRF and DGETRS subroutines)
 *  
 *  [ A  B ][ u ]   [ fu ]
 *  [      ][   ] = [    ]
 *  [ C  D ][ v ]   [ fv ]
 
 *  where A, B, C, D are JxJ block  matrices,
 *  and u, v, f1, f2 are J-vectors
 *
 * NOTE: LAPACK also has subroutines for banded matrices,
 *       e.g. DGBTRF and DGBTRS
 */
void lablock(int J, double *A, double *B, double *C, double *D, 
             double *fu, double *fv, double *u, double *v){
	int i, j, ii, jj;
	int iblk, jblk;
	int J2 = J*2;
	int J2J2 = J2*J2;
	double mat[J2J2], rhs[J2], sln[J2];
	
	// assemble linear system mat*sln = rhs
	for (iblk = 0; iblk < 2; iblk++){
		for (jblk = 0; jblk < 2; jblk++){
			for (i = 0; i < J; i++){
				ii = i + J;
				for (j = 0; j < J; j++){
					jj = j + J;

					if (iblk == 0 && jblk == 0) mat[i *J2 + j ] = A[i*J + j];
					if (iblk == 0 && jblk == 1) mat[i *J2 + jj] = B[i*J + j];
					if (iblk == 1 && jblk == 0) mat[ii*J2 + j ] = C[i*J + j];
					if (iblk == 1 && jblk == 1) mat[ii*J2 + jj] = D[i*J + j];
				}
			}
		}
	}

	for (iblk = 0; iblk < 2; iblk++){
		for (i = 0; i < J; i++){
			ii = i + J;
			if (iblk == 0) rhs[i ] = fu[i];
			if (iblk == 1) rhs[ii] = fv[i];
		}
	}
	
	for (iblk = 0; iblk < 2; iblk++){
		for (i = 0; i < J; i++){
			ii = i + J;
			if (iblk == 0) sln[i ] = u[i];
			if (iblk == 1) sln[ii] = v[i];
		}
	}

	/* solve linear system by Gaussian elimination
	 * (LAPACK's DGETRF and DGETRS subroutines) */
	char trans = 'N';		// form of the system of equations
	int  ldmat =  J2;		// leading dimension of mat
	int  ldrhs =  1 ;		// leading dimension of rhs
	int  nrhs  =  1 ;		// number of right-hand sides
	int  ipiv[J2]   ;		// array of pivot indices 
	int  info       ;		// error flag

	// LU-factorize
	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, J2, J2, mat, ldmat, ipiv);
	if (info != 0){
		printf("Error: dgetrf did not return successfully.\n");
		return;
	}

	// solve
	info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, trans, J2, nrhs, mat, ldmat, ipiv, rhs, ldrhs);
	if (info != 0){
		printf("Error: dgetrs did not return successfully.\n");
		return;
	}

	// store solution
	for (i = 0; i < J2; i++){
		sln[i] = rhs[i];
	}

	// copy solution
	for (i = 0; i < J; i++){
		ii = i + J;
		u[i] = sln[i ];
		v[i] = sln[ii];
	}
}


/* bi-pentadiagonal algorithm
 *  
 *  [ A  B ][ u ]   [ fu ]
 *  [      ][   ] = [    ]
 *  [ C  D ][ v ]   [ fv ]
 *
 *  where A, B, C, D are JxJ block pentadiagonal matrices,
 *  and u, v, f1, f2 are J-vectors
 */
void labipent(int J, double *A, double *B, double *C, double *D, 
              double *fu, double *fv, double *u, double *v){
	int i, j, m, n;
	double a[J][4], b[J][4], c[J][4], d[J][4], e[J][4], f[J][2];
	double beta   [4], tau    [4], zeta   [2];
	double sig    [4], delt   [4], eta    [2];
	double mu;
	double lam [J][4], kapp[J][4], gamm[J][2];

	// initialize
	for (m = 0; m < 4; m++){
		beta   [m] = 0.0;
		tau    [m] = 0.0;
		sig    [m] = 0.0;
		delt   [m] = 0.0;
		for (i = 0; i < J; i++){
			a   [i][m] = 0.0;
			b   [i][m] = 0.0;
			c   [i][m] = 0.0;
			d   [i][m] = 0.0;
			e   [i][m] = 0.0;
			lam [i][m] = 0.0;
			kapp[i][m] = 0.0;
		}
	}
	for (n = 0; n < 2; n++){
		zeta   [n] = 0.0;
		eta    [n] = 0.0;
		for (i = 0; i < J; i++){
			f   [i][n] = 0.0;
			gamm[i][n] = 0.0;
		}
	}
	mu = 0.0;

	// store matrix coefficients
	for (i = 0; i < J; i++){
		if (i > 1){
			a[i][0] = A[i*J + (i-2)];
			a[i][1] = B[i*J + (i-2)];
			a[i][2] = C[i*J + (i-2)];
			a[i][3] = D[i*J + (i-2)];
		}
		
		if (i > 0){
			b[i][0] = A[i*J + (i-1)];
			b[i][1] = B[i*J + (i-1)];
			b[i][2] = C[i*J + (i-1)];
			b[i][3] = D[i*J + (i-1)];
		}
		
		c[i][0] = A[i*J +  i   ];
		c[i][1] = B[i*J +  i   ];
		c[i][2] = C[i*J +  i   ];
		c[i][3] = D[i*J +  i   ];
	
		if (i < J-1){
			d[i][0] = A[i*J + (i+1)];
			d[i][1] = B[i*J + (i+1)];
			d[i][2] = C[i*J + (i+1)];
			d[i][3] = D[i*J + (i+1)];
		}
		
		if (i < J-2){
			e[i][0] = A[i*J + (i+2)];
			e[i][1] = B[i*J + (i+2)];
			e[i][2] = C[i*J + (i+2)];
			e[i][3] = D[i*J + (i+2)];
		}
		
		f[i][0] = fu[i];
		f[i][1] = fv[i];
	}

	/*--- FORWARD ELIMINATION ---*/

	// initialize
	for (m = 0; m < 4; m++){
		sig    [m] = c[0][m];
		delt   [m] = d[0][m];
		beta   [m] = b[1][m];
		tau    [m] = c[1][m];
	}
	for (n = 0; n < 2; n++){
		eta    [n] = f[0][n];
		zeta   [n] = f[1][n];
	}

	// forward-compute coefficients (only need to store lam, kapp, gamm)
	for (i = 0; i < J; i++){
		if (i > 1){
			beta   [0] = b   [i][0] - a   [i][0]*lam [i-2][0] - a   [i][1]*lam [i-2][2]; 
			beta   [1] = b   [i][1] - a   [i][0]*lam [i-2][1] - a   [i][1]*lam [i-2][3]; 
			beta   [2] = b   [i][2] - a   [i][2]*lam [i-2][0] - a   [i][3]*lam [i-2][2]; 
			beta   [3] = b   [i][3] - a   [i][2]*lam [i-2][1] - a   [i][3]*lam [i-2][3]; 

			tau    [0] = c   [i][0] - a   [i][0]*kapp[i-2][0] - a   [i][1]*kapp[i-2][2]; 
			tau    [1] = c   [i][1] - a   [i][0]*kapp[i-2][1] - a   [i][1]*kapp[i-2][3]; 
			tau    [2] = c   [i][2] - a   [i][2]*kapp[i-2][0] - a   [i][3]*kapp[i-2][2]; 
			tau    [3] = c   [i][3] - a   [i][2]*kapp[i-2][1] - a   [i][3]*kapp[i-2][3]; 

			zeta   [0] = f   [i][0] - a   [i][0]*gamm[i-2][0] - a   [i][1]*gamm[i-2][1]; 
			zeta   [1] = f   [i][1] - a   [i][2]*gamm[i-2][0] - a   [i][3]*gamm[i-2][1]; 
		}
		if (i > 0){
			sig    [0] = tau    [0] - beta   [0]*lam [i-1][0] - beta   [1]*lam [i-1][2]; 
			sig    [1] = tau    [1] - beta   [0]*lam [i-1][1] - beta   [1]*lam [i-1][3]; 
			sig    [2] = tau    [2] - beta   [2]*lam [i-1][0] - beta   [3]*lam [i-1][2]; 
			sig    [3] = tau    [3] - beta   [2]*lam [i-1][1] - beta   [3]*lam [i-1][3]; 

			delt   [0] = d   [i][0] - beta   [0]*kapp[i-1][0] - beta   [1]*kapp[i-1][2]; 
			delt   [1] = d   [i][1] - beta   [0]*kapp[i-1][1] - beta   [1]*kapp[i-1][3]; 
			delt   [2] = d   [i][2] - beta   [2]*kapp[i-1][0] - beta   [3]*kapp[i-1][2]; 
			delt   [3] = d   [i][3] - beta   [2]*kapp[i-1][1] - beta   [3]*kapp[i-1][3]; 

			eta    [0] = zeta   [0] - beta   [0]*gamm[i-1][0] - beta   [1]*gamm[i-1][1]; 
			eta    [1] = zeta   [1] - beta   [2]*gamm[i-1][0] - beta   [3]*gamm[i-1][1]; 
		}
	
		mu         =  sig    [0]*sig    [3] - sig    [2]*sig    [1]    ;
		
		lam [i][0] = (sig    [3]*delt   [0] - sig    [1]*delt   [2])/mu;
		lam [i][1] = (sig    [3]*delt   [1] - sig    [1]*delt   [3])/mu;
		lam [i][2] = (sig    [0]*delt   [2] - sig    [2]*delt   [0])/mu;
		lam [i][3] = (sig    [0]*delt   [3] - sig    [2]*delt   [1])/mu;
		
		kapp[i][0] = (sig    [3]*e   [i][0] - sig    [1]*e   [i][2])/mu;
		kapp[i][1] = (sig    [3]*e   [i][1] - sig    [1]*e   [i][3])/mu;
		kapp[i][2] = (sig    [0]*e   [i][2] - sig    [2]*e   [i][0])/mu;
		kapp[i][3] = (sig    [0]*e   [i][3] - sig    [2]*e   [i][1])/mu;
		
		gamm[i][0] = (sig    [3]*eta    [0] - sig    [1]*eta    [1])/mu;
		gamm[i][1] = (sig    [0]*eta    [1] - sig    [2]*eta    [0])/mu;
	}


	/*--- BACKWARD SUBSTITUTION ---*/
	
	// initialize
	u[J-1] = gamm[J-1][0];
	v[J-1] = gamm[J-1][1];
	u[J-2] = gamm[J-2][0] - lam [J-2][0]*u[J-1] - lam [J-2][1]*v[J-1];
	v[J-2] = gamm[J-2][1] - lam [J-2][2]*u[J-1] - lam [J-2][3]*v[J-1];

	// backward-compute solution vector
	for (i = J-3; i > -1; i--){
		u[i] = gamm[i][0] - lam [i][0]*u[i+1] - lam [i][1]*v[i+1] 
		                  - kapp[i][0]*u[i+2] - kapp[i][1]*v[i+2];
		v[i] = gamm[i][1] - lam [i][2]*u[i+1] - lam [i][3]*v[i+1] 
		                  - kapp[i][2]*u[i+2] - kapp[i][3]*v[i+2];
	}
}

/* compute condition number of JxJ matrix A
 * (employs LAPACK's DGETRF and DGECON subroutines) 
 *
 * NOTE: LAPACK also has subroutines for banded matrices,
 *       e.g. DGBTRF and DGBCON
 */
void lacond(int J, double *A, double& rcond){
	int i, j;
	char   norm  = '1';		// specifies which norm is required
												// ('1','O' = 1-norm, 'I' = infinity=norm
	int    lda	 =  J ;		// leading dimension of mat
	int    ipiv[J]    ;		// array of pivot indices 
	double anorm      ;		// 1-norm or infinity-norm of A
	int    info       ;		// error flag

	// compute the norm of A 
	anorm = 0.0;
	if (norm == '1' || norm == 'O'){	// maximum absolute column sum
		for (j = 0; j < J; j++){
			double sum = 0.0;
			for (i = 0; i < J; i++){
				sum	+= fabs(A[i*J + j]);
			}
			anorm = max(anorm,sum);
		}
	}
	else if (norm == 'I'){						// maximum absolute row sum
		for (i = 0; i < J; i++){
			double sum = 0.0;
			for (j = 0; j < J; j++){
				sum	+= fabs(A[i*J + j]);
			}
			anorm = max(anorm,sum);
		}
	}

	// LU-factorize
	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, J, J, A, lda, ipiv);
	if (info != 0){
		printf("Error: dgetrf did not return successfully.\n");
		return;
	}

	// solve
	info = LAPACKE_dgecon(LAPACK_ROW_MAJOR, norm, J, A, lda, anorm, &rcond);
	if (info != 0){
		printf("Error: dgecon did not return successfully.\n");
		return;
	}
}



#endif
