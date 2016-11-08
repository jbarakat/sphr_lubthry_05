/* FINITE DIFFERENCE METHOD
 *  Finite difference method for the one-dimensional evolution of film depth and insoluble surfactant
 *  concentration on an initially planar free surface penetrated by a rigid sphere from below.
 *
 *  The space and time domains are discretized on a uniform stencil with grid spacing dr and dt,
 *  respectively. Cylindrical polar coordinates are used in space. Time advancement is carried out
 *  using a trapezoidal method with centered-difference spatial operators (Crank-Nicholson scheme).
 *
 * REFERENCES
 *  von Rosenberg, Methods for the Numerical Solution of Partial Differential Equations (Elsevier, 1969)
 *  
 * PARAMETERS
 *  J		[input]			number of space segments
 *  N		[input]			number of time segments calculated
 *  M		[input]			number of time segments recorded
 *  dr  [input]			spatial grid spacing
 *  dt  [input]			temporal grid spacing
 *
 *  h   [output]		film thickness
 *  g   [output]		excess surface concentration
 *  f   [output]		sphere position
 *  p   [output]		dynamic pressure
 *  q   [output]		volume flux per unit length
 *  vs  [output]		surface velocity
 */

#ifndef FD_H
#define FD_H


/* HEADER FILES */
#include "la.h"

using namespace std;

/* PROTOTYPES */
void fdgrid(int, int, int, double, double, double *, double *);
void fdmat (int, double *, double *, double *, double *, double *, double *);
void fddiff(int, double, double, double *, double *, double *, double *, double *, double *, double *);
void fdstep(int, int, double, double, double *, double *, double *, double *, double *, double *, double *);
void fdevol(int, int, int, double, double, double *, double *, double *, double *);

/* IMPLEMENTATIONS */
// space and time stencils
void fdgrid(int J, int N, int M, double dr, double dt, double *R, double *T){
	int i, n, m;
	int J1 = J+1;
	int N1 = N+1;
	int M1 = M+1;
	int nrec = N/M;
	double r, t;

	m = 0;
	for (n = 0; n < N1; n++){
		t = n*dt;
		if (n % nrec == 0){
			//printf("ts = %d\n",m);
			for (i = 0; i < J1; i++){
				r = i*dr;
				R[m*J1 + i] = r;
				T[m*J1 + i] = t;
			}
			m++;
		}
	}
}

// unit and derivative matrices
void fdmat(int J, double *U, double *DF, double *DB, double *DC, double *DA, double *DD){
	int i, j;
	int JJ = J*J;

	// initialize
	for (i = 0; i < JJ; i++){
		U [i] = 0.0;
		DF[i] = 0.0;
		DB[i] = 0.0;
		DC[i] = 0.0;
		DA[i] = 0.0;
		DD[i] = 0.0;
	}

	// identity matrix
	for (i = 0; i < J; i++){
		U[i*J + i] = 1;
	}

	// derivative matrices
	for (i = 1; i < J; i++){
		// forward
		DF[ i   *J +  i   ] = -1;
		DF[(i-1)*J +  i   ] =  1;

		// backward
		DB[ i   *J +  i   ] =  1;
		DB[ i   *J + (i-1)] = -1;

		// centered
		DC[ i   *J + (i-1)] = -0.5;
		DC[(i-1)*J +  i   ] =  0.5;

		// backward-gradient
		DA[ i   *J +  i   ] =  1 + 0.5/i;
		DA[ i   *J + (i-1)] = -1 + 0.5/i;
	}
	DF[0] = -1; DF[1] =  1;
	DB[0] =  1; DB[1] = -1;
	DC[1] =  0; DC[(J-1)*J + (J-1)] =  1.5; DC[(J-1)*J + (J-2)] = -2; DC[(J-1)*J + (J-3)] = 0.5;
	DA[0] =  4; 
		
	// centered-Laplacian
//	// METHOD #1: input entries manually
//	for (i = 1; i < J-1; i++){
//		DD[ i   *J +  i   ] = -2;
//		DD[ i   *J + (i-1)] =  1 - 0.5/i;
//		DD[ i   *J + (i+1)] =  1 + 0.5/i;
//	}
//	DD[0] = -4; DD[1] = 4; DD[(J-1)*J + (J-1)] = -2; DD[(J-1)*J + (J-2)] = 1 - 0.5/(J-1);
//
//	// METHOD #2: matrix multiplication
//	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, J, J, J, 1.0, DA, J, DF, J, 0.0, DD, J);
	// METHOD #3: banded matrix multiplication
	laband(J, DA, 0, 1, DF, 1, 0, DD);	
}

// nonlinear diffusion coefficients
//  k = diffusion of h-stuff
//  m = diffusion of g-stuff
void fddiff(int J, double k, double m,
            double *h, double *g, double *f,
            double *kh, double *kg, double *mh, double *mg){
	int i;

	for (i = 0; i < J; i++){
		kh[i] = k*pow((h[i] + f[i]),3.0)/3.0;
		kg[i] = k*pow((h[i] + f[i]),2.0)/2.0;
		mh[i] = m*pow((h[i] + f[i]),2.0)/2.0;
		mg[i] = m*    (h[i] + f[i])         ;
	}

	for (i = 0; i < J; i++){
		kg[i] *= (g[i] + 1.0);
		mg[i] *= (g[i] + 1.0);
	}
}

// advance +1 timestep
//  h0, g0, h1, g1 are (J+1)-vectors
void fdstep(int n, int J, double dr, double dt, double *params,
            double *h0, double *g0, double *f0,
					  double *h1, double *g1, double *f1){
	int i, j;
	int JJ = J*J;
	double dr2 = dr*dr;
	double dth = 0.5*dt;
	double t   = n*dt - dth;

	// allocate memory and initialize to zero
	double *hm = (double*) calloc(J , sizeof(double));
	double *gm = (double*) calloc(J , sizeof(double));
	double *fm = (double*) calloc(J , sizeof(double));
	double *h  = (double*) calloc(J , sizeof(double));
	double *g  = (double*) calloc(J , sizeof(double));
	double *sh = (double*) calloc(J , sizeof(double));
	double *sg = (double*) calloc(J , sizeof(double));
	double *s1 = (double*) calloc(J , sizeof(double));
	double *s2 = (double*) calloc(J , sizeof(double));
	double *kh = (double*) calloc(J , sizeof(double));
	double *mh = (double*) calloc(J , sizeof(double));
	double *kg = (double*) calloc(J , sizeof(double));
	double *mg = (double*) calloc(J , sizeof(double));
	double *Kh = (double*) calloc(JJ, sizeof(double));
	double *Mh = (double*) calloc(JJ, sizeof(double));
	double *Kg = (double*) calloc(JJ, sizeof(double));
	double *Mg = (double*) calloc(JJ, sizeof(double));
	double *U  = (double*) calloc(JJ, sizeof(double));
	double *DF = (double*) calloc(JJ, sizeof(double));
	double *DB = (double*) calloc(JJ, sizeof(double));
	double *DC = (double*) calloc(JJ, sizeof(double));
	double *DA = (double*) calloc(JJ, sizeof(double));
	double *DD = (double*) calloc(JJ, sizeof(double));
	double *P  = (double*) calloc(JJ, sizeof(double));
	double *DP = (double*) calloc(JJ, sizeof(double));
	double *A  = (double*) calloc(JJ, sizeof(double));
	double *B  = (double*) calloc(JJ, sizeof(double));
	double *C  = (double*) calloc(JJ, sizeof(double));
	double *D  = (double*) calloc(JJ, sizeof(double));

	// physical parameters
	double Ca    = params[0]; double k0 = 1.0/Ca;
	double Bo    = params[1];
	double Ma    = params[2]; double m0 = Ma;
	double tstop = params[3];
	
	// initialize source vectors
	if (t < tstop){
		for (i = 0; i < J; i++){
			sh[i] = 1.0;
			sg[i] = 0.0;
		}
	}
	else { // add evaporation
		double H = 0.001; // fitting parameter
		for (i = 0; i < J; i++){
			sh[i] += -H/pow(t,0.5);
		}
	}

	for (i = 0; i < J; i++){
		sh[i] *= dt;
	}

	// derivative matrices
	fdmat(J, U, DF, DB, DC, DA, DD);
	for (i = 0; i < JJ; i++){
		DF[i] /= dr;
		DB[i] /= dr;
		DC[i] /= dr;
		DA[i] /= dr;
		DD[i] /= dr2;
	}
	
	// pressure and pressure gradient operator matrices
	for (i = 0; i < JJ; i++) P[i] = Bo*U[i] - DD[i];
	laband(J, DF, 1, 0, P, 1, 1, DP);	

	// evaluate functions h, g, f at (j+1/2,n)
	for (i = 0; i < J; i++){
		hm[i] = 0.5*(h0[i] + h0[i+1]);
		gm[i] = 0.5*(g0[i] + g0[i+1]);
		fm[i] = 0.5*(f0[i] + f0[i+1]);
	}

	// diffusion coefficients kh, kg, mh, mg at (j,n)
	// (NO PROJECTION TO THE 1/2 TIME LEVEL - NOT TECHNICALLY CORRECT...)
	fddiff(J, k0, m0, hm, gm, fm, kh, kg, mh, mg);
	for (i = 0; i < J; i++){
		Kh[i*J + i] = -dth*kh[i];
		Kg[i*J + i] = -dth*kg[i];
		Mh[i*J + i] = -dth*mh[i];
		Mg[i*J + i] = -dth*mg[i];
	}
	
	// operator matrices A, B, C, D
	laband(J, DA, 0, 1, Kh, 0, 0, A); laband(J, A , 0, 1, DP, 2, 1, A);
	laband(J, DA, 0, 1, Mh, 0, 0, B); laband(J, B , 0, 1, DF, 1, 0, B);
	laband(J, DA, 0, 1, Kg, 0, 0, C); laband(J, C , 0, 1, DP, 2, 1, C);
	laband(J, DA, 0, 1, Mg, 0, 0, D); laband(J, D , 0, 1, DF, 1, 0, D);

	// update source functions at (j,n)
	for (i = 0; i < J; i++) { s1[i] = 0.0; s2[i] = 0.0; }
	laband(J, A, 2, 2, h0, s1);
	laband(J, B, 2, 2, g0, s2);
	for (i = 0; i < J; i++){
		sh[i] += h0[i] - s1[i] - s2[i];
	}
	
	for (i = 0; i < J; i++) { s1[i] = 0.0; s2[i] = 0.0; }
	laband(J, C, 2, 2, h0, s1);
	laband(J, D, 2, 2, g0, s2);
	for (i = 0; i < J; i++){
		sg[i] += g0[i] - s1[i] - s2[i];
	}
	
	// add the identity matrix to A and D
	for (i = 0; i < J; i++){
		A[i*J + i] += 1.0;
		D[i*J + i] += 1.0;
	}
	
//	/*--- for debugging ---*/
//	for (i = 0; i < J; i++){
//		for (j = 0; j < J; j++){
//			printf("%9.4f ", A[i*J + j]);
//		}
//		printf("\n");
//	}
//
//	double rcond;
//	double A2[JJ];
//	for (i = 0; i < JJ; i++) A2[i] = A[i];
//	lacond(J, A2, rcond);
//	printf("\n%0.2e\n",1.0/rcond);
//
//	/*---*/
	
	// solve linear system
	//lablock (J, A, B, C, D, sh, sg, h, g);
	labipent(J, A, B, C, D, sh, sg, h, g);
	
	// copy solution and assign right BCs
	for (i = 0; i < J; i++){
		h1[i] = h[i];
		g1[i] = g[i];
	}
	h1[J] = 0.0;
	g1[J] = 0.0;

	// advance f
	for (i = 0; i < J+1; i++){
		if (t < tstop)  f1[i] = f0[i] - dt;
		else            f1[i] = f0[i]     ;
	}

	// free memory
	free(hm);
	free(gm);
	free(fm);
	free(h );
	free(g );
	free(sh);
	free(sg);
	free(s1);
	free(s2);
	free(kh);
	free(mh);
	free(kg);
	free(mg);
	free(Kh);
	free(Mh);
	free(Kg);
	free(Mg);
	free(U );
	free(DF);
	free(DB);
	free(DC);
	free(DA);
	free(DD);
	free(P );
	free(DP);
	free(A );
	free(B );
	free(C );
	free(D );
}

// time evolution
void fdevol(int J, int N, int M, double dr, double dt, double *params,
					  double *H, double *G, double *F){
	int i, j, n, m, nrec;
	int J1 = J+1;
	int N1 = N+1;
	int M1 = M+1;
	double t, r;
	
	// allocate memory and initialize to zero
	double *h0 = (double*) calloc(J1, sizeof(double));
	double *g0 = (double*) calloc(J1, sizeof(double));
	double *f0 = (double*) calloc(J1, sizeof(double));
	double *h1 = (double*) calloc(J1, sizeof(double));
	double *g1 = (double*) calloc(J1, sizeof(double));
	double *f1 = (double*) calloc(J1, sizeof(double));

	// initialize record counter and wait time
	m = 0;
	nrec = N/M;

	// initial conditions
	for (i = 0; i < J1; i++){
		r = i*dr;
		H[0*J1 + i] = 0.0;
		G[0*J1 + i] = 0.0;
		F[0*J1 + i] = 1.0 + 0.5*r*r;
		
		h0[i] = H[0*J1 + i];
		g0[i] = G[0*J1 + i];
		f0[i] = F[0*J1 + i];
	}
	
	// time-evolve
	m = 0;
	for (n = 1; n < N1; n++){
		printf("n = %d / %d\n", n, N);
		fdstep(n, J, dr, dt, params, h0, g0, f0, h1, g1, f1);
		
		for (i = 0; i < J1; i++){
			h0[i] = h1[i];
			g0[i] = g1[i];
			f0[i] = f1[i];
		}
		
		if (n % nrec == 0){
			printf("Time step recorded.\n");
			m++;

			for (i = 0; i < J1; i++){
				H[m*J1 + i] = h1[i];
				G[m*J1 + i] = g1[i];
				F[m*J1 + i] = f1[i];
			}
		}
	}
	
	// free memory
	free(h0);
	free(g0);
	free(f0);
	free(h1);
	free(g1);
	free(f1);
}

// back-calculate p, q, and vs
void fdaux(int J, int N, int M, double dr, double dt, double *params,
           double *H, double *G, double *F, double *P, double *Q, double *Vs){
	int i, j, n, m, nrec;
	int ii, jj;
	int JJ = J*J;
	int J1 = J+1;
	int N1 = N+1;
	int M1 = M+1;
	int J1J1 = J1*J1;
	double dr2 = dr*dr;
	double dth = 0.5*dt;

	// allocate memory and initialize to zero
	double *h  = (double*) calloc(J , sizeof(double));
	double *g  = (double*) calloc(J , sizeof(double));
	double *f  = (double*) calloc(J , sizeof(double));
	double *p  = (double*) calloc(J , sizeof(double));
	double *q  = (double*) calloc(J , sizeof(double));
	double *vs = (double*) calloc(J , sizeof(double));
	double *s1 = (double*) calloc(J , sizeof(double));
	double *s2 = (double*) calloc(J , sizeof(double));
	double *kh = (double*) calloc(J , sizeof(double));
	double *mh = (double*) calloc(J , sizeof(double));
	double *kg = (double*) calloc(J , sizeof(double));
	double *mg = (double*) calloc(J , sizeof(double));
	double *Kh = (double*) calloc(JJ, sizeof(double));
	double *Mh = (double*) calloc(JJ, sizeof(double));
	double *Kg = (double*) calloc(JJ, sizeof(double));
	double *Mg = (double*) calloc(JJ, sizeof(double));
	double *U  = (double*) calloc(JJ, sizeof(double));
	double *DF = (double*) calloc(JJ, sizeof(double));
	double *DB = (double*) calloc(JJ, sizeof(double));
	double *DC = (double*) calloc(JJ, sizeof(double));
	double *DA = (double*) calloc(JJ, sizeof(double));
	double *DD = (double*) calloc(JJ, sizeof(double));
	double *A  = (double*) calloc(JJ, sizeof(double));
	double *B  = (double*) calloc(JJ, sizeof(double));
	double *C  = (double*) calloc(JJ, sizeof(double));
	double *D  = (double*) calloc(JJ, sizeof(double));
	double *E  = (double*) calloc(JJ, sizeof(double));

	// initialize record counter and wait time
	m = 0;
	nrec = N/M;

	// physical parameters
	double Ca    = params[0]; double k0 = 1.0/Ca;
	double Bo    = params[1];
	double Ma    = params[2]; double m0 = Ma;
	double tstop = params[3];
	
	// derivative matrices
	fdmat(J, U, DF, DB, DC, DA, DD);
	for (i = 0; i < JJ; i++){
		DF[i] /= dr;
		DB[i] /= dr;
		DC[i] /= dr;
		DA[i] /= dr;
		DD[i] /= dr2;
	}

	// pressure operator
	for (i = 0; i < JJ; i++)
		E[i] = Bo*U[i] - DD[i];

	// advance in time
	m = 0;
	for (n = 0; n < N1; n++){
		if (n % nrec == 0){
			// call solution vectors
			for (i = 0; i < J; i++){
				h[i] = H[m*J1 + i];
				g[i] = G[m*J1 + i];
				f[i] = F[m*J1 + i];
			}
			
			// diffusion coefficients kh, kg, mh, mg at (j+1/2,n)
			fddiff(J, k0, m0, h, g, f, kh, kg, mh, mg);
			for (i = 0; i < J; i++){
				Kh[i*J + i] = -kh[i];
				Kg[i*J + i] = -kg[i];
				Mh[i*J + i] = -mh[i];
				Mg[i*J + i] = -mg[i];
			}
	
			// operator matrices
			laband(J, Kh, 0, 0, DC, 1, 1, A);
			laband(J, Mh, 0, 0, DC, 1, 1, B);
			laband(J, Kg, 0, 0, DC, 1, 1, C);
			laband(J, Mg, 0, 0, DC, 1, 1, D);
			
			// dynamic pressure
			laband(J, E, 1, 1, h, p);
			
			// volume flux per unit length
			laband(J, A, 1, 1, p, s1);
			laband(J, B, 1, 1, g, s2);
			for (i = 0; i < J; i++) q[i] = s1[i] + s2[i];
			
			// meridional surface velocity
			laband(J, C, 1, 1, p, s1);
			laband(J, D, 1, 1, g, s2);
			for (i = 0; i < J; i++) vs[i] = s1[i] + s2[i];
			
			// store results
			for (i = 0; i < J; i++){
				P [m*J1 + i] = p [i];
				Q [m*J1 + i] = q [i];
				Vs[m*J1 + i] = vs[i];
			}
			P [m*J1 + J] = 0.0;
			Q [m*J1 + J] = 0.0;
			Vs[m*J1 + J] = 0.0;
			
			m++;
		}
	}
	
	// free memory
	free(h );
	free(g );
	free(f );
	free(p );
	free(q );
	free(vs);
	free(s1);
	free(s2);
	free(kh);
	free(mh);
	free(kg);
	free(mg);
	free(Kh);
	free(Mh);
	free(Kg);
	free(Mg);
	free(U );
	free(DF);
	free(DB);
	free(DC);
	free(DA);
	free(DD);
	free(A );
	free(B );
	free(C );
	free(D );
	free(E );
}

#endif
