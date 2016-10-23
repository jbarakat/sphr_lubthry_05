/* MAIN PROGRAM */

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <math.h>
#include <string>
#include "../include/fd.h"
#include "../include/rw.h"

using namespace std;

int main(){
	int i, j, m, n;
	char fn[1024];
	double params[50] = {0.0};

	// read parameters
	sprintf(fn, "./params.in" );
	readInput(fn, params);

	double Ca    = params[0];
	double Bo    = params[1];
	double Ma    = params[2];
	double tstop = params[3];
	double r1    = params[4];
	double t1    = params[5];
	double dr    = params[6];
	double dt    = params[7];
	double dtrec = params[8];
	
	int J  = r1/dr;
	int N  = t1/dt + 1;
	int M = N*dt/dtrec;
//	N = 1000;
//	M = 100;
	
	int J1 = J+1;
	int M1 = M+1;
	int M1J1 = M1*J1;
	int M1J  = M1*J;

	/* NOTE: The CFL number for this problem is given by
	 *
	 *         dt D
	 *   CFL = ----
	 *         dr^4
	 *
	 * where D is a fourth-order diffusion coefficient. Typically,
	 *
	 *   D ~ O(r1^6 / Ca)    or    D ~ O(r1^4 Ma)
	 *
	 * Having CFL = O(1) is desirable for numerical accuracy, but is not 
	 * required for numerical stability.
	 */
	
	// allocate memory and initialize to zero
	double *R  = (double*) calloc(M1J1, sizeof(double));
	double *T  = (double*) calloc(M1J1, sizeof(double));
	double *H  = (double*) calloc(M1J1, sizeof(double));
	double *G  = (double*) calloc(M1J1, sizeof(double));
	double *F  = (double*) calloc(M1J1, sizeof(double));
	double *P  = (double*) calloc(M1J1, sizeof(double));
	double *Q  = (double*) calloc(M1J1, sizeof(double));
	double *Vs = (double*) calloc(M1J1, sizeof(double));
	
	// time evolution of h, g, and f
	fdevol(J, N, M, dr, dt, params, H, G, F);
	
	// back-calculate p, q, and vs
	fdaux (J, N, M, dr, dt, params, H, G, F, P, Q, Vs);

	// space and time domains
	fdgrid(J, N, M, dr, dt, R, T);

	// write to file
	sprintf(fn, "r" ); write(J, M, dr, dt, params, R , fn);
	sprintf(fn, "t" ); write(J, M, dr, dt, params, T , fn);
	sprintf(fn, "h" ); write(J, M, dr, dt, params, H , fn);
	sprintf(fn, "g" ); write(J, M, dr, dt, params, G , fn);
	sprintf(fn, "f" ); write(J, M, dr, dt, params, F , fn);
	sprintf(fn, "p" ); write(J, M, dr, dt, params, P , fn);
	sprintf(fn, "q" ); write(J, M, dr, dt, params, Q , fn);
	sprintf(fn, "vs"); write(J, M, dr, dt, params, Vs, fn);

	// free memory
	free(R );
	free(T );
	free(H );
	free(G );
	free(F );
	free(P );
	free(Q );
	free(Vs);
	
	return(0);
}







int test(){
	int i, j;
	int J = 10;
	int JJ = J*J;
	double U[JJ], DF[JJ], DB[JJ], DC[JJ], DA[JJ], DD[JJ];
	double h0[J+1], g0[J+1], f0[J+1];
	double h1[J+1], g1[J+1], f1[J+1];
	double params[3];
	double Ca = 0.001;
	double Bo = 1.0;
	double Ma = 0.0;
	params[0] = Ca; 
	params[1] = Bo; 
	params[2] = Ma; 
	double dr = 0.05;
	double dt = dr*dr*Ca/100.0;

	
	// initialize
	for (i = 0; i < J+1; i++){
		f0[i] = 1.0 - 0.5*(i*dr)*(i*dr);
	}

	fdmat(J, U, DF, DB, DC, DA, DD);
//	vector<double> DF(JJ,0), DB(JJ,0), DC(JJ,0), DA(JJ,0), DD(JJ,0);
//	fdmat(J, DF.data(), DB.data(), DC.data(), DA.data(), DD.data());

	fdstep(0, J, dr, dt, params, h0, g0, f0, h1, g1, f1);	

//	for (i = 0; i < J; i++){
//		for (j = 0; j < J; j++){
//			double f = DD[i*J+j];
//			if (f < 0) printf( "%.4f  ",f);
//			else       printf(" %.4f  ",f);
//		}
//		printf("\n");
//	}
	
	// test bi-pentadiagonal algorithm
	double A[JJ], B[JJ], C[JJ], D[JJ];
	double f[J], g[J], u[J], v[J];

	for (i = 0; i < JJ; i++){
		A[i] = 0.0;
		B[i] = 0.0;
		C[i] = 0.0;
		D[i] = 0.0;
	}

	for (i = 0; i < J; i++){
		f[i] = 1.0;
		g[i] = 1.0;
		u[i] = 0.0;
		v[i] = 0.0;
	}
	
	for (i = 0; i < J; i++){
		if (i > 1){
			A[i*J + (i-2)] = 1.0;
			B[i*J + (i-2)] = 2.0;
			C[i*J + (i-2)] = 3.0;
			D[i*J + (i-2)] = 4.0;
		}
		if (i > 0){
			A[i*J + (i-1)] = 5.0;
			B[i*J + (i-1)] = 6.0;
			C[i*J + (i-1)] = 7.0;
			D[i*J + (i-1)] = 8.0;
		}

		A[i*J + i] = 9.1;
		B[i*J + i] = 1.0;
		C[i*J + i] = 2.0;
		D[i*J + i] = 3.0;
		
		if (i < J-1){
			A[i*J + (i+1)] = 4.0;
			B[i*J + (i+1)] = 5.0;
			C[i*J + (i+1)] = 6.0;
			D[i*J + (i+1)] = 7.0;
		}
		if (i < J-2){
			A[i*J + (i+2)] = 8.0;
			B[i*J + (i+2)] = 9.0;
			C[i*J + (i+2)] = 1.0;
			D[i*J + (i+2)] = 2.0;
		}
	}

//	laband(J, DD, 1, 1, DD, 1, 1, A);
//	laband(J, DD, 1, 1, DD, 1, 1, B);
//	laband(J, DD, 1, 1, DD, 1, 1, C);
//	laband(J, DD, 1, 1, DD, 1, 1, D);
//	for (i = 0; i < JJ; i++){
//		A[i] *= 0.1;
//		B[i] *= 0.2;
//		C[i] *= 0.05;
//		D[i] *= 0.3;
//	}

//	for (i = 0; i < J; i++){
//		for (j = 0; j < J; j++){
//			double a = A[i*J+j];
//			if (f < 0) printf( "%.4f  ",a);
//			else       printf(" %.4f  ",a);
//		}
//		printf("\n");
//	}

	double u2[J], v2[J];
	labipent(J,A,B,C,D,f,g,u,v);
	lablock(J,A,B,C,D,f,g,u2,v2);
//	for (i = 0; i < J; i++){
//		printf("%8.4f  %8.4f  \n", u[i], u2[i]);
//	}
//	for (i = 0; i < J; i++){
//		printf("%8.4f  %8.4f  \n", v[i], v2[i]);
//	}
//
//	printf("\n\n");
//
	for (i = 0; i < J; i++){
		for (j = 0; j < J; j++){
			double a = A[i*J+j];
			printf("%8.4f  ",a);
		}
		printf("\n");
	}
//
//	printf("\n\n");

	for (i = 0; i < J; i++){
//		u[i] = 0.0 + i;
		v2[i] = 0.0;
	}

	cblas_dgemv(CblasRowMajor, CblasNoTrans, J, J, 1.0, A, J, u, 1, 1.0, v2, 1);
	laband(J, A, 2, 2, u, v);

	for (i = 0; i < J; i++){
		printf("%8.4f  %8.4f  \n", v[i], v2[i]);
	}

	// compare to LAPACK LU-decomposition

	//fdgrid(J, N, M, dr, dt, R, T);
	//
	//int i, m;
	//for (m = 0; m < M1; m++){
	//	for (i = 0; i < J1; i++){
	//		printf("%4.2f ",T[m*J1 + i]);
	//	}
	//	printf("\n");
	//}


	
	//double h0[J+1], g0[J+1], f0[J+1];
	//double h1[J+1], g1[J+1], f1[J+1];
	//for (i = 0; i < J+1; i++){
	//	h0[i] = 0.0;
	//	g0[i] = 0.0;
	//	f0[i] = 1.0 + 0.5*(i*dr)*(i*dr);
	//	h1[i] = 0.0;
	//	g1[i] = 0.0;
	//	f1[i] = 0.0;
	//}
	//fdstep(J, dr, dt, params, h0, g0, f0, h1, g1, f1);

}
