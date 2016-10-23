/* READ / WRITE
 *
 * REFERENCES
 *  
 * PARAMETERS
 *  J		[input]			number of space segments
 */

#ifndef RW_H
#define RW_H


/* HEADER FILES */
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <algorithm>

using namespace std;

/* PROTOTYPES */

/* IMPLEMENTATIONS */

void write(int J, int M, double dr, double dt, double *params,
           double *R, double *T, 
					 double *H, double *G, double *F,
					 double *P, double *Q, double *Vs){
	int i, j, n, m;
	int J1 = J+1;
	int M1 = M+1;
	const int CHRLEN = 1024;
	double *ptr;
	char *cptr;

	// store parameters in character arrays
	char Ca[CHRLEN], Bo[CHRLEN], Ma[CHRLEN], Dr[CHRLEN], Dt[CHRLEN];
	sprintf(Ca, "_Ca%.0e", params[0]);
	sprintf(Bo, "_Bo%.0e", params[1]);
	sprintf(Ma, "_Ma%.0e", params[2]);
	sprintf(Dr, "_dr%.0e", dr       );
	sprintf(Dt, "_dt%.0e", dt       );

	// replace unwanted characters
	for (n = 0; n < 5; n++){
		if (n == 0) cptr = &Ca[0];
		if (n == 1) cptr = &Bo[0];
		if (n == 2) cptr = &Ma[0];
		if (n == 3) cptr = &Dr[0];
		if (n == 4) cptr = &Dt[0];
		for (i = 0; cptr[i] != '\0'; i++){
			if (cptr[i] == '-') cptr[i] = 'm'; 
			if (cptr[i] == '+') cptr[i] = 'p';
		}
	}

	// concatenate into suffix
	char sfx[CHRLEN];
	sprintf(sfx, "%s%s%s%s%s.txt", Ca, Bo, Ma, Dr, Dt);
	
	// write data to files
	FILE *fp;
	char fn[CHRLEN];
	char hdr[CHRLEN];

	for (n = 0; n < 8; n++){
		if (n == 0){ sprintf(hdr,"../output/r" ); ptr = &R [0]; }
		if (n == 1){ sprintf(hdr,"../output/t" ); ptr = &T [0]; }
		if (n == 2){ sprintf(hdr,"../output/h" ); ptr = &H [0]; }
		if (n == 3){ sprintf(hdr,"../output/g" ); ptr = &G [0]; }
		if (n == 4){ sprintf(hdr,"../output/f" ); ptr = &F [0]; }
		if (n == 5){ sprintf(hdr,"../output/p" ); ptr = &P [0]; }
		if (n == 6){ sprintf(hdr,"../output/q" ); ptr = &Q [0]; }
		if (n == 7){ sprintf(hdr,"../output/vs"); ptr = &Vs[0]; }
		
		//sprintf(fn, "%s%s",hdr,sfx);
		sprintf(fn, "%s.txt",hdr,sfx);
		fp = fopen(fn,"w");
		for (m = 0; m < M1; m++){
			for (i = 0; i < J1; i++){
				fprintf(fp, "%12.10f ", ptr[m*J1 + i]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
		
		printf("Data written to %s.\n",fn);
	}


	//sprintf(fn, "%s%4.4f", "./test", Ca);
	//sprintf(fn, "%s%s", fn, ".txt");

	//printf("%s\n",fn);

//	fp = fopen(fn,"w");
//	fprintf(fp,"testing.\n");
//	fclose(fp);

//	// file name
//	fn = "data.txt";
//
//	ofstream file;
//	file.open(fn.c_str());
//	for (m = 0; m < M1; m++){
//		for (i = 0; i < J1; i++){
//			
//		}
//	}

}



#endif
