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

// write to file
void write(int J, int M, double dr, double dt, double *params,
           double *F, char *f){
	int i, j, n, m;
	int J1 = J+1;
	int M1 = M+1;
	const int CHRLEN = 1024;
	char *cptr;
	double *ptr = &F[0];
	double t1, tstop;
	char Ca[CHRLEN], Bo[CHRLEN], Ma[CHRLEN], Dr[CHRLEN], Dt[CHRLEN], Tstop[CHRLEN];
	char pfx[CHRLEN], sfx[CHRLEN], fn[CHRLEN];
	FILE *fp;

	// store parameters in character arrays
	sprintf(Ca   , "Ca%.0e"   , params[0]);
	sprintf(Bo   , "Bo%.0e"   , params[1]);
	sprintf(Ma   , "Ma%.0e"   , params[2]);
	sprintf(Dr   , "dr%.0e"   , dr       );
	sprintf(Dt   , "dt%.0e"   , dt       );
	sprintf(Tstop, "tstop%.1f", params[3]);

	// get tstop and t1
	tstop = params[3];
	t1    = params[5];

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

	cptr = &Tstop[0];
	for (i = 0; cptr[i] != '\0'; i++){
		if (cptr[i] == '.') cptr[i] = '-'; 
	}

	// create filename
	if (tstop > t1)
		sprintf(pfx, "../output/go/%s/%s/%s/%s/%s/", Ca, Bo, Ma, Dr, Dt);
	else
		sprintf(pfx, "../output/stop/%s/%s/%s/%s/%s/%s/", Ca, Bo, Ma, Tstop, Dr, Dt);
	sprintf(sfx, "%s.txt", f);
	sprintf(fn , "%s%s", pfx, sfx);

	// write data to files
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


// read input file
void readInput(char *fn, double *params){
	FILE *fp;
	const int CHRLEN = 1024;
	char pch[CHRLEN];
	double pdb;
	int i, ch;

	fp = fopen(fn, "r");
	if (fp == NULL)
		perror("Error opening file.");
	else {
		i  = 0;
		ch = fscanf(fp, "%s", pch);
		while (ch != EOF){
			fscanf(fp, "%s", pch);
			fscanf(fp, "%s", pch);
			pdb = atof(pch);
			params[i] = pdb;
			
			i++;
			ch = fscanf(fp, "%s", pch);
		}
		printf("Input parameters read from %s.\n", fn);
	}
	fclose(fp);
}



#endif
