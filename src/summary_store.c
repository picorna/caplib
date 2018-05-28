/*
 *             CAPLIB - CAPsid LIBrary (beta-version)
 *
 *  For Calculations on Icosahedrally Symmetric Virus Capsid Structures
 *
 *  The source codes reported in the following article is contained in this directory,
 *  Shigetaka Yoneda, Yukina Hara-Yamada, Aya Kosugi, Maiko Nanao, Takami Saito, Shunsuke Sato, Nozomu Yamada, and Go Watanabe,
 *  "CAPLIB: A New Program Library for the Modeling and Analysis of Icosahedrally Symmetric Viral Capsids",
 *  ACS Omega 2018, 3, 4458−4465.
 *
 *  The purpose of this library is to analyze directions of rotation axes, calculate cell numbers, 
 *  generate the entire structure of capsid from protomer struture, etc.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, version 2.
 *
 *  Febrary, 2018.
 */
/*
 *	capaxis - calculations on icosahedrally symmetric virus capsid structures
 *            analyse directions of rotation axes, calculate cell numbers, and generate the entire
 *            structure of a capsid from a protomer struture
 *  copyright 2004, 2011, Shigetaka Yoneda
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, version 2.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *	Shigetaka Yoneda, September, 2011.
 */

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "summary_store.h"

void error_store(int, char *);

#define ELIMIT_ROTTYPES 1.0
#define ROTDIFF_LIMIT 10.0

int summary_store(int nfiles,double theta[60],double axis[60][3],char *pdbname, char id[4],
				  int top[3],float top_percent[3],int f_verbose) {
	int naxis_0,naxis_72,naxis_120,naxis_144,naxis_180;
	int i,j;
	double diff,d,d0,d1,d2;
	char message[81];
		
	strncpy(store[nfiles].pdbname,pdbname,81);
	store[nfiles].id[0] = id[0];
	store[nfiles].id[1] = id[1];
	store[nfiles].id[2] = id[2];
	store[nfiles].id[3] = id[3];
	store[nfiles].pdb_type = nfiles;
	for (i=0; i<60; i++) {
		store[nfiles].theta[i] = theta[i];
		for (j=0; j<3; j++) {
			store[nfiles].axis[i][j] = axis[i][j];		}
		
	}

	for (i=0; i<3; i++) {
		store[nfiles].top[i] = top[i];
		store[nfiles].top_percent[i] = top_percent[i];
	}
	
	naxis_0   = 0; // number of unit matrix
	naxis_72  = 0; // number of 5fold rotation matrices (1)
	naxis_120 = 0; // number of 3fold rotation matrices
	naxis_144 = 0; // number of 5fold rotation matrices (2)
	naxis_180 = 0; // number of 2fold rotation matrices
	for (i=0; i<60; i++) {
		if (fabs(theta[i]) < ELIMIT_ROTTYPES) {
			naxis_0++;
		} else if (fabs(theta[i]-72.0) < ELIMIT_ROTTYPES) {
			naxis_72++;
		} else if (fabs(theta[i]-120.0) < ELIMIT_ROTTYPES) {
			naxis_120++;
		} else if (fabs(theta[i]-144.0) < ELIMIT_ROTTYPES) {
			naxis_144++;
		} else if (fabs(theta[i]-180.0) < ELIMIT_ROTTYPES) {
			naxis_180++;
		}
	}
	if (naxis_0 != 1 || naxis_72 != 12 || naxis_120 != 20 || naxis_144 != 12 || naxis_180 != 15) {
		sprintf(message, "No.%d %c%c%c%c errorneous numbers of axes for 0, 72, 120, 144, & 180\n",
				nfiles,id[0],id[1],id[2],id[3]);
		fprintf(stderr, "%s\n",message);
		printf("%s",message);
		error_store(nfiles,message);
		if (f_verbose) {
			printf("axis 0 axis:%1d, 72 axis:%2d, 120 axis:%2d, 144 axis:%2d 180 axis:%2d\n",
				   naxis_0,naxis_72,naxis_120,naxis_144,naxis_180);

		}
		return 1;
	}
	
	/* find the PDB file with the most similar set of rotations */
	for (j=0; j<nfiles; j++) {
		diff = 0.0;
		for (i=0; i<60; i++) {
			d  = theta[i]   - store[j].theta[i];
			d0 = axis[i][0] - store[j].axis[i][0];
			d1 = axis[i][1] - store[j].axis[i][1];
			d2 = axis[i][2] - store[j].axis[i][2];
			diff += d*d + (d0*d0 + d1*d1 + d2*d2)*10000.0;
		}
		if (diff < ROTDIFF_LIMIT) {
//			printf("is equal to the file No. %4d %c%c%c%c\n",j,store[j].id[0],id[1],id[2],id[3]);
			store[nfiles].pdb_type = j;
			return 0;
		} else {
			
		}

	}
	if (f_verbose == 1) {
		printf("No.%d %c%c%c%c is a firstly appeared rotation type.\n",nfiles+1,id[0],id[1],id[2],id[3]);
	}
	return 0;
}