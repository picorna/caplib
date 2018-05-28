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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int cell_number(double *);

// calculate cell number for all the atoms and analyse which cell contains the most many atoms
void calc_cell(int f_verbose,int natoms,double *pdb_xyz,double rot[60][3][3],
			   int top[3],float top_percent[3]) {
	
	int i,j,k;
	double *xyz;
	int cell_plot[60],tmp[3];
	
	if (natoms <= 0) {
		top[0] = -1;
		top[1] = -1;
		top[2] = -1;
		top_percent[0] = -1.;
		top_percent[1] = -1.;
		top_percent[2] = -1.;
		return;
	}

	// determine the cell number of the PDB coordinates
	for (i=0; i<60; i++) {
		cell_plot[i] = 0;
	}
	xyz = pdb_xyz;
	for (i=0; i<natoms; i++) {
		//printf("cell number %d",cell_number(xyz));
		cell_plot[cell_number(xyz)]++;
		/*if (cell_number(xyz) == 58) {
			printf("%d %8.3f %8.3f %8.3f\n",i,*xyz,*(xyz+1),*(xyz+2));
		}*/
		xyz += 3;
	}
	tmp[0] = 0;
	tmp[1] = 0;
	tmp[2] = 0;
	for (i=0; i<60; i++) {
		for (j=0; j<3; j++) {
			if (tmp[j] < cell_plot[i]) {
				for (k=1; j<=k; k--) {
					tmp[k+1] = tmp[k];
					top[k+1] = top[k];
				}
				tmp[j] = cell_plot[i];
				top[j] = i;
				break;
			}
		}
	}
	top_percent[0] = (float) tmp[0];
	top_percent[1] = (float) tmp[1];
	top_percent[2] = (float) tmp[2];
	top_percent[0] /= natoms;
	top_percent[1] /= natoms;
	top_percent[2] /= natoms;
}
