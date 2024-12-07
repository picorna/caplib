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
#include <string.h>

int cell_number(double *);

/* Although this function was firstly made to genearte the entire structure by repeating 60 rotations 
 of BIOMET lines in a PDB file, this function is used now to print the PDB file with and without rotations, 
 adding some informations on BIOMT, etc, . */

void print_a_pdbs(int f_verbose,int nlines,char *pdb_line,double *pdb_xyz,
				 double rot[60][3][3],double xyzcent[60][3]) {
	
	int j,k;
	double *xyz;
	char *line;
	int nrot,nbiomt;
	int cell;
	
	
	line = pdb_line;
	xyz  = pdb_xyz;
	nbiomt = 0;
	for (j=0; j<nlines; j++) {
		if (strncmp(line,"ATOM",4) == 0 || strncmp(line,"HETATM",6) == 0){
			if(f_verbose == 1) {
				cell = cell_number(xyz);
			}
			for (k=0; k<30; k++) {
			    printf("%c",*(line+k));
			}
			printf("%8.3f%8.3f%8.3f",*xyz,*(xyz+1),*(xyz+2));
			for (k=54; k<60; k++) {
			    printf("%c",*(line+k));
			}
			if(f_verbose == 1) {
				printf("%6.2f",(double) cell);
				printf("%s\n",line+66);
			} else {
				printf("%s\n",line+60);
			}
			xyz += 3;
		} else if ((strncmp(line,"REMARK",6) == 0) && (strncmp(line+10,"   BIOMT",8) == 0)) {
			for (k=0; k<24; k++) {
				printf("%c",*(line+k));
			}
/* printing the axes and angles for the rotation matrices at the BIOMT lines*/
			nrot = nbiomt/3;
			k    = nbiomt%3;
//01234567890123456789012012345678901234567890123456789012345678901234
//REMARK 350   BIOMT2  11  0.809017  0.309017 -0.500000     -208.52413
			if(k==0) {
				printf("%9.6f %9.6f %9.6f       %9.6f\n",
				rot[nrot][k][0],rot[nrot][k][1],rot[nrot][k][2],xyzcent[nrot][k]);
			} else {
				printf("%9.6f %9.6f %9.6f       %9.6f\n",
				rot[nrot][k][0],rot[nrot][k][1],rot[nrot][k][2],xyzcent[nrot][k]);
			}
			nbiomt++;
		} else {
			printf("%s\n",line);
		}
		line += 81;
	}
}
