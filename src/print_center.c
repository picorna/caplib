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
#include <stdlib.h>
#include <string.h>

// print the center of atomic coordinates
void print_center(int f_verbose,int nlines,char *pdb_line,double *pdb_xyz,char id[4],double xyzcent[60][3]) {
	
	int j;
	double *xyz;
	char *line;
	double xyzsum[3],xyzdip;
	int natoms;
	double d0,d1,d2;
	
	printf("<< Short Summary >>\nCell numbers of atoms are written at the B-factor position for each ATOM line\n");
	// print the center
	line = pdb_line;
	xyz  = pdb_xyz;
	xyzsum[0] = 0.0;
	xyzsum[1] = 0.0;
	xyzsum[2] = 0.0;
	natoms = 0;
	for (j=0; j<nlines; j++) {
		if (strncmp(line,"ATOM",4) == 0 || strncmp(line,"HETATM",6) == 0){
			xyzsum[0] += *xyz;
			xyzsum[1] += *(xyz+1);
			xyzsum[2] += *(xyz+2);
			xyz += 3;
		}
		line += 81;
		natoms++;
	}
	xyzsum[0] /= natoms;
	xyzsum[1] /= natoms;
	xyzsum[2] /= natoms;
	printf("%c%c%c%c center of atomic coordinates = %10.4f%10.4f%10.4f\n",id[0],id[1],id[2],id[3],
		   xyzsum[0],xyzsum[1],xyzsum[2]);
	
	// print the moment of inertia (mass of all the atoms assumed to be 1)
	line = pdb_line;
	xyz  = pdb_xyz;
	xyzdip = 0.0;
	natoms = 0;
	for (j=0; j<nlines; j++) {
		if (strncmp(line,"ATOM",4) == 0 || strncmp(line,"HETATM",6) == 0){
			d0 = *xyz     - xyzsum[0];
			d1 = *(xyz+1) - xyzsum[1];
			d2 = *(xyz+2) - xyzsum[2];
			xyzdip += d0*d0 + d1*d1 + d2*d2;
			xyz += 3;
		}
		line += 81;
		natoms++;
	}
	xyzdip /= natoms;
	printf("                   inertia moment = %10.4f\n\n\n",xyzdip);
	
}