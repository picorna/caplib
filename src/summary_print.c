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

#include "summary_store.h"

void print_rot(int,int,int);

void summary_print(int nfiles) {
	int i,j,k,m10,m10end,mt;
	int mt10[]={1,12,20,12,15};
	int nerrors,npdbtypes;
	
	nerrors = 0;
	npdbtypes = 0;
	for (i=0; i<nfiles; i++) {
		if (store[i].pdb_type == i) {
			npdbtypes++;			
		}
		if (&(store[i].error) == NULL) {
			nerrors++;
		}
	}
	
	printf("Total %4d files, %4d types of rotations, %4d too severe error files\n\n", nfiles, npdbtypes, nerrors);
	printf("< Classification >\n");
	printf("(Rotation type and PDB ID (sequential number of the file))\n");
	
	npdbtypes = 0;
	for (i=0; i<nfiles; i++) {
		if (store[i].pdb_type == i) {
			npdbtypes++;
			printf("\n< Type %3d >\n",npdbtypes);
			m10 = 0;
			for (mt=0; mt<5; mt++) {
				m10end = m10+mt10[mt];
				// print rotation angles and axes
				print_rot(i,m10,m10end);
				m10 = m10end;
			}
			for (j=i; j<nfiles; j++) {
				if (store[j].pdb_type == i) {
					printf("%c%c%c%c(%3d) ",store[j].id[0],store[j].id[1],store[j].id[2],store[j].id[3],j+1);
					printf("1st, 2nd, 3rd major cells = %2d(%2.0f%%)",
						   store[j].top[0],store[j].top_percent[0]*100.);
					for (k=1; k<3; k++) {
						if (store[j].top_percent[k] > 0.) {
							printf(" %2d(%2.0f%%)",
								   store[j].top[k],store[j].top_percent[k]*100.0);
						}
					}
					printf("\n");
				}
			}
			printf("\n");
		}

	}
		
	return;
}
