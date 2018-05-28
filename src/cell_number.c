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

#include "planes.h"
#include "ptable.h"

//
// cell number is determined for the coordinates, xyz[3]
// The 2PLV convention is adopted.
//
// I decided to the pentamer (1/12 of entire structure) number from the cell (1/60 of entire structure) numbers,
//      pentamer number = cell number % 5
// for a atom position very easily. Thus I do not write any functions to calculate the pentamer number.
//
int cell_number(double *xyz) {
	int updown[15],cell,return_cell;
	int i,not,p;
	double up;
	int nequal[60];
	
	//printf("updown ");
	for (i=0;i<15;i++) {
		up = (*xyz)*(rsbc_planes.axis2[i][0]) + *(xyz+1)*(rsbc_planes.axis2[i][1]) + *(xyz+2)*(rsbc_planes.axis2[i][2]);
		if (up < 0.0) {
			updown[i] = -1;
		} else {
			updown[i] = 1;
		}
	}

	// The 12 pentermer cells of icosahedral symmetry
	/*for(cell5=0; cell5<12; cell5++) {
		not = 0;
		nequal5[cell5] = 0;
		for (i=0; i<15; i++) {
			p = rsbc.p12table[cell5][i];
			if (p == 0) {
				//printf("z5 ");
				nequal5[cell5]++;
			} else if (p == updown[i]) {
				nequal5[cell5]++;
			} else {
				//printf("i = %d  ",i);
				not = 1;
				//break;
			}
		}
		printf("%2d",nequal5[cell5]);
		if (not == 1) {
			continue;
		} else {
			printf("\nHell5!\n\n");
			break;
		}
	}
	printf("END5\n");*/
	
	// The standard 60 cells of icosahedral symmetry
	for(cell=0; cell<60; cell++) {
		not = 0;
//		nequal[cell] = 0;
		/*for (i=0; i<15; i++) {
			if ((p = rsbc_ptable.ptable[cell][i]) == 0) {
				continue;
			}
			if (p == updown[i]) {
				continue;
			} else {
				//printf("i = %d  ",i);
				not = 1;
				break;
			}
		}*/
		for (i=0; i<15; i++) {
			p = rsbc.ptable[cell][i];
			if (p == 0) {
				continue;
//				nequal[cell]++;
			} else if (p != updown[i]) {
				not = 1;
				break;
//				nequal[cell]++;
			}
		}
//		printf("%2d",nequal[cell]);
		if (not == 1) {
			nequal[cell] = -1;
			continue;
		} else {
//			printf("\nHell! cell = %d\n",cell);
			return_cell = cell;
			nequal[cell] = 1;
		}
	}
	return return_cell;
}
