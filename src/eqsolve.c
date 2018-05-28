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

// determine a rotation angle (return value) and a rotation axis (axis) for a rotation matrix (rot0)
double eqsolv(double rot0[3][3],double axis[3]) {
	double rot[3][3],work[3][3];
	double axinor[3],axino2[3],outer[3];
	int maxcol,maxrow,nexcol,nexrow,ione;
	double valmax,valnex,div,tmp,div1,div2,costhe;
	double small;
	int ismall;
	double return_val;
	int i,j;
	double deg = 180.0/3.14159265358979323;
	
	/* (0) adding -1 at the diagonal elements */	
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			rot[i][j] = rot0[i][j];
		}
		rot[i][i] = rot0[i][i] - 1.0;
	}

	/* (1) finding the larges element */
	valmax = 0.0;
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			if(fabs(rot[i][j]) > valmax) {
				maxcol = i;
				maxrow = j;
				valmax = fabs(rot[i][j]);
			}
		}
	}

	/* (2) sweeping by the largest column (maxcol) */
	for (i=0; i<3; i++) {
		div = rot[i][maxrow] / rot[maxcol][maxrow];
		for (j=0; j<3; j++) {
			work[i][j] = rot[i][j] - div*rot[maxcol][j];
		}
		if(i != maxcol) {work[i][maxrow] = 0.0;}
	}
	div = 1.0 / rot[maxcol][maxrow];
	for (j=0; j<3; j++) {
		work[maxcol][j] = div*rot[maxcol][j];
	}
	work[maxcol][maxrow] = 1.0;
	
	/* (3) finding the largest except the element at the maxcol column */
	nexcol = 0;
	nexrow = 0;
	valnex = 0.0;
	for (i=0; i<3; i++) {
		if(i != maxcol) {
			for (j=0; j<3; j++) {
				if (fabs(work[i][j]) > valnex) {
					nexcol = i;
					nexrow = j;
					valnex = fabs(work[i][j]);
				}
			}
		}
	}

	/* (4) sweeping by the 2nd largest column (nexcol) */
	for (i=0; i<3; i++) {
		div = work[i][nexrow] / work[nexcol][nexrow];
		for (j=0; j<3; j++) {
			rot[i][j] = work[i][j] - div*work[nexcol][j];
		}
		if(i != nexcol) {rot[i][nexrow] = 0.0;}
	}
	div = 1.0 / work[nexcol][nexrow];
	for (j=0; j<3; j++) {
		rot[nexcol][j] = div*work[nexcol][j];
	}
	rot[nexcol][nexrow] = 1.0;	
		
	/* (5) determining the solution of the eigenvector (normalized) */
	for (i=0; i<3; i++) {
		if((i != maxcol) && (i != nexcol)) {
			ione = i;
			break;
		}
	}
	axis[maxrow] = -rot[maxcol][ione];
	axis[nexrow] = -rot[nexcol][ione];
	axis[ione]   = 1.0;
	div = 1.0 / sqrt(axis[0]*axis[0]+axis[1]*axis[1]+axis[2]*axis[2]);
	for (i=0; i<3; i++) {
		axis[i] = div*axis[i];
	}

	/* (6) determining the rotation angle */
	/* (6.1) find the smallest element of axis */
	small = 1.0e30;
	for (i=0; i<3; i++) {
		if(axis[i] < small) {
			ismall = i;
			small = axis[i];
		}
	}
	for (i=0; i<3; i++) {
		axinor[i] = 0.0;
	}
	axinor[ismall] = 1.0;
	/* the Gram-Shmidt orthogonalization */
	tmp = axis[0]*axinor[0]+axis[1]*axinor[1]+axis[2]*axinor[2];
	for (i=0; i<3; i++) {
		axinor[i] = axinor[i] - tmp*axis[i];
	}
	/* This axinor is a normal vector to axis */
	for (i=0; i<3; i++) {
		axino2[i] = rot0[i][0]*axinor[0] + rot0[i][1]*axinor[1] + rot0[i][2]*axinor[2];
	}
	div1 = axinor[0]*axinor[0]+axinor[1]*axinor[1]+axinor[2]*axinor[2];
	div2 = axino2[0]*axino2[0]+axino2[1]*axino2[1]+axino2[2]*axino2[2];
	costhe = (axinor[0]*axino2[0]+axinor[1]*axino2[1]+axinor[2]*axino2[2]) / sqrt(div1*div2);
	if(costhe <= -1.0) {
		return_val = 180.0;
	} else if (costhe >= +1.0) {
		return_val= 0.0;
	} else {
		return_val = deg*acos(costhe);
		outer[0] = axinor[1]*axino2[2]-axinor[2]*axino2[1];
		outer[1] = axinor[2]*axino2[0]-axinor[0]*axino2[2];
		outer[2] = axinor[0]*axino2[1]-axinor[1]*axino2[0];
		tmp = axis[0]*outer[0]+axis[1]*outer[1]+axis[2]*outer[2];
		if(tmp < 0.0) {
			axis[0] = -axis[0];
			axis[1] = -axis[1];
			axis[2] = -axis[2];
		}
	}
	return return_val;
}