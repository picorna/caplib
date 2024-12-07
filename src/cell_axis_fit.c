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
#include <math.h>
#include "planes.h"

extern void make_orthogonal_matrix(double [3][3],double [3],double [3]);

/* 1) find the center of all the atoms
 2) find the nearest 5fold, 3fold, two 2fold axes to the center determined above
 3) find a rotation matrix to overlay the 5fold, 3fold and one 2fold axes determined above onto 
 the equivalent axes of the 2PLV protomer
 4) rotate the atomic coordinates by the rotation matrix determined above
 5) calculate the cell numbers of the atoms for the rotated atom (this 5) is to be calculated in calc_cell.c
 */
int cell_axis_fit(int f_verbose,int natoms,double *pdb_xyz,double axis[60][3],
				   double theta[60],double rot_fit[3][3]){
// axis[60][3] is rotation axes determined by analyse_biomt.c from the BIOMT lines
// in the PDB file
	
	double *xyz;
	double center[3];
	double d,dmax[5];
	int dmaxaxis[5];
	int i,j,k;
	int iaxis5,iaxis3;
	double amatrix[3][3],bmatrix[3][3];
	int debug = 0;
//	
// calculate the center of protein, center[3]
	xyz = pdb_xyz;
	center[0] = 0.;
	center[1] = 0.;
	center[2] = 0.;
	for (i=0; i<natoms; i++) {
		center[0] += *(xyz++);
		center[1] += *(xyz++);
		center[2] += *(xyz++);
	}
	center[0] /= natoms;
	center[1] /= natoms;
	center[2] /= natoms;
//	printf("center = %8.3f %8.3f %8.3f\n",center[0],center[1],center[2]);
// calculate the nearest 5 rotation axes to the center, center[3]
	dmax[0] = 0.;
	dmax[1] = 0.;
	dmax[2] = 0.;
	dmax[3] = 0.;
	dmax[4] = 0.;
	for (i=0; i<60; i++) {
		if(fabs(theta[i]) < 0.001) continue;
		d = center[0]*axis[i][0] + center[1]*axis[i][1] + center[2]*axis[i][2];
		d /= sqrt(axis[i][0]*axis[i][0] + axis[i][1]*axis[i][1] + axis[i][2]*axis[i][2]);
		if (debug == 1) {
			printf("i=%d  theta=%f  d=%f\n",i,theta[i],d);
		}
		for (j=0; j<5; j++) {
			if (d > dmax[j]) {
				for (k=3; k>=j; k--) {
					dmax[k+1] = dmax[k];
					dmaxaxis[k+1] = dmaxaxis[k];
				}
				dmax[j] = d;
				dmaxaxis[j] = i;
				break;
			}
		}
	}
//
	iaxis5 = -1;
	iaxis3 = -1;
	// find the nearest 5fold axes
	for (j=0; j<5; j++) {
		i = dmaxaxis[j];
		if (debug == 1) {
			printf("5fold loop: i = %d   theta = %f\n",i,theta[i]);
		}
		if (fabs(theta[i]-72.) < 1.0 || fabs(theta[i]-144.) < 1.0) {
			iaxis5 = i;
			break;
		}
	}
	if(iaxis5 == -1) {
		fprintf(stderr,"no nearest 5fold axis in axis fitting in cell_axis.c\n");
		return EXIT_FAILURE;
	}	
	// find the nearest 3fold axes
	for (j=0; j<5; j++) {
		i = dmaxaxis[j];
		if (debug == 1) {
			printf("3fold loop: i = %d   theta = %f\n",i,theta[i]);
		}
		if (fabs(theta[i]-120.) < 1.0 || fabs(theta[i]-240.) < 1.0) {
			iaxis3 = i;
			break;
		}
	}
	if(iaxis3 == -1) {
		fprintf(stderr,"no nearest 3fold axis fitting in cell_axis.c\n");
		return EXIT_FAILURE;
	}	
	if (debug == 1) {
		printf("");
	}

/*
 By the above, the nearest 5fold and 3fold axes (iaxis5 and iaxis3) were determined.
 The cell (5fold and 3fold) from the BIOMT lines of PDB file is overalaid onto the
 2PLV one by a rotation matrix, rot_fit[3][3]. rot_fit is determined by the following
 equation,
		A = rot_fit * B
 where A and B are orthogonal matrices, and * denotes multiplication of matrices.
 A is built from the two axes (iaxis5 and iaxis3) of 2PLV cell zero and B is built from
 the three axes of BIOMT lines similarly. Thus,
		rot_fit = A * tB
where tB is transpose of B, because the inverse matrix of B is the transpose of B.
*/
	/*double xtmp1[3] = {1.,0.,0.};
	double xtmp2[3] = {0.,1.,0.};
	double ytmp1[3] = {0.,1.,0.};
	double ytmp2[3] = {0.,0.,1.};
	double ztmp[3] = {2.,1.,0.};
	double ztmp2[3];
	make_orthogonal_matrix(amatrix,xtmp1,xtmp2);
	make_orthogonal_matrix(bmatrix,ytmp1,ytmp2);

	for (i=0; i<3; i++) {
	for (j=0; j<3; j++) {
	rot_fit[i][j] = amatrix[i][0]*bmatrix[j][0] + amatrix[i][1]*bmatrix[j][1] + amatrix[i][2]*bmatrix[j][2];
	}
	}
	for (i=0; i<3; i++) {
	ztmp2[i] = rot_fit[i][0]*ztmp[0] + rot_fit[i][1]*ztmp[1] + rot_fit[i][2]*ztmp[2];
	}
	printf("ztmp2 = %8.3f %8.3f %8.3f\n",ztmp2[0],ztmp2[1],ztmp2[2]);
	exit(1);
	*/
//	printf("iaxis5 = %d   iaxis3 = %d\n",iaxis5,iaxis3);

	make_orthogonal_matrix(bmatrix,&axis[iaxis5][0],&axis[iaxis3][0]);
//
// cell 0 has 5fold axis, 0, and 3 fold axis, 1-2-3 (cdh, 0-1-2)
//
	make_orthogonal_matrix(amatrix,&rsbc_planes.axis5[0][0],&rsbc_planes.axis3[4][0]);
//
	for (i=0; i<3; i++) {
	for (j=0; j<3; j++) {
	rot_fit[i][j] = amatrix[i][0]*bmatrix[j][0] + amatrix[i][1]*bmatrix[j][1] + amatrix[i][2]*bmatrix[j][2];
	}
//	if (debug == 2) {
//	printf("rot_fit = %8.3f%8.3f%8.3f\n",rot_fit[i][0],rot_fit[i][1],rot_fit[i][2]);
//	}
	}
	return EXIT_SUCCESS;
}
