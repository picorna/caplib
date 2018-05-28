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
#include <math.h>
#ifndef UNITMATRIX_LIMIT
#define UNITMATRIX_LIMIT 0.001
#endif
#ifndef ANGLE_PRECISION
#define ANGLE_PRECISION 1.0
#endif

void sort_double(int, double [], int i[]);
double eqsolv(double [3][3],double [3]);

int analyse_biomt(int f_verbose,double rot[60][3][3],double theta[60],double axis[60][3],int do_sort) {
	int i,isort[60],should_return = 0;
	double determinant;
	double val[60];
	double work_theta[60],work_axis[60][3];
	// if do_sort == 1, then sorting is applied for axes and theta.
	
	/* loop over the 60 matrices */
	should_return = 0;
	for(i=0;i<60;i++) {
		/* (1.9) examine whether the determinant is one */
		determinant = rot[i][0][0]*rot[i][1][1]*rot[i][2][2] 
			+ rot[i][0][1]*rot[i][1][2]*rot[i][2][0] 
			+ rot[i][0][2]*rot[i][1][0]*rot[i][2][1]
			- rot[i][0][2]*rot[i][1][1]*rot[i][2][0] 
			- rot[i][0][0]*rot[i][1][2]*rot[i][2][1] 
			- rot[i][0][1]*rot[i][1][0]*rot[i][2][2];
		if(fabs(determinant-1.0) > 0.01) {
			fprintf(stderr,"Error: determinant of a BIOMT matrix is %6.2f (not 1.0)\n",determinant);
			should_return = 1;
		}
	}
	if (should_return==1) {
		return 1;
	}
	
	//	int k;
	for(i=0;i<60;i++) {
		/* (2) calculating the rotation angle and axis for rot (eigenvector and eigenvalue) */
		//		for(k=0;k<3;k++) {
		//			printf("%f %f %f\n",rot[i][k][0],rot[i][k][1],rot[i][k][2]);
		//		}
		if ((fabs(rot[i][0][0]-1.0) < UNITMATRIX_LIMIT) && 
			(fabs(rot[i][1][1]-1.0) < UNITMATRIX_LIMIT) && 
			(fabs(rot[i][2][2]-1.0) < UNITMATRIX_LIMIT)) {
			axis[i][0] = 0.0;
			axis[i][1] = 0.0;
			axis[i][2] = 0.0;
			theta[i]   = 0.0;
		} else {
			theta[i] = fabs(eqsolv((double (*)[3])&(rot[i][0][0]),(double (*))&(axis[i][0])));
		}
	}
	
	// Diction of eigenvector is made inverse to make the X coordinate to be positive when the rotation angle is 180.
	for(i=0;i<60;i++) {
		if (fabs(theta[i]-180.0) < ANGLE_PRECISION) {
			if (axis[i][0] < 0.0) {
				axis[i][0] = -axis[i][0];
				axis[i][1] = -axis[i][1];
				axis[i][2] = -axis[i][2];
			}
		}
	}
	
	if(do_sort != 1) return 0;
	
	/* (3) sorting the rotation angles and axes */
	for(i=0;i<60;i++) {
		val[i] = theta[i]*100. + axis[i][0]*100.0 + axis[i][1]*1.0 + axis[i][2]*0.01;
	}
	sort_double(60, val, isort);
	for(i=0;i<60;i++) {
		work_theta[i]   = theta[isort[i]];
		work_axis[i][0] = axis[isort[i]][0];
		work_axis[i][1] = axis[isort[i]][1];
		work_axis[i][2] = axis[isort[i]][2];
	}
	for(i=0;i<60;i++) {
		theta[i]   = work_theta[i];
		axis[i][0] = work_axis[i][0];
		axis[i][1] = work_axis[i][1];
		axis[i][2] = work_axis[i][2];
	}
	
	return 0;
}
