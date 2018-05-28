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

double matrix_determinant(double [3][3]);

void make_orthogonal_matrix(double matrix[3][3],double vec5[3],double vec3[3]){
	double uvec5[3],ovec3[3],normal[3];
	double d;
	int i;
	int debug = 0;

	d = 1.0/sqrt(vec5[0]*vec5[0] + vec5[1]*vec5[1] + vec5[2]*vec5[2]);
	uvec5[0] = vec5[0]*d;
	uvec5[1] = vec5[1]*d;
	uvec5[2] = vec5[2]*d;
	d = vec3[0]*uvec5[0] + vec3[1]*uvec5[1] + vec3[2]*uvec5[2];
	ovec3[0] = vec3[0] - uvec5[0]*d;
	ovec3[1] = vec3[1] - uvec5[1]*d;
	ovec3[2] = vec3[2] - uvec5[2]*d;
	d = 1.0/sqrt(ovec3[0]*ovec3[0] + ovec3[1]*ovec3[1] + ovec3[2]*ovec3[2]);
	ovec3[0] *= d;
	ovec3[1] *= d;
	ovec3[2] *= d;
	normal[0] = ovec3[1]*uvec5[2] - ovec3[2]*uvec5[1];
	normal[1] = ovec3[2]*uvec5[0] - ovec3[0]*uvec5[2];
	normal[2] = ovec3[0]*uvec5[1] - ovec3[1]*uvec5[0];
	for (i=0; i<3; i++) {
		matrix[i][0] = uvec5[i];
		matrix[i][1] = ovec3[i];
		matrix[i][2] = -normal[i];
	}
	if (debug == 1) {
		printf("determinant in make_orthogonal_matrix = %10.5f\n",matrix_determinant(matrix));
	}
}