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
//
// matrix_inverse.c
// to calculate the inverse matrix of a 3 * 3 matrix
//  Created by yoneda on 2014/07/30.
//
//

#include <stdio.h>
#define TOOSMALL 1.0e-30

void matrix_inverse(double a[3][3],double inv[3][3]) {
    double work[3][3];
	int i,j,i1,i2,j1,j2;
    
// The determinant of the matrix, a, is assumed to be 1.0;
// This assumption is confirmed in set_2plv_rotmat.c before calling this function.

    for(i=0;i<3;i++) {
        i1 = (i+1)%3;
        i2 = (i+2)%3;
		for(j=0;j<3;j++) {
            j1 = (j+1)%3;
            j2 = (j+2)%3;
            work[i][j] = a[i1][j1]*a[i1][j2] - a[i1][j2]*a[i2][j1];
		}
	}
	for(i=0;i<3;i++) {
		for(j=0;j<3;j++) {
           inv[i][j] = work[j][i];
		}
	}
    return;
}