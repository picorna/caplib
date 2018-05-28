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

void datashift(int [], int, int);

void sort_double(int m, double x[], int isort[]) {
/*
     Sorting double-precision data in the ascending order
     This program was written by myself, Shigetaka Yoneda, 1983.Feb.09   

     x is the data to be sorted. But x is not changed in this function.
     m is the number of the data. m is not changed in this function.
     isort stores the result of the calculation in this function (the ordered indices of x).
     work is a temporary storage.
*/
	double xi,xs,xe,xmid;
	int i,is,ie,imid;
	int length;

	isort[0] = 0;
	for (i = 1; i < m; i++) {
		xi = x[i];
		is = 0;
		ie = i-1;
		length = ie+1; /* This is equivalent to length = ie-is+1;*/
		xs = x[isort[is]];
		xe = x[isort[ie]];
		if(xi <= xs) {
			datashift(isort,1,i);
//		} else if(xe >= xi) {
		} else if(xe <= xi) {
			isort[i] =i;
		} else {
NEXTDATUM:
			if(length <= 2) {
				datashift(isort,ie,i);
			} else {
				imid = is+(length/2);
				xmid = x[isort[imid]];        
				if(xi < xmid) {
					ie = imid;
					xe = xmid;
					length = ie-is+1;
					goto NEXTDATUM;
				} else if(xi > xmid) {
					is = imid;
					xs = xmid;
					length = ie-is+1;
					goto NEXTDATUM;
				} else {
					datashift(isort,imid,i);
				}
			}
		}
	}
	return;
}