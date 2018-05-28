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

#include "summary_store.h"

void print_rot(int i,int m10,int m10end) {
	int k,n,nini,nend,nprint;
	
	nprint = m10end-m10;
	for (k=0; k<(nprint-1)/10+1; k++) {
		nini = m10 + 10*k;
		nend = nini + 10;
		if (nend>m10end) {
			nend = m10end;
		}
		printf("angle = ");
		for (n=nini; n<nend; n++) {
			printf("%7.1f",store[i].theta[n]);
		}
		printf("\n");
		printf("x     = ");
		for (n=nini; n<nend; n++) {
			printf("%7.3f",store[i].axis[n][0]);
		}
		printf("\n");
		printf("y     = ");
		for (n=nini; n<nend; n++) {
			printf("%7.3f",store[i].axis[n][1]);
		}
		printf("\n");
		printf("z     = ");
		for (n=nini; n<nend; n++) {
			printf("%7.3f",store[i].axis[n][2]);
		}
		printf("\n");			
	}
}