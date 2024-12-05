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
 *             CAPLIB - CAPsid LIBrary
 *
 *  For calculations on icosahedrally symmetric virus capsid structures
 *
 *  To analyse directions of rotation axes, calculate cell numbers,
 *  and generate the entire structure of capsid from protomer struture
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, version 2.
 *
 *    Shigetaka Yoneda, Febrary, 2018.
 */

#include <stdio.h>
#include <math.h>
#include "planes.h"

#ifndef DEBUG
#define DEBUG 0
#endif

void set_planes() {	
	double lambda  = (1.0+sqrt(5.0))/4.0;
	double zero    = 0.0;
	double half    = 1.0/2.0;
	double one     = 1.0;
	double three   = 3.0;
	double six     = 6.0;
	double quarter = 1.0/4.0;
	double eighth  = 1.0/8.0;
	double root5   = sqrt(5.0);
	double root3   = sqrt(3.0);
	
	// variables used for debugging
	int i;
	double d0,d1,d2,d;
	double e0,e1,e2,e,emin[5];
	int j,jmin[5],k,l;
	int debug = 5;
	
	rsbc_planes.axis2[0][0] = -lambda;
	rsbc_planes.axis2[0][1] = zero;
	rsbc_planes.axis2[0][2] = zero;

	rsbc_planes.axis2[1][0] = -(one+root5)*eighth;
	rsbc_planes.axis2[1][1] = quarter;
	rsbc_planes.axis2[1][2] = (three+root5)*eighth;

	rsbc_planes.axis2[2][0] = quarter;
	rsbc_planes.axis2[2][1] = (three+root5)*eighth;
	rsbc_planes.axis2[2][2] = (one+root5)*eighth;

	rsbc_planes.axis2[3][0] = quarter;
	rsbc_planes.axis2[3][1] = (three+root5)*eighth;
	rsbc_planes.axis2[3][2] = -(one+root5)*eighth;

	rsbc_planes.axis2[4][0] = -(one+root5)*eighth;
	rsbc_planes.axis2[4][1] = quarter;
	rsbc_planes.axis2[4][2] = -(three+root5)*eighth;

	rsbc_planes.axis2[5][0] = (three+root5)*eighth;
	rsbc_planes.axis2[5][1] = (one+root5)*eighth;
	rsbc_planes.axis2[5][2] = -quarter;

	rsbc_planes.axis2[6][0] = zero;
	rsbc_planes.axis2[6][1] = zero;
	rsbc_planes.axis2[6][2] = -lambda;

	rsbc_planes.axis2[7][0] = -(three+root5)*eighth;
	rsbc_planes.axis2[7][1] = -(one+root5)*eighth;
	rsbc_planes.axis2[7][2] = -quarter;

	rsbc_planes.axis2[8][0] = -(one+root5)*eighth;
	rsbc_planes.axis2[8][1] = -quarter;
	rsbc_planes.axis2[8][2] = (three+root5)*eighth;

	rsbc_planes.axis2[9][0] = (one+root5)*eighth;
	rsbc_planes.axis2[9][1] = quarter;
	rsbc_planes.axis2[9][2] = (three+root5)*eighth;

	rsbc_planes.axis2[10][0] = -quarter;
	rsbc_planes.axis2[10][1] = (three+root5)*eighth;
	rsbc_planes.axis2[10][2] = (one+root5)*eighth;

	rsbc_planes.axis2[11][0] = -(three+root5)*eighth;
	rsbc_planes.axis2[11][1] = (one+root5)*eighth;
	rsbc_planes.axis2[11][2] = -quarter;

	rsbc_planes.axis2[12][0] = zero;
	rsbc_planes.axis2[12][1] = lambda;
	rsbc_planes.axis2[12][2] = zero;

	rsbc_planes.axis2[13][0] = -(three+root5)*eighth;
	rsbc_planes.axis2[13][1] = (one+root5)*eighth;
	rsbc_planes.axis2[13][2] = quarter;

	rsbc_planes.axis2[14][0] = -quarter;
	rsbc_planes.axis2[14][1] = (three+root5)*eighth;
	rsbc_planes.axis2[14][2] = -(one+root5)*eighth;

	for(int i=0;i<15;i++){
		rsbc_planes.axis2[15+i][0] = -rsbc_planes.axis2[i][0];
		rsbc_planes.axis2[15+i][1] = -rsbc_planes.axis2[i][1];
		rsbc_planes.axis2[15+i][2] = -rsbc_planes.axis2[i][2];
	}
	
	if(DEBUG > 4) {
		for (i=0; i<15; i++) {
			printf("2fold axis %d:%8.4f%8.4f%8.4f\n",i,
			rsbc_planes.axis2[i][0],
			rsbc_planes.axis2[i][1],
			rsbc_planes.axis2[i][2]);
		}
	}
	if(DEBUG == 1) {
		// The following debugging lines were used to confirm the lengths of the 2fold axes are equal to 0.809 
		for (i=0; i< 15; i++) {
			d0 = rsbc_planes.axis2[i][0];
			d1 = rsbc_planes.axis2[i][1];
			d2 = rsbc_planes.axis2[i][2];
			d = d0*d0 + d1*d1 +d2*d2;
			printf("%8.4f ",sqrt(d));
		}
		printf("\n");
		
		// Following debugging lines are used to confirm the 2fold axis vectors have the 4 most neighbor ones each. */
		for (i=0; i< 15; i++) {
			emin[0] = -1.0e30;
			emin[1] = -1.0e30;
			emin[2] = -1.0e30;
			emin[3] = -1.0e30;
			jmin[0] = -1;
			jmin[1] = -1;
			jmin[2] = -1;
			jmin[3] = -1;
			d0 = rsbc_planes.axis2[i][0];
			d1 = rsbc_planes.axis2[i][1];
			d2 = rsbc_planes.axis2[i][2];
			
			for (j=0; j< 15; j++) {
				if(j==i) {continue;}
				e0 = rsbc_planes.axis2[j][0] * d0;
				e1 = rsbc_planes.axis2[j][1] * d1;
				e2 = rsbc_planes.axis2[j][2] * d2;
				e = e0 + e1 +e2;
				if(e<0.0) e = -e;
				for (k=0; k<5; k++) {
					if(e > emin[k]) {
						for (l=3; l>=k; l--){
							emin[l+1] = emin[l];
							jmin[l+1] = jmin[l];
						}
						emin[k] = e;
						jmin[k] = j;
						break;
					}
				}
			}
			printf("closest for %d = ",i);
			for (j=0;j<5;j++){
				printf(" %8.4f %d",emin[j],jmin[j]);
			}
			printf("\n");	
		}
	}

	//3fold rotation axes
	//abc 5-12-1 4-11-0
	rsbc_planes.axis3[0][0] = zero;
	rsbc_planes.axis3[0][1] = (one + root5)*root3/six;
	rsbc_planes.axis3[0][2] = (root5 - one)*root3/six;
	//abc' 5-12-11  4-11-10
	rsbc_planes.axis3[1][0] = zero;
	rsbc_planes.axis3[1][1] = (one + root5)*root3/six;
	rsbc_planes.axis3[1][2] = -(root5 - one)*root3/six;
	// def 2-9-8
	rsbc_planes.axis3[2][0] = zero;
	rsbc_planes.axis3[2][1] = -(one + root5)*root3/six;
	rsbc_planes.axis3[2][2] = (root5 - one)*root3/six;
	//d'ef 7-9-8
	rsbc_planes.axis3[3][0] = zero; 
	rsbc_planes.axis3[3][1] = -(one + root5)*root3/six;
	rsbc_planes.axis3[3][2] = -(root5 - one)*root3/six;
	//cdh 1-2-3
	rsbc_planes.axis3[4][0] = (root5 - one)*root3/six;
	rsbc_planes.axis3[4][1] = zero;
	rsbc_planes.axis3[4][2] = (one + root5)*root3/six;
	//cdg 1-2-4
	rsbc_planes.axis3[5][0] = -(root5 - one)*root3/six;
	rsbc_planes.axis3[5][1] = zero;
	rsbc_planes.axis3[5][2] =  (one + root5)*root3/six;
	//c'd'h' 11-7-10
	rsbc_planes.axis3[6][0] = (root5 - one)*root3/six;
	rsbc_planes.axis3[6][1] = zero;
	rsbc_planes.axis3[6][2] = -(one + root5)*root3/six;
	//c'd'g' 11-7-6
	rsbc_planes.axis3[7][0] = -(root5 - one)*root3/six;
	rsbc_planes.axis3[7][1] = zero;
	rsbc_planes.axis3[7][2] = -(one + root5)*root3/six;
	//agg' 5-4-6
	rsbc_planes.axis3[8][0] = -(one + root5)*root3/six; 
	rsbc_planes.axis3[8][1] = (root5 - one)*root3/six;
	rsbc_planes.axis3[8][2] = zero;
	//bhh' 12-3-10
	rsbc_planes.axis3[9][0] = (one + root5)*root3/six; 
	rsbc_planes.axis3[9][1] = (root5 - one)*root3/six;
	rsbc_planes.axis3[9][2] = zero;
	//egg' 9-4-6
	rsbc_planes.axis3[10][0] = -(one + root5)*root3/six; 
	rsbc_planes.axis3[10][1] = -(root5 - one)*root3/six;
	rsbc_planes.axis3[10][2] = zero;
	//fhh' 8-3-10
	rsbc_planes.axis3[11][0] = (one + root5)*root3/six; 
	rsbc_planes.axis3[11][1] = -(root5 - one)*root3/six;
	rsbc_planes.axis3[11][2] = zero;
	//acg 5-1-4
	rsbc_planes.axis3[12][0] = -one/root3;
	rsbc_planes.axis3[12][1] = one/root3;
	rsbc_planes.axis3[12][2] = one/root3;
	//bch 12-1-3
	rsbc_planes.axis3[13][0] = one/root3;
	rsbc_planes.axis3[13][1] = one/root3;
	rsbc_planes.axis3[13][2] = one/root3;
	//deg 2-9-4
	rsbc_planes.axis3[14][0] = -one/root3;
	rsbc_planes.axis3[14][1] = -one/root3;
	rsbc_planes.axis3[14][2] = one/root3;
	//dfh 2-8-3  <--- This line was previously 2-8-4 (error). Thus corrected.
	rsbc_planes.axis3[15][0] = one/root3;
	rsbc_planes.axis3[15][1] = -one/root3;
	rsbc_planes.axis3[15][2] = one/root3;
	//ac'g' 5-11-6
	rsbc_planes.axis3[16][0] = -one/root3;
	rsbc_planes.axis3[16][1] = one/root3;
	rsbc_planes.axis3[16][2] = -one/root3;
	//bc'h'	12-11-10
	rsbc_planes.axis3[17][0] = one/root3;
	rsbc_planes.axis3[17][1] = one/root3;
	rsbc_planes.axis3[17][2] = -one/root3;
	//d'eg' 7-9-6
	rsbc_planes.axis3[18][0] = -one/root3;
	rsbc_planes.axis3[18][1] = -one/root3;
	rsbc_planes.axis3[18][2] = -one/root3;
	//d'fh' 7-8-10
	rsbc_planes.axis3[19][0] = one/root3;
	rsbc_planes.axis3[19][1] = -one/root3;
	rsbc_planes.axis3[19][2] = -one/root3;
	
	if(DEBUG == 2) {
		// The following debugging lines were used to confirm the lengths of the 3fold axes are equal to 1.0 
		for (i=0; i< 20; i++) {
			d0 = rsbc_planes.axis3[i][0];
			d1 = rsbc_planes.axis3[i][1];
			d2 = rsbc_planes.axis3[i][2];
			d = d0*d0 + d1*d1 +d2*d2;
			printf("%8.4f ",sqrt(d));
		}
		printf("\n");
		
		// Following debugging lines are used to confirm the 3fold axis vectors have the 4 most neighbor ones each. */
		for (i=0; i< 20; i++) {
			emin[0] = -1.0e30;
			emin[1] = -1.0e30;
			emin[2] = -1.0e30;
			emin[3] = -1.0e30;
			jmin[0] = -1;
			jmin[1] = -1;
			jmin[2] = -1;
			jmin[3] = -1;
			d0 = rsbc_planes.axis3[i][0];
			d1 = rsbc_planes.axis3[i][1];
			d2 = rsbc_planes.axis3[i][2];
			
			for (j=0; j< 20; j++) {
				if(j==i) {continue;}
				e0 = rsbc_planes.axis3[j][0] * d0;
				e1 = rsbc_planes.axis3[j][1] * d1;
				e2 = rsbc_planes.axis3[j][2] * d2;
				e = e0 + e1 +e2;
				if(e<0.0) e = -e;
				for (k=0; k<5; k++) {
					if(e > emin[k]) {
						for (l=3; l>=k; l--){
							emin[l+1] = emin[l];
							jmin[l+1] = jmin[l];
						}
						emin[k] = e;
						jmin[k] = j;
						break;
					}
				}
			}
			printf("closest for %d = ",i);
			for (j=0;j<5;j++){
				printf(" %8.4f %d",emin[j],jmin[j]);
			}
			printf("\n");	
		}
	}
	
	//5fold rotation axis
	//c 1 0
	rsbc_planes.axis5[0][0] = zero;
	rsbc_planes.axis5[0][1] = half;
	rsbc_planes.axis5[0][2] = lambda;
	//d 2 1
	rsbc_planes.axis5[1][0] = zero;
	rsbc_planes.axis5[1][1] = -half;
	rsbc_planes.axis5[1][2] = lambda;	
	//h 3 2
	rsbc_planes.axis5[2][0] = lambda;
	rsbc_planes.axis5[2][1] = zero;
	rsbc_planes.axis5[2][2] = half;
	//g 4 3
	rsbc_planes.axis5[3][0] = -lambda;
	rsbc_planes.axis5[3][1] = zero;
	rsbc_planes.axis5[3][2] = half;
	//a 5 4
	rsbc_planes.axis5[4][0] = -half;
	rsbc_planes.axis5[4][1] = lambda;
	rsbc_planes.axis5[4][2] = zero;
	//g' 6 5
	rsbc_planes.axis5[5][0] = -lambda;
	rsbc_planes.axis5[5][1] = zero;
	rsbc_planes.axis5[5][2] = -half;
	//d' 7 6
	rsbc_planes.axis5[6][0] = zero;
	rsbc_planes.axis5[6][1] = -half;
	rsbc_planes.axis5[6][2] = -lambda;
	//f 8 7
	rsbc_planes.axis5[7][0] = half;
	rsbc_planes.axis5[7][1] = -lambda;
	rsbc_planes.axis5[7][2] = zero;
	//e 9 8
	rsbc_planes.axis5[8][0] = -half;
	rsbc_planes.axis5[8][1] = -lambda;
	rsbc_planes.axis5[8][2] = zero;
	//h' 10 9
	rsbc_planes.axis5[9][0] = lambda;
	rsbc_planes.axis5[9][1] = zero;
	rsbc_planes.axis5[9][2] = -half;	
	//c' 11 10
	rsbc_planes.axis5[10][0] = zero;
	rsbc_planes.axis5[10][1] = half;
	rsbc_planes.axis5[10][2] = -lambda;
	//b 12 11
	rsbc_planes.axis5[11][0] = half;
	rsbc_planes.axis5[11][1] = lambda;
	rsbc_planes.axis5[11][2] = zero;
	
	if(DEBUG == 3) {
		// The following debugging lines were used to confirm the lengths of the 5fold axes are equal to 0.951
		for (i=0; i< 12; i++) {
			d0 = rsbc_planes.axis5[i][0];
			d1 = rsbc_planes.axis5[i][1];
			d2 = rsbc_planes.axis5[i][2];
			d = d0*d0 + d1*d1 +d2*d2;
			printf("%8.4f ",sqrt(d));
		}
		printf("\n");
		
		// Following debugging lines are used to confirm the 5fold axis vectors have the 4 most neighbor ones each. */
		for (i=0; i< 12; i++) {
			emin[0] = -1.0e30;
			emin[1] = -1.0e30;
			emin[2] = -1.0e30;
			emin[3] = -1.0e30;
			jmin[0] = -1;
			jmin[1] = -1;
			jmin[2] = -1;
			jmin[3] = -1;
			d0 = rsbc_planes.axis5[i][0];
			d1 = rsbc_planes.axis5[i][1];
			d2 = rsbc_planes.axis5[i][2];
			
			for (j=0; j< 20; j++) {
				if(j==i) {continue;}
				e0 = rsbc_planes.axis5[j][0] * d0;
				e1 = rsbc_planes.axis5[j][1] * d1;
				e2 = rsbc_planes.axis5[j][2] * d2;
				e = e0 + e1 +e2;
				if(e<0.0) e = -e;
				for (k=0; k<5; k++) {
					if(e > emin[k]) {
						for (l=3; l>=k; l--){
							emin[l+1] = emin[l];
							jmin[l+1] = jmin[l];
						}
						emin[k] = e;
						jmin[k] = j;
						break;
					}
				}
			}
			printf("closest for %d = ",i);
			for (j=0;j<5;j++){
				printf(" %8.4f %d",emin[j],jmin[j]);
			}
			printf("\n");	
		}
	}
}
