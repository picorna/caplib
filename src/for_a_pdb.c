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
 *  forapdb.c         calculation for the 60 matrices of a PDB file 
 *  capaxis
 *
 *  Created by yoneda on 11/04/22.
 *  Copyright 2011 Shigetaka Yoneda. All rights reserved.
 *
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#ifndef UNITMATRIX_LIMIT
#define UNITMATRIX_LIMIT 0.001
#endif

int line80_read(FILE *,char [],char,int *);
double eqsolv(double [3][3],double [3]);
void sort_double(int, double [], int i[]);
void error_store(int,char *);

int for_a_pdb(FILE *fin,int nfiles,double theta[60],double axis[60][3],char id[4]){
	double xyz[3],xyzcent[3];
	double rot[60][3][3];
	double val[60],work_theta[60],work_axis[60][3];
	int isort[60];
	int nbiomt,natom,nrot,final,length;
	int i,n;
	double determinant;
	char line[81];
	char char_in_prev_line;
	int should_return;
	char e_message1[81],e_message2[81];
		
	/* (0) reading lines in the PDB file */
	nbiomt = 0;
	natom  = 0;
	final = 0;
	char_in_prev_line = '\n';
	id[0] = ' ';
	id[1] = ' ';
	id[2] = ' ';
	id[3] = ' ';
	nrot = 0;
	for(;final==0;) {
		/* reading a line in the PDB file */
		if(fgets(line, 81, fin) == NULL) {
			final = 1;
		} else {
			if (strchr(line,'\n') != NULL) {
				length = strlen(line);
				line[length - 1] = '\0';
			} else {
				while (fgetc(fin) != '\n');
			}
		}

		/* for an ATOM line */
		if (strncmp(line,"ATOM",4)==0){
			/* reading the XYZ coordinates of the ATOM lines (HETATM excluded) */
			n = sscanf(line+30, "%8lf%8lf%8lf", &xyz[0], &xyz[1], &xyz[2]);
			xyzcent[0] += xyz[0];
			xyzcent[1] += xyz[1];
			xyzcent[2] += xyz[2];
			natom++;
		/* for a BIOMT line */
		} else if (line[0]=='R' && line[1]=='E' && line[2]=='M'
				   && line[13]=='B' && line[14]=='I'
				   && line[15]=='O' && line[16]=='M' && line[17]=='T'
				   && (line[18]=='1' || line[18]=='2' || line[18]=='3')) {
			/* eading the rotation matrix */
			if (nbiomt<180) {
				nrot = nbiomt / 3;
				i    = nbiomt % 3;
				n = sscanf(line+24, "%10lf%10lf%10lf", &rot[nrot][i][0], &rot[nrot][i][1], &rot[nrot][i][2]);
				nbiomt++;
			}
		} else if (strncmp(line,"HEADER",6)==0) {
			//123456789012345678901234567890123456789012345678901234567890123
			//HEADER    VIRUS                                   25-JAN-88   4RHV 
			id[0] = line[62];
			id[1] = line[63];
			id[2] = line[64];
			id[3] = line[65];
		}
	}
	printf("%3d %c%c%c%c\n",nfiles+1,id[0],id[1],id[2],id[3]);
	
	if (natom == 0) {
		error_store(nfiles,"No atoms in this PDB file");
	} else {
		xyzcent[0] /= natom;
		xyzcent[1] /= natom;
		xyzcent[2] /= natom;
	}

	if(nbiomt != 60*3) {
		sprintf(e_message1,"The number of REMARK BIOMT lines is %d (not 60*3)",nbiomt);
		error_store(nfiles,e_message1);
		return 1;
	}

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
			sprintf(e_message2,"Determinant of a BIOMT matrix is %6.2f (not 1.0)",determinant);
			error_store(nfiles,e_message2);
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
			theta[i] = fabs(eqsolv((double (*)[3])&rot[i][0][0],(double (*))&axis[i][0]));
		}
	}

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
		theta[i]          = work_theta[i];
		axis[i][0] = work_axis[i][0];
		axis[i][1] = work_axis[i][1];
		axis[i][2] = work_axis[i][2];
	}

	return 0;
}
