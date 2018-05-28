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
#include <string.h>

void read_biomt_file(char *filename,double *rot_array,double *xyzcent_t_array) {
	int record_finish,finish;
	int i,j,k,l,nrot;
	int ch;
	char line[82];
    int nbiomt_t;
    double rot[60][3][3];
    double xyzcent_t[60][3];
    FILE *fp;
    
    /* open the PDB file */
    if ((fp = fopen(filename, "r+")) == NULL) {
        printf("file open error\n");
        exit(1);
    }
    
	/* read lines in the PDB file */
	nbiomt_t = 0;
	finish  = 0;
	while(fgets(line, 81, fp) != NULL) {
		if ((strncmp(line,"REMARK",6) == 0) ) {
			if((strncmp(line+10,"   BIOMT",8) == 0)) {
				/* eading the rotation matrix */
				if (nbiomt_t<180) {
					nrot = (nbiomt_t) / 3;
					i    = (nbiomt_t) % 3;
					//01234567890123456789012012345678901234567890123456789012345678901234
					//REMARK 350   BIOMT2  11  0.809017  0.309017 -0.500000     -208.52413
					sscanf(line+24, "%10lf%10lf%10lf%15lf",&rot[nrot][i][0],&rot[nrot][i][1],&rot[nrot][i][2],&xyzcent_t[nrot][i]);
				}
				nbiomt_t++;
			}
		} 
	}
    fclose (fp);
    
    l = 0;
    for (i = 0; i < 60; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                rot_array[l] = rot[i][j][k];
                l++;
            }
        }
    }
    
    l = 0;
    for (i = 0; i < 60; i++) {
        for (j = 0; j < 3; j++) {
            xyzcent_t_array[l] = xyzcent_t[i][j];
            l++;
        }
    }
}



void gen_entire(double *rot_array,double *xyzcent_t_array,double *xyz_x,double *xyz_y,double *xyz_z,double *xyz_x2,double *xyz_y2,double *xyz_z2,int n,int sn) {
    
    int i,j,k,l;
    double rot[60][3][3];
    double xyzcent_t[60][3];
    double xyz_rotated[3];
    l = 0;
    
    for (i = 0; i < 60; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                rot[i][j][k] = rot_array[l];
                l++;
            }
        }
    }
    
    l = 0;
    
    for (i = 0; i < 60; i++) {
        for (j = 0; j < 3; j++) {
            xyzcent_t[i][j] = xyzcent_t_array[l];
            l++;
        }
    }
    
    for (i = 0; i < n; i++) {
        xyz_rotated[0] = rot[sn][0][0]*(xyz_x[i]) + rot[sn][0][1]*(xyz_y[i]) + rot[sn][0][2]*(xyz_z[i]) + xyzcent_t[sn][0];
        xyz_rotated[1] = rot[sn][1][0]*(xyz_x[i]) + rot[sn][1][1]*(xyz_y[i]) + rot[sn][1][2]*(xyz_z[i]) + xyzcent_t[sn][1];
        xyz_rotated[2] = rot[sn][2][0]*(xyz_x[i]) + rot[sn][2][1]*(xyz_y[i]) + rot[sn][2][2]*(xyz_z[i]) + xyzcent_t[sn][2];
        xyz_x2[i] = xyz_rotated[0];
        xyz_y2[i] = xyz_rotated[1];
        xyz_z2[i] = xyz_rotated[2];
    }
}
