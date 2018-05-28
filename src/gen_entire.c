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
 *	CAPLIB - calculations on icosahedrally symmetric virus capsid structures
 *            analyse directions of rotation axes, calculate cell numbers, and generate the entire
 *            structure of a capsid from a protomer struture
 *  copyright 2004, 2011, Shigetaka Yoneda
 *
 */

#include <stdio.h>
#include <string.h>

int cell_number(double *);

/* Although this function was firstly made to genearte the entire structure by repeating 60 rotations 
 of BIOMET lines in a PDB file, this function is used now to print the PDB file with and without rotations, 
 adding some informations on BIOMT, etc, . */

void gen_entire(int f_verbose,int make_entire,int get_center_of_entire,double xyz_cent_entire[3],
				 int nlines,char *pdb_line,double *pdb_xyz,double rot[60][3][3],double xyzcent[60][3],
				 double axis[60][3],double theta[60]) {
	
	int i,j,k;
	double *x,*y,*z;
	double xyz_rotated[3];
	char *line;
	int first_struct;
	int nrot,nbiomt;
	int cell,natoms60;
	
	//printf("f_verbose = %d\n",f_verbose);
	
	// When get_center_of_entire == 0, the center of the entire capsid is calculated, but no printing.
	first_struct = 1;
	nbiomt   = 0;
	natoms60 = 0;
	xyz_cent_entire[0] = 0.;
	xyz_cent_entire[1] = 0.;
	xyz_cent_entire[2] = 0.;
	for (i=0; i<60; i++) {
		//printf("******* %d ****\n",i);
		line = pdb_line;
		x    = pdb_xyz;
		for (j=0; j<nlines; j++) {
			if (strncmp(line,"ATOM",4) == 0 || strncmp(line,"HETATM",6) == 0){
				natoms60++;
				if ((make_entire == 1) && (first_struct == 1) && (get_center_of_entire != 1)) {
					printf("MODEL     %4d\n",i+1);
					first_struct = 0;
				}
				y = x+1;
				z = x+2;
				if(make_entire == 1) {
					xyz_rotated[0] = rot[i][0][0]*(*x) + rot[i][0][1]*(*y) + rot[i][0][2]*(*z) + xyzcent[i][0];
					xyz_rotated[1] = rot[i][1][0]*(*x) + rot[i][1][1]*(*y) + rot[i][1][2]*(*z) + xyzcent[i][1];
					xyz_rotated[2] = rot[i][2][0]*(*x) + rot[i][2][1]*(*y) + rot[i][2][2]*(*z) + xyzcent[i][2];
				} else {
					xyz_rotated[0] = *x;
					xyz_rotated[1] = *y;
					xyz_rotated[2] = *z;
				}
				if (get_center_of_entire == 1) {
					xyz_cent_entire[0] += xyz_rotated[0];				
					xyz_cent_entire[1] += xyz_rotated[1];				
					xyz_cent_entire[2] += xyz_rotated[2];
				} else {
					if (make_entire == 1) {
						if(f_verbose == 1) {
							cell = cell_number(xyz_rotated);
							for (k=0; k<30; k++) {
								printf("%c",*(line+k));
							}
							if (strlen(line)>54) {
								printf("%8.3f%8.3f%8.3f",
									   xyz_rotated[0],xyz_rotated[1],xyz_rotated[2]);
								for (k=54; k<60; k++) {
									printf("%c",*(line+k));
								}
								printf("%3d.00",cell);
								for (k=66; k<80; k++) {
									printf("%c",*(line+k));
								}
								printf("\n");
							} else {
								printf("%8.3f%8.3f%8.3f                      %3d\n",
									   xyz_rotated[0],xyz_rotated[1],xyz_rotated[2],cell);
							}
						} else {
							for (k=0; k<30; k++) {
								printf("%c",*(line+k));
							}
							if (strlen(line)>54) {
								printf("%8.3f%8.3f%8.3f%s\n",
									   xyz_rotated[0],xyz_rotated[1],xyz_rotated[2],line+54);
							} else {
								printf("%8.3f%8.3f%8.3f\n",
									   xyz_rotated[0],xyz_rotated[1],xyz_rotated[2]);
							}					
						}
					} else {
						if(f_verbose == 1) {
							cell = cell_number(xyz_rotated);
							for (k=0; k<60; k++) {
								printf("%c",*(line+k));
							}
							printf("%6.2f",(double) cell);
							printf("%s\n",line+66);
							
						} else {
							printf("%s\n",line);
						}
					}
				}
				x += 3;
			} else if (i==0) {
				if (get_center_of_entire != 1) {
					if ((f_verbose == 1) && (strncmp(line,"REMARK",6) == 0) && (strncmp(line+10,"   BIOMT",8) == 0)) {
						for (k=0; k<24; k++) {
							printf("%c",*(line+k));
						}
						/* printing the axes and angles for the rotation matrices at the BIOMT lines*/
						nrot = nbiomt/3;
						k    = nbiomt%3;
						//01234567890123456789012012345678901234567890123456789012345678901234
						//REMARK 350   BIOMT2  11  0.809017  0.309017 -0.500000     -208.52413
						if(k==0) {
							printf("%7.3f%7.3f%7.3f%7.3f axis=%6.3f angle=%6.1f\n",
								   rot[nrot][k][0],rot[nrot][k][1],rot[nrot][k][2],xyzcent[nrot][k],
								   axis[nrot][k],theta[nrot]);
						} else {
							printf("%7.3f%7.3f%7.3f%7.3f axis=%6.3f\n",
								   rot[nrot][k][0],rot[nrot][k][1],rot[nrot][k][2],xyzcent[nrot][k],
								   axis[nrot][k]);
						}
						nbiomt++;
					} else if ((make_entire != 1) || (j < nlines-1 || strncmp(line,"END",3) != 0)) {
						printf("%s\n",line);
					}
				}
			}
			line += 81;
		}
		if(make_entire != 1 && get_center_of_entire != 1) {
			break;
		}
		if(get_center_of_entire != 1) {
			printf("ENDMDL\n");
			if (i == 59) {
				printf("END\n");
			} else {
				printf("MODEL     %4d\n",i+2);
			}
		}
	}
	if (get_center_of_entire == 1) {
		xyz_cent_entire[0] /= natoms60;
		xyz_cent_entire[1] /= natoms60;
		xyz_cent_entire[2] /= natoms60;
		//printf("%8.3f %8.3f %8.3f %d\n",xyz_cent_entire[0],xyz_cent_entire[1],xyz_cent_entire[2],natoms60);
	}
}
