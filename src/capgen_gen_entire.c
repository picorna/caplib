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
#include <string.h>

int cell_number(double *);

/* Memo(SY):Although this function was firstly made to genearte the entire structure 
 by repeating 60 rotations of BIOMET lines in a PDB file, this is used now to print
 the PDB file with and without rotations, adding some informations on BIOMT, etc, . */

void gen_entire(int cells[60],int f_verbose,
        int nlines,char *pdb_line,double *pdb_xyz,
        double rot[60][3][3],double xyzcent[60][3]) {
    
    int i,j,k;
    double *x,*y,*z;
    double xyz_rotated[3];
    char *line;
    int first_struct,first_atom;
    int nrot,nbiomt;
    int cell,natoms60;
    int make_entire=1;
    int get_center_of_entire=1;
    
    first_struct = 1;
    nbiomt   = 0;
    natoms60 = 0;
    for (i=0; i<60; i++) {
        if(cells[i] == 1){
        line = pdb_line;
        x    = pdb_xyz;
        first_atom = 1;
        for (j=0; j<nlines; j++) {
            if (strncmp(line,"ATOM",4) == 0 || strncmp(line,"HETATM",6) == 0){
                natoms60++;
                if (first_atom == 1) {
                    printf("MODEL     %4d\n",i);
                    first_atom = 0;
                }
                y = x+1;
                z = x+2;
                xyz_rotated[0] =
                    rot[i][0][0]*(*x) + rot[i][0][1]*(*y) + rot[i][0][2]*(*z) + xyzcent[i][0];
                xyz_rotated[1] =
                    rot[i][1][0]*(*x) + rot[i][1][1]*(*y) + rot[i][1][2]*(*z) + xyzcent[i][1];
                xyz_rotated[2] =
                    rot[i][2][0]*(*x) + rot[i][2][1]*(*y) + rot[i][2][2]*(*z) + xyzcent[i][2];
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
                x += 3;
            } else if (first_struct == 1) {
                if(!(*line=='E' && *(line+1)=='N' && *(line+2)=='D')) {
                    printf("%s\n",line);
                }
            }
            line += 81;
        }
        printf("ENDMDL    %4d\n",i);
        first_struct = 0;
        }
    }
    printf("END\n");
}
