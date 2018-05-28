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
//
//  load_wat216.c
//  capgen
//
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
//  This file was created by yoneda on 2014/07/23.
//
//

#include <stdio.h>
#include <stdlib.h>
#include "wat216.h"

int read_wat216(FILE *fin) {
    //
    // This function reads the widely-known WAT216.dat file of water coordinates.
    //
    int record_finish;
    int i,j;
    int finish,ladd,mol,in_a_line;
    char ch,line[80];
    
    // reading a line with length of 41 characters
    if(fgets(line, 42, fin) == NULL) {
        return 1;
    } else {
        record_finish = 0;
        for(i=0;i<41;i++) {
            if (*(line+i) == '\n') {
                *(line+i) = '\0';
                record_finish = 1;
                break;
            }
        }
        if (record_finish == 0) {
            *(line+41) = '\0';
            while ((ch = fgetc(fin)) != '\n') {
                if (ch == EOF) {
                    finish = 1;
                    break;
                }
            }
        }
    }
    // read nwater_mols and size_wat216 in %5d and %12.3f
    sscanf(line,"%5d%12lf",&wat216.nmols,&wat216.cubic_size);
    if(wat216.nmols > MAX_WAT_MOLS){
        fprintf(stderr,"Error: number of mols is %d, greater than MAX_WAT_MOLS in read_wat216.c\n",wat216.nmols);
		exit(EXIT_FAILURE);
    }
    if(wat216.nmols != 216){
        fprintf(stderr,"The number of water molecules is not the standard number of 216\n");
    }

    in_a_line = 0;
    ladd = 0;
    // Although a water molecule in wat216.dat contains 4 atoms, only the coordinate of 3 atoms are read, because the coordinates are used only for printing water atom coordinates in a pdb file.
    // wat216.xyz[3 or 4 atoms][x,y,z][216 mols]
    for(i=0;i<3;i++){
        for(j=0;j<3;j++){
            for(mol=0;mol<wat216.nmols;mol++){
                // reading a line with length of 78 characters
                if(in_a_line==0){
                    if(fgets(line, 79, fin) == NULL) {
                        return 1;
                    } else {
                        record_finish = 0;
                        for(i=0;i<78;i++) {
                            if (*(line+i) == '\n') {
                                *(line+i) = '\0';
                                record_finish = 1;
                                break;
                            }
                        }
                        if (record_finish == 0) {
                            *(line+78) = '\0';
                            while ((ch = fgetc(fin)) != '\n') {
                                if (ch == EOF) {
                                    finish = 1;
                                    break;
                                }
                            }
                        }
                    }
                    in_a_line = 1;
                    ladd = 0;
                }
                // reading data in a line in the format of 6 * f13
                sscanf(line+ladd, "%13lf",&wat216.xyz[i][j][mol]);
                ladd += 13;
                if(ladd>77){
                    in_a_line = 0;
                }
            }
        }
    }
    return 0;
}