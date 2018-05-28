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

/* Memo(SY): This function was a copy from read_a_pdb.c, and many functions were deleted.*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void read_biomt_file(FILE *fin_matrix,double rot[60][3][3],double xyzcent_t[60][3],
                int *nbiomt_t) {
    int record_finish,finish;
    int i,nrot;
    int ch;
    char line[82];
    
    /* read lines in the PDB file */
    *nbiomt_t = 0;
    finish  = 0;
    while(finish == 0) {
    // read a line in the PDB file 
    /* fgets(line, 81, fin_matrix) reads 80 characters and adds \0 at the 81th position
       When \n appears from 0 to 79, \0 is added after \n */
        if(fgets(line, 81, fin_matrix) == NULL) {
            return;
        } else {
            record_finish = 0;
            for(i=0;i<80;i++) {
                if (*(line+i) == '\n') {
                    *(line+i) = '\0';
                    record_finish = 1;
                    break;
                }
            }
            if (record_finish == 0) {
                *(line+80) = '\0';
                while ((ch = fgetc(fin_matrix)) != '\n') {
                    if (ch == EOF) {
                        finish = 1;
                        break;
                    }
                }
            }
        }
        
        if ((strncmp(line,"REMARK",6) == 0) ) {
            if((strncmp(line+10,"   BIOMT",8) == 0)) {
                /* eading the rotation matrix */
                if (*nbiomt_t<180) {
                    nrot = (*nbiomt_t) / 3;
                    i    = (*nbiomt_t) % 3;
                    //01234567890123456789012012345678901234567890123456789012345678901234
                    //REMARK 350   BIOMT2  11  0.809017  0.309017 -0.500000     -208.52413
                    sscanf(line+24, "%10lf%10lf%10lf%15lf",
                               &rot[nrot][i][0],&rot[nrot][i][1],&rot[nrot][i][2],
                               &xyzcent_t[nrot][i]);
                }
                (*nbiomt_t)++;
            }
        } 
    }
}
