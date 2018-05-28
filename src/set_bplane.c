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

#include "wat216.h"
#include "planes.h"

void set_bplane(vodi){
    int i;
    
    // Another notation of planes and axes for cell 0 (1/60 and 1/12)
    // Border plane for cell 0 (1/60) are the 2fold axes, 0 (down), 3 (down), 12 (up), 13 (up)
    if(wat216.type_of_rsbc == 60) {
        for(i=0;i<3;i++) {
            // The cell 0 (1/60) has 4 side planes: 2 of them are wide (long), and the other 2 are narrow (short)
            rsbc_planes.bplane[0][i] = -rsbc_planes.axis2[0][i]; // long plane
            rsbc_planes.bplane[1][i] = -rsbc_planes.axis2[0][i]; // long plane
            rsbc_planes.bplane[2][i] =  rsbc_planes.axis2[12][i]; // short plane
            rsbc_planes.bplane[3][i] =  rsbc_planes.axis2[12][i]; // short plane
        }
        
        // In the previous program, apricot, bcenter[][] was subtracted in cell_number.c and front.c.
        // However, I try to omit the subtraction.
        /*
         for(k=0;k<4;k++){
         k1 = (k+1) % 4;
         for(i=0;i<3;i++){
         rsbc_planes.bcenter[k][i] = rsbc_planes.bplane[k][i] + rsbc_planes.bplane[k1][i];
         }
         d = rsbc_planes.bcenter[k][0]*rsbc_planes.bcenter[k][0]
         + rsbc_planes.bcenter[k][1]*rsbc_planes.bcenter[k][1] + rsbc_planes.bcenter[k][2]*rsbc_planes.bcenter[k][2];
         d = 0.5/sqrt(d);
         for(i=0;i<3;i++){
         rsbc_planes.bcenter[k][i] *= d*(wat216.inner + wat216.outer);
         }
         }*/
    } else {
        // Border plane for cell 0 (1/12) are the 2fold axes, 13 (up), 12 (up), 8 (up), 7 (down), 4 (down)
        for(i=0;i<3;i++) {
            rsbc_planes.bplane12[0][i] =  rsbc_planes.axis2[13][i]; // for cell 0 and 1 (1/60)
            rsbc_planes.bplane12[1][i] =  rsbc_planes.axis2[8][i]; // for cell 1 and 2 (1/60)
            rsbc_planes.bplane12[2][i] = -rsbc_planes.axis2[4][i]; // for cell 2 and 3 (1/60)
            rsbc_planes.bplane12[3][i] = -rsbc_planes.axis2[7][i]; // for cell 3 and 4 (1/60)
            rsbc_planes.bplane12[4][i] =  rsbc_planes.axis2[12][i]; // for cell 0 and 4 (1/60)
        }
        // In the previous program, apricot, bcenter[][] was subtracted in cell_number.c and front.c.
        // However, I try to omit the subtraction.
        /*
         for(k=0;k<5;k++){
         k1 = (k+1) % 5;
         for(i=0;i<3;i++){
         rsbc_planes.bcenter12[k][i] = rsbc_planes.bplane12[k][i] + rsbc_planes.bplane12[k1][i];
         }
         d = rsbc_planes.bcenter12[k][0]*rsbc_planes.bcenter12[k][0]
         + rsbc_planes.bcenter12[k][1]*rsbc_planes.bcenter12[k][1] + rsbc_planes.bcenter12[k][2]*rsbc_planes.bcenter12[k][2];
         d = 0.5/sqrt(d);
         for(i=0;i<3;i++){
         rsbc_planes.bcenter12[k][i] *= d*(wat216.inner + wat216.outer);
         }
         }
         */
    }
    return;
}
