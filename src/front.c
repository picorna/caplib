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
//
//  front.c
//  capgen
//
//  Created by yoneda on 2014/08/04.
//
//

#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include "wat216.h"
#include "planes.h"

int front(double *xyz,int border[5],double near){
    int ip,ipmax;
    double distance;
    double cutoff_cell;
    int near_to_plane;
    const unsigned short ibflag[] = {1,2,4,8,16,32,64,128,256,512};
    /*
     Distance from a plane that is defined with its center and its orthogonal vector.
     In my simulation program, the distance from a restricted plane was calculated previously. That is, the plane was defined with
     the orthogonal vector, the center, and the two additional vectors that denote the two terminal lines of the plane.
     However, I adopt the distance as a simple definition between an atom and a infinetely wide plane herein, because the complicated
     previous definition is not related the calculational speed nor precision.
     */
    
    // cutoff_cell < 0   when the atom is located at the back region of the plane
    
    if(wat216.type_of_rsbc == 60){
        ipmax = 4;
    } else {
        ipmax = 5;
    }
    near_to_plane = 0;
    for(ip=0;ip<ipmax;ip++){
        distance = (*xyz)*(rsbc_planes.bplane[ip][0]) + *(xyz+1)*(rsbc_planes.bplane[ip][1]) + *(xyz+2)*(rsbc_planes.bplane[ip][2]);
        // In the previous program, apricot, bcenter[][] was subtracted. However, I try to omit the subtraction.
        //        dx = (xyz[0] - rsbc_planes.bcenter[ip][0])*rsbc_planes.bplane[ip][0];
        //        dy = (xyz[1] - rsbc_planes.bcenter[ip][1])*rsbc_planes.bplane[ip][0];
        //        dz = (xyz[2] - rsbc_planes.bcenter[ip][2])*rsbc_planes.bplane[ip][0];
        //        d2 = dx*dx + dy*dy + dz*dz;
        if(distance < cutoff_cell) {
            near_to_plane += ibflag[ip];
        }
    }
    return near_to_plane;
}