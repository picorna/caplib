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
#ifndef capgen_planes_h
#define capgen_planes_h

struct _rsbc_planes {
	double axis2[15][3]; //vectn
	//number of axis2 is half, because the inverse vector is considered to be the same for axis2.
	double axis3[20][3]; 
	double axis5[12][3];
	
	// rotation matrices used for 2PLV and 4RHV
	// rotmat_2plv[ ][ ]: rotation matrix for 2PLV
	// rotmat_4rhv[ ][ ]] rotation matrix for 4RHV
	double rotmat_2plv[60][3][3];
	double rotmat_4rhv[60][3][3];
    // rotmat_2plv_12[12][][]: rotate the cell 0 (1/12 cell) to other (1/12 cells)
    // invmat_2plv_12[12][][]: the inverse matrix of rotmat_2plv_12
    double rotmat_2plv_12[12][3][3];
    double invmat_2plv_12[12][3][3];
    // invmat_2plv and invmat_4rhv are the inverse matrices of rotmat_2plv and rotmat_4rhv, respectively.
    // Although invmat_4rhv is not used anywhere in practice, the matrix is prepared for perfection of the program.
    // invmat_2plv and invmat_4rhv are calculated in set_2plv_rotmat.c
    double invmat_2plv[60][3][3];
    double invmat_4rhv[60][3][3];
    
    // border planes
    double bplane[4][3],bplane12[5][3];
    //double bcenter[4][3],bcenter12[5][3];
    // the border atom flag;
    unsigned short *front;
};

struct _rsbc_planes rsbc_planes;

#endif
