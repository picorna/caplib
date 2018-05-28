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
#ifndef capgen_wat216_h
#define capgen_wat216_h

#define MAX_WAT_MOLS 216

struct _wat216 {
    int nmols; // usually 216
    // int natoms_in_a_mol=4; // 4 in the wat216.dat, but only 3 is used in capaxis
    // The 1st atom of a water molecule in wat216.dat is assumed to be oxygen.
    double cubic_size; // size of the water cubic box
    // A water molecule in wat216.dat contains 4 atoms.
    // However, only the coordinate of 3 atoms are read in capaxis, because the coordinates are used only for printing water atom coordinates in a pdb file.
    //double xyz[4][MAX_WAT_MOLS][3];};
    double xyz[3][MAX_WAT_MOLS][3];
    // The inner and outer radii of water region in rsbc
    double inner,outer;
    double short_contact,border_contact;
    // The type of rotational symmetry boundary condition
    // type_of_rsbc == 60  -> the ordinary 1/60 cell is used.
    // type of_rsbc == 12  -> the new 1/12 (pentamer) cell is used.
    int type_of_rsbc;
};

struct _wat216 wat216;

#endif
