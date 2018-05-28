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
#ifndef MAX_PDB_NAME_LENGTH
#define MAX_PDB_NAME_LENGTH 81
#endif
#ifndef MAX_PDB_FILES
#define MAX_PDB_FILES 50000
#endif
#ifndef MAX_ERRORS
#define MAX_ERRORS 3
#endif

struct _store_ {
	char pdbname[MAX_PDB_NAME_LENGTH];
	double theta[60];
	double axis[60][3];
	int nerror;
	char error[MAX_ERRORS][81];
	char id[4];
	int pdb_type;
	int top[3];
	float top_percent[3];
};

struct _store_ store[MAX_PDB_FILES];
