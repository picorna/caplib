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

// transform XYZ coordinates from the 2PLV system to 4RHV system, or inversely transform
// if from2plvto4rhv == 1, then rotate 2PLV to 4RHV
//                   != 1, then rotate 4RHV to 2PLV
void trans_2plvto4rhv(int from2plvto4rhv, int natoms, double *pdb_xyz) {
	
	double mat_p2r[3][3] = {0., 1., 0., -1., 0., 0., 0., 0., 1.};
	double mat_r2p[3][3] = {0., -1., 0., 1., 0., 0., 0., 0., 1.};
	double *xyz[3];
	double xyz_rot[3];
	int i;
	
	xyz[0] = pdb_xyz;
	xyz[1] = pdb_xyz+1;
	xyz[2] = pdb_xyz+1;
	
	if (from2plvto4rhv == 1) {
		for (i=0; i< natoms; i++) {
			xyz_rot[0] = mat_p2r[0][0]*(*xyz[0]) + mat_p2r[0][1]*(*xyz[1]) + mat_p2r[0][2]*(*xyz[2]);
			xyz_rot[1] = mat_p2r[1][0]*(*xyz[0]) + mat_p2r[1][1]*(*xyz[1]) + mat_p2r[1][2]*(*xyz[2]);
			xyz_rot[2] = mat_p2r[2][0]*(*xyz[0]) + mat_p2r[2][1]*(*xyz[1]) + mat_p2r[2][2]*(*xyz[2]);
			*xyz[0] = xyz_rot[0];
			*xyz[1] = xyz_rot[1];
			*xyz[2] = xyz_rot[2];
			pdb_xyz += 3;
		}
	} else {
		for (i=0; i< natoms; i++) {
			xyz_rot[0] = mat_r2p[0][0]*(*xyz[0]) + mat_r2p[0][1]*(*xyz[1]) + mat_r2p[0][2]*(*xyz[2]);
			xyz_rot[1] = mat_r2p[1][0]*(*xyz[0]) + mat_r2p[1][1]*(*xyz[1]) + mat_r2p[1][2]*(*xyz[2]);
			xyz_rot[2] = mat_r2p[2][0]*(*xyz[0]) + mat_r2p[2][1]*(*xyz[1]) + mat_r2p[2][2]*(*xyz[2]);
			*xyz[0] = xyz_rot[0];
			*xyz[1] = xyz_rot[1];
			*xyz[2] = xyz_rot[2];
			pdb_xyz += 3;
		}
	}
}