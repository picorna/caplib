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
 *  read_a_pdb - read a Protein Databank file. This function is a part of the program capaxis
 *
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef UNITMATRIX_LIMIT
#define UNITMATRIX_LIMIT 0.001
#endif
#ifndef PDB_ADD_NLINES
#define PDB_ADD_NLINES 500000
#endif

void read_a_pdb(FILE *fin,char id[4],int *max_lines,int *nlines,char *pdb_line,
   int *max_atoms,int *natoms,double *pdb_xyz,double rot[60][3][3],double xyzcent[60][3],int *nbiomt) {	
	
	int record_finish,finish;
	int i,nrot;
	int ch;
	double *xyz;
	char *line;
			
	/* read lines in the PDB file */
	id[0] = ' ';
	id[1] = ' ';
	id[2] = ' ';
	id[3] = ' ';
	*nlines = 0;
	*natoms = 0;
	*nbiomt = 0;
	finish  = 0;
	while(finish == 0) {
		line = pdb_line + 81*(*nlines);
		// read a line in the PDB file 
		/* fgets(line, 81, fin) reads 80 characters and adds \0 at the 81th position
		   When \n appears from 0 to 79, \0 is added after \n */
		if(fgets(line, 81, fin) == NULL) {
			return;
		} else {
			if ((*nlines)+1 > *max_lines) {
				fprintf(stderr,"Memory reallocation to store large PDB files\n");
				*max_lines += PDB_ADD_NLINES;
				if((pdb_line = (char *)realloc(pdb_line, (*max_lines)*81*sizeof(char))) == 0) {
					fprintf(stderr,"Error: memory exhausted for reading a PDB file\n");
					free(pdb_line);
					free(pdb_xyz);					
					exit(EXIT_FAILURE);
				}
			}
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
				while ((ch = fgetc(fin)) != '\n') {
					if (ch == EOF) {
						finish = 1;
						break;
					}
				}
			}
		}
		(*nlines)++;


		/* for an ATOM or HETATM line */
		if (strncmp(line,"ATOM",4) == 0 || strncmp(line,"HETATM",6) == 0){
			xyz = pdb_xyz+3*(*natoms);
			sscanf(line+30, "%8lf%8lf%8lf",xyz,xyz+1,xyz+2);
			if ((*natoms)+1 > *max_atoms) {
				*max_atoms += PDB_ADD_NLINES;
				if((pdb_xyz = realloc(pdb_xyz, (*max_atoms)*3*sizeof(double))) == 0) {
					fprintf(stderr,"Error: memory exhausted for reading XYZ in a PDB file\n");
					free(pdb_line);
					free(pdb_xyz);
					exit(EXIT_FAILURE);
				}
			}
			(*natoms)++;
		/* for a BIOMT line */
		} else if ((strncmp(line,"REMARK",6) == 0) ){
			if((strncmp(line+10,"   BIOMT",8) == 0)) {
				/* eading the rotation matrix */
				if (*nbiomt<180) {
				nrot = (*nbiomt) / 3;
				i    = (*nbiomt) % 3;
				//01234567890123456789012012345678901234567890123456789012345678901234
				//REMARK 350   BIOMT2  11  0.809017  0.309017 -0.500000     -208.52413
				sscanf(line+24, "%10lf%10lf%10lf%15lf",&rot[nrot][i][0],&rot[nrot][i][1],
						   &rot[nrot][i][2],&xyzcent[nrot][i]);
				}
				(*nbiomt)++;
			}
		} else if (strncmp(line,"HEADER",6)==0) {
			//123456789012345678901234567890123456789012345678901234567890123
			//HEADER    VIRUS                                   25-JAN-88   4RHV 
			id[0] = line[62];
			id[1] = line[63];
			id[2] = line[64];
			id[3] = line[65];
		}
	}
}
