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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "summary_store.h"
// vectors, matrices, and planes or rotational symmetry
#include "ptable.h"
#include "planes.h"
#include "wat216.h"
//#include "functions.h"

#ifndef MAX_PDB_NAME_LENGTH
#define MAX_PDB_NAME_LENGTH 81
#endif

void read_a_pdb(FILE *,char [4],int *,int *,char *,int *,int *,double *,double [60][3][3],double [60][3],int *);
void gen_entire(int,int,int,double [3],int,char *,double *,double [60][3][3],double [60][3],double [60][3],double [60]);
struct _rsbc set_ptable();
void set_planes();
void set_2plv_rotmat(double [60][3][3],double [60][3][3],double [12][3][3],double [12][3][3],double [60][3][3],double [60][3][3]);
void read_biomt_file(FILE *,double [60][3][3],double [60][3],int *);
int analyse_biomt(int,double [60][3][3],double [60],double [60][3],int);
int summary_store(int,double [60],double [60][3],char *, char [4],int [3],float [3],int);
void summary_print(int);
void error_print(int);
void print_a_pdb(int,int,char *,double *,double [60][3][3],double [60][3],double [60][3],double [60]);
void print_center(int,int,char *,double *,char [4],double [60][3]);
void print_biomt(int,char *,double [60][3][3],double [60][3],double [60][3],double [60]);
int call_getopt(int,char *[],char [MAX_PDB_NAME_LENGTH],char [MAX_PDB_NAME_LENGTH]);
int cell_axis_fit(int,int,double *,double [60][3],double [60],double rot_fit[3][3]);
int cell_axis_fit_any(int,double *,double [60][3],double [60],double [3][3]);
void calc_cell(int,int,double *,double [60][3][3],int [3],float [3]);
void rotate_atoms(int,double *,double [3][3]);
int cell_number(double *);
void usage(int status);
int read_wat216(FILE *);
void set_bplane();

int main (int argc,char *argv[]) {	
	int f_verbose;	// output loudly (-v, f_verbose  == 1) 
					// Axes and angles of rotation matrices are printed at the BIOMT lines
	int f_input_stdin = 0;  // when no files are asigned in the command line, PDB data are read from stdin (f_input_stdin == 1)
	int f_option_m;		// The BIOMT lines are read from a separate file.
	int f_option_t;		// f_option_t == 1 --> 2PLV matrices are used.
						// f_option_t == 2 --> 4RHV matrices are used.
						// f_option_t == 0 --> the -t option is not used.
	//int f_option_rot_type; // When rotated coordinates are printed, the rotation matrix is selected 

	FILE *fin,*fin_matrix,*fin_wat216;
	
	char id[4],id_r[4];
	int finish,ierr;
	int ch,i,j,k;
	
	int nlines,nfiles,nlines_r;
	char pdbname1[MAX_PDB_NAME_LENGTH],pdbname2[MAX_PDB_NAME_LENGTH];
	char *pdb_line,*pdb_line_r;
	int natoms,natoms_r;
	double *pdb_xyz,*xtmp,*pdb_xyz_r;
	double theta[60],theta_t[60],axis[60][3],axis_t[60][3];
	double rot[60][3][3],rot_t[60][3][3];
	double xyzcent[60][3],xyzcent_t[60][3];
#ifndef PDB_INIT_NLINES
#define PDB_INIT_NLINES 500000
#endif
	int max_lines = PDB_INIT_NLINES;
	int max_atoms = PDB_INIT_NLINES;
	int nbiomt,nbiomt_t;
	int record_finish;
	int opt;
	int top[3];
	float top_percent[3];
	double rot_fit[3][3],bmatrix[3][3],amatrix[3][3];
	int there_is_xyzcent;
	double xyz_cent_entire[3];
	
	// set the external structures, rsbc_ptable and rsbc_planes
	rsbc = set_ptable();
	// set the partition planes
	set_planes();	
	// set the rotation matrices of 2PLV and 4RHV systems
	set_2plv_rotmat(rsbc_planes.rotmat_2plv,rsbc_planes.rotmat_4rhv,rsbc_planes.rotmat_2plv_12,rsbc_planes.invmat_2plv_12,rsbc_planes.invmat_2plv,rsbc_planes.invmat_4rhv);

	// read a run type and options from the command line
	opt = call_getopt(argc,argv,pdbname1,pdbname2);

	// allocate initial memory to store lines in a PDB file
	if((pdb_line  = (char *)malloc(sizeof(char)*81*PDB_INIT_NLINES)) == 0) {
		fprintf(stderr,"Error: memory exhausted for reading lines from a PDB file\n");
		exit(EXIT_FAILURE);
	}
	if((pdb_xyz = (double *)malloc(sizeof(double)*3*PDB_INIT_NLINES)) == 0) {
		fprintf(stderr,"Error: memory exhausted for reading xyz from a PDB file\n");
		exit(EXIT_FAILURE);
	}
	
	// allocate initial memory to store lines in a PDB file
	if((pdb_line_r  = (char *)malloc(sizeof(char)*81*PDB_INIT_NLINES)) == 0) {
		fprintf(stderr,"Error: memory exhausted for reading lines from -r PDB file\n");
		exit(EXIT_FAILURE);
	}
	if((pdb_xyz_r = (double *)malloc(sizeof(double)*3*PDB_INIT_NLINES)) == 0) {
		fprintf(stderr,"Error: memory exhausted for reading xyz from -r PDB file\n");
		exit(EXIT_FAILURE);
	}
	if ((fin_matrix = fopen(pdbname2, "r")) == NULL) {
		fprintf(stderr,"File open error for reading a rotation matrices in PDB with the -r option, %s\n",pdbname2);	
		exit(EXIT_FAILURE);
	}

	// read reference PDB file
	read_a_pdb(fin_matrix,id_r,&max_lines,&nlines_r,pdb_line_r,&max_atoms,&natoms_r,pdb_xyz_r,rot_t,xyzcent_t, &nbiomt_t);
	if (nbiomt_t != 180) {
		fprintf(stderr,"Error in reading -r file\n");
		exit(EXIT_FAILURE);
	}
	// determine axes and angles with sorting
	ierr = analyse_biomt(f_verbose,rot_t,theta_t,axis_t,0);
	if(ierr != 0) {
		fprintf(stderr,"Error in reading -r file\n");
		exit(EXIT_FAILURE);
	}
	// Because some PDB file has the "center of rotation" term in the BIOMT lines, once the 
	// entire structure is built to find the center or rotation. After finding the center,
	// the 5fold and 3fold axes closest to the protomer 0 are found.
	there_is_xyzcent = 0;
	for (i=0; i<60; i++) {
		if (fabs(xyzcent_t[i][0]) > 0.001 || fabs(xyzcent_t[i][1]) > 0.001 || fabs(xyzcent_t[i][2]) > 0.001) {
			there_is_xyzcent = 1;
			break;
		}
	}
	if (there_is_xyzcent == 1) {
		gen_entire(0,1,1,xyz_cent_entire,nlines_r,pdb_line_r,pdb_xyz_r,rot_t,xyzcent_t,axis_t,theta_t);
		xtmp = pdb_xyz_r;
		for (i=0; i<natoms; i++) {
			*(xtmp++) -= xyz_cent_entire[0];
			*(xtmp++) -= xyz_cent_entire[1];
			*(xtmp++) -= xyz_cent_entire[2];
		}
	}
	
        // read the PDB file that you rotate onto the reference PDB file
	if ((fin = fopen(pdbname1, "r")) == NULL) {
		fprintf(stderr,"file open error for PDB file, %s\n",pdbname1);
		free(pdb_line);
		free(pdb_xyz);
		exit(EXIT_FAILURE);
	}
	read_a_pdb(fin,id,&max_lines,&nlines,pdb_line,&max_atoms,&natoms,pdb_xyz,rot,xyzcent,&nbiomt);
	fclose(fin);

        // rotate a PDB file to make the rotation axes to match the axes of the other PDB file 
	// determine axes and angles with no sorting
	// print a rotated structure from a PDB file 
	// function, gen_entire, do not apply rotation if the second parameter is 0.
	// Thus, gen_entire is usde only for printing out.
	ierr = analyse_biomt(f_verbose,rot,theta,axis,0);
	if(ierr != 0) {
		fprintf(stderr,"Error in BIOMT lines\n");
		exit(EXIT_FAILURE);
	}
	// Because some PDB file has the "center of rotation" term in the BIOMT lines, once the 
	// entire structure is built to find the center or rotation. After finding the center,
	// the 5fold and 3fold axes closest to the protomer 0 are found.
	there_is_xyzcent = 0;
	for (i=0; i<60; i++) {
		if (fabs(xyzcent[i][0]) > 0.001 || fabs(xyzcent[i][1]) > 0.001 || fabs(xyzcent[i][2]) > 0.001) {
			there_is_xyzcent = 1;
			break;
		}
	}
	if (there_is_xyzcent == 1) {
		gen_entire(0,1,1,xyz_cent_entire,nlines,pdb_line,pdb_xyz,rot,xyzcent,axis,theta);
		xtmp = pdb_xyz;
		//printf("xyz_cent_entire = %8.3f %8.3f %8.3f\n",xyz_cent_entire[0],xyz_cent_entire[1],xyz_cent_entire[2]);
		for (i=0; i<natoms; i++) {
			*(xtmp++) -= xyz_cent_entire[0];
			*(xtmp++) -= xyz_cent_entire[1];
			*(xtmp++) -= xyz_cent_entire[2];
		}
	}
	ierr = cell_axis_fit_any(natoms,  pdb_xyz,  axis,  theta,  bmatrix);
	ierr = cell_axis_fit_any(natoms_r,pdb_xyz_r,axis_t,theta_t,amatrix);
	//fprintf(stderr,"%d %d\n",natoms,natoms_r);
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			rot_fit[i][j] = amatrix[i][0]*bmatrix[j][0] + amatrix[i][1]*bmatrix[j][1] + amatrix[i][2]*bmatrix[j][2];
		}
		//fprintf(stderr, "%8.3f %8.3f %8.3f\n",rot_fit[i][0],rot_fit[i][1],rot_fit[i][2]);
	}
	rotate_atoms(natoms,pdb_xyz,rot_fit);
	print_a_pdb(f_verbose,nlines,pdb_line,pdb_xyz,rot,xyzcent,axis,theta);

	free(pdb_line);
	free(pdb_xyz);
	free(pdb_line_r);
	free(pdb_xyz_r);
	
	exit (EXIT_SUCCESS);
}
