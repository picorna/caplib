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

#include "ptable.h"
#include "planes.h"
#include "wat216.h"

#ifndef MAX_PDB_NAME_LENGTH
#define MAX_PDB_NAME_LENGTH 81
#endif

void read_a_pdb(FILE *,char [4],int *,int *,char *,int *,
                int *,double *,double [60][3][3],double [60][3],int *);
void gen_entire(int [60],int,int,char *,double *,
                double [60][3][3],double [60][3]);
struct _rsbc set_ptable();
void set_planes();
void set_2plv_rotmat(double [60][3][3],double [60][3][3],
                double [12][3][3],double [12][3][3],
                double [60][3][3],double [60][3][3]);
void read_biomt_file(FILE *,double [60][3][3],double [60][3],int *);
int call_getopt(int,char *[],int *,int *,int *,
                char [MAX_PDB_NAME_LENGTH],int [60],
                char [MAX_PDB_NAME_LENGTH],char [MAX_PDB_NAME_LENGTH]);
int cell_number(double *);
void usage(int status);

int main(int argc,char *argv[]) {
#ifndef DEBUG
#define DEBUG 0
#endif
    int f_option_m;        // The BIOMT lines are read from a separate file.
    int f_verbose = 0;
    int f_option_t = 1;
                // f_option_t == 1 --> 2PLV matrices are used.
                // f_option_t == 2 --> 4RHV matrices are used.
                // f_option_t == 0 --> the -t option is not used.
    int cells[60];
    FILE *fin,*fin_matrix;
    char id[4],id_r[4];
    int ierr;
    int ch,i,j,k;
    
    int nlines,nlines_r;
#ifndef MAX_PDB_NAME_LENGTH
#define MAX_PDB_NAME_LENGTH 81
#endif
    char pdbname[MAX_PDB_NAME_LENGTH];
    char fname_m[MAX_PDB_NAME_LENGTH],fname_r[MAX_PDB_NAME_LENGTH];
    char *pdb_line;
    int natoms;
    double *pdb_xyz;
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
    double xyz_cent_entire[3];
    
    // set rsbc_ptable and rsbc_planes
    rsbc = set_ptable();
    // set the partition planes
    set_planes();    
    // set the rotation matrices of 2PLV and 4RHV systems
    set_2plv_rotmat(rsbc_planes.rotmat_2plv,
                    rsbc_planes.rotmat_4rhv,
                    rsbc_planes.rotmat_2plv_12,
                    rsbc_planes.invmat_2plv_12,
                    rsbc_planes.invmat_2plv,
                    rsbc_planes.invmat_4rhv);

    // read options from the command line
    opt = call_getopt(argc,argv,&f_verbose,&f_option_m,&f_option_t,
             fname_r,cells,fname_m,pdbname);

    if(DEBUG>=1) {for(k=0;k<60;k++){printf("cells =%3d\n",cells[k]);}} 
    
    if (f_option_m == 1) {
// The rotation matrices are read from the seperate file written in the PDB format
        if ((fin_matrix = fopen(fname_m, "r")) == NULL) {
            fprintf(stderr,
                "file open error for reading a rotation matrices with the -m option, %s\n",
                fname_m);    
            exit(EXIT_FAILURE);
        }
        read_biomt_file(fin_matrix,rot_t,xyzcent_t,&nbiomt_t);
    }
    
    // allocate initial memory to store lines in a PDB file
    if((pdb_line  = (char *)malloc(sizeof(char)*81*PDB_INIT_NLINES)) == 0) {
        fprintf(stderr,"Error: memory exhausted for reading lines from a PDB file\n");
        exit(EXIT_FAILURE);
    }
    if((pdb_xyz = (double *)malloc(sizeof(double)*3*PDB_INIT_NLINES)) == 0) {
        fprintf(stderr,"Error: memory exhausted for reading xyz from a PDB file\n");
        exit(EXIT_FAILURE);
    }

// read a PDB file
    if ((fin = fopen(pdbname, "r")) == NULL) {
        fprintf(stderr,"file open error for PDB file, %s\n",argv[0]);
        fprintf(stderr,"file name is %s\n",pdbname);
        free(pdb_line);
        free(pdb_xyz);
        exit(EXIT_FAILURE);
    }
    read_a_pdb(fin,id,&max_lines,&nlines,pdb_line,&max_atoms,&natoms,pdb_xyz,rot,xyzcent,&nbiomt);
    fclose(fin);
    free(pdb_line);
    free(pdb_xyz);

//Finally, the selected protein units with cells[60] are generated
    if(f_option_t == 1) {
        printf("2PLV rotation matrices used\n");
    } else if(f_option_t == 2) {
        printf("4RHV rotation matrices used\n");
    }
    if(f_option_t == 1 || f_option_t == 2) {
        for(i=0;i<60;i++) {
            if(cells[i]==1){
                for(j=0;j<3;j++) {
                for(k=0;k<3;k++) {
                    if(f_option_t == 1) {
                        rot_t[i][j][k] = rsbc_planes.rotmat_2plv[i][j][k];
                    } else if(f_option_t == 2) {
                        rot_t[i][j][k] = rsbc_planes.rotmat_4rhv[i][j][k];
                    }
                }
                }
            }
        }
        gen_entire(cells,f_verbose,nlines,pdb_line,pdb_xyz,rot_t,xyzcent);
    } else {
        if (f_option_m == 1) {
            if (nbiomt_t != 180) {
                fprintf(stderr, "The entire capsid structure can not be built, \n");
                fprintf(stderr, "because the number of BIOMT lines for -m option is %d, not equal to 180\n",nbiomt);
                fprintf(stderr, "Try the -t option, by which you can use the rotation matrices of 2PLV or 4RHV\n");
            } else {
                gen_entire(cells,f_verbose,nlines,pdb_line,pdb_xyz,rot_t,xyzcent_t);
            }
    	} else {
            if (nbiomt != 180) {
                fprintf(stderr, "The entire capsid structure can not be built, \n");
                fprintf(stderr, "because the number of BIOMT lines is %d, not equal to 180\n",nbiomt);
                fprintf(stderr, "Try the -t or -m option, by which you can use the rotation matrices of 2PLV or 4RHV\n");
            } else {
                gen_entire(cells,f_verbose,nlines,pdb_line,pdb_xyz,rot,xyzcent);
            }
        }
    }
    exit (EXIT_SUCCESS);
}
