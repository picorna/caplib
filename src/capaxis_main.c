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
int call_getopt(int,char *[],int *,int *,int *,char [MAX_PDB_NAME_LENGTH],char [MAX_PDB_NAME_LENGTH]);
int cell_axis_fit(int,int,double *,double [60][3],double [60],double rot_fit[3][3]);
int cell_axis_fit_any(int,double *,double [60][3],double [60],double [3][3]);
void calc_cell(int,int,double *,double [60][3][3],int [3],float [3]);
void rotate_atoms(int,double *,double [3][3]);
int cell_number(double *);
void usage(int status);
int read_wat216(FILE *);
void set_bplane();

int main (int argc,char *argv[]) {    
    int f_option_v = 0;
        // if 1, output loudly (-v)
        // Axes and angles of rotation matrices are printed at the BIOMT lines
    int f_option_l;        // A pdb file list is read from a separate file, if 1.
    
    FILE *fin,*fin_matrix,*fin_pdblist,*fin_wat216;
    
    char id[4],id_r[4];
    int finish,ierr;
    int ch,i,j,k;
    
    int nlines,nfiles,nlines_r;
#ifndef MAX_PDB_NAME_LENGTH
#define MAX_PDB_NAME_LENGTH 81
#endif
    char pdbname[MAX_PDB_NAME_LENGTH];
    char fname_pdblist[MAX_PDB_NAME_LENGTH];
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
    int argc_return;
    
    // set the external structures, rsbc_ptable and rsbc_planes
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

    // read a run type and options from the command line
    opt = call_getopt(argc,argv,&argc_return,&f_option_v,&f_option_l,fname_pdblist,pdbname);

    // allocate initial memory to store lines in a PDB file
    if((pdb_line  = (char *)malloc(sizeof(char)*81*PDB_INIT_NLINES)) == 0) {
        fprintf(stderr,"Error: memory exhausted for reading lines from a PDB file\n");
        exit(EXIT_FAILURE);
    }
    if((pdb_xyz = (double *)malloc(sizeof(double)*3*PDB_INIT_NLINES)) == 0) {
        fprintf(stderr,"Error: memory exhausted for reading xyz from a PDB file\n");
        exit(EXIT_FAILURE);
    }
    
    if (f_option_l == 1) {
        if ((fin_pdblist = fopen(fname_pdblist, "r")) == NULL) {
            fprintf(stderr,"file open error for reading a PDB file list with the -l option, %s\n",fname_pdblist);
            exit(EXIT_FAILURE);
        }
    }
    
    nfiles = 0;
    finish = 0;
    while (finish == 0) {
    if(f_option_l == 1) {
        // A list of PDB file names are read from a list file determined with the -l option.
        if(fgets(pdbname, 81, fin_pdblist) == NULL) {
            finish = 1;
            fclose(fin_pdblist);
            break;
        } else {
            record_finish = 0;
            for(i=0;i<MAX_PDB_NAME_LENGTH;i++) {
                if (pdbname[i] == '\n') {
                    pdbname[i] = '\0';
                    record_finish = 1;
                    break;
                }
            }
            if (record_finish == 0) {
                pdbname[MAX_PDB_NAME_LENGTH-1] = '\0';
                while ((ch = fgetc(fin_pdblist)) != '\n') {
                    if (ch == EOF) {
                        finish = 1;
                        break;
                    }
                }
            }
            if ((fin = fopen(pdbname, "r")) == NULL) {
                fprintf(stderr,"file open error for PDB file, %s\n",pdbname);
                free(pdb_line);
                free(pdb_xyz);
                fclose(fin_pdblist);
                exit(EXIT_FAILURE);
            }
            read_a_pdb(fin,id,&max_lines,&nlines,pdb_line,&max_atoms,&natoms,pdb_xyz,rot,xyzcent,
                       &nbiomt);
            fclose(fin);
        }
    } else {
        argc_return -= 1;
        if (argc_return<=0) {
            argv++;
            finish = 1;
        } else {
            argv++;
        }
        strncpy(pdbname,argv[0],MAX_PDB_NAME_LENGTH);
        pdbname[MAX_PDB_NAME_LENGTH-1] = '\0';
        if ((fin = fopen(pdbname, "r")) == NULL) {
            fprintf(stderr,"file open error for PDB file, %s\n",argv[0]);
            free(pdb_line);
            free(pdb_xyz);
            exit(EXIT_FAILURE);
        }
        read_a_pdb(fin,id,&max_lines,&nlines,pdb_line,&max_atoms,&natoms,pdb_xyz,rot,xyzcent,&nbiomt);
        fclose(fin);
    }

    // analyse the rotation axes from the rotation matrices in the BIOMT lines
    if (nbiomt == 180) {
        // determine axes and angles with sorting
        ierr = analyse_biomt(f_option_v,rot,theta,axis,1);
        if(ierr == 0) {
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
            cell_axis_fit(f_option_v,natoms,pdb_xyz,axis,theta,rot_fit);
            if (ierr == 0) {
                rotate_atoms(natoms,pdb_xyz,rot_fit);
                //printf("%c%c%c%c\n",id[0],id[1],id[2],id[3]);
                calc_cell(f_option_v,natoms,pdb_xyz,rot,top,top_percent);
            }
        }
        ierr = summary_store(nfiles,theta,axis,pdbname,id,top,top_percent,f_option_v);
        nfiles++;
    } else {
        fprintf(stderr, "No.%d %c%c%c%c has the number of BIOMT lines of not 180, but %d\n",
                nfiles,id[0],id[1],id[2],id[3],nbiomt);
    }
    }
    free(pdb_line);
    free(pdb_xyz);
    
    // print the final results of analysis
    fprintf(stdout,"\n< Summary of Analysis on BIOMT lines >\n");
    summary_print(nfiles);
    error_print(nfiles);
    exit (EXIT_SUCCESS);
}
