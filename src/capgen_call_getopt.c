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
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifndef MAX_PDB_NAME_LENGTH
#define MAX_PDB_NAME_LENGTH 81
#endif
#ifndef MAX_CELLLIST_LENGTH
#define MAX_CELLLIST_LENGTH 131
#endif

void usage(int status);

int call_getopt(int argc,char *argv[],
            int *f_verbose,int *f_option_m,int *f_option_t,
            char fname_r[MAX_PDB_NAME_LENGTH],int cells[60],
            char fname_m[MAX_PDB_NAME_LENGTH],char pdbname[MAX_PDB_NAME_LENGTH]) {

    char ch;
    extern char *optarg;
    extern int  optind, opterr;
    char temp[MAX_CELLLIST_LENGTH], *cell_name, *cell_hyphen;

    *f_verbose  = 0;
    *f_option_m = 0;
    *f_option_t = 0;
    int k,k1,k2;
    int if_cell_list = 0;
    
    // read parameters from argc and argv
    while ((ch = getopt(argc,argv,"vht:m:r:c:")) != -1) {
        switch (ch){
            case 'c': // list of cells
                if_cell_list = 1; 
                for(k=0;k<60;k++) {cells[k]=0;}
                strncpy(temp, optarg, MAX_CELLLIST_LENGTH);
                if((cell_name = strtok(temp, ",")) != NULL) {
                    if(cell_name == '\0') {
                        exit(EXIT_FAILURE);
                    }
                    if((cell_hyphen = strstr(cell_name,"-"))){
                        k2 = atoi(cell_hyphen+1);
                        *cell_hyphen = '\0';
                        k1 = atoi(cell_name);
                        for(k=k1;k<=k2;k++){ cells[k] = 1; }    
                    } else {
                        k = atoi(cell_name);
                        cells[k] = 1;
                    }
                } else {
                        k = atoi(cell_name);
                        break;
                }
                while((cell_name!=NULL)) {
                    if((cell_name = strtok(NULL, ",")) != NULL) {
                        if((cell_hyphen = strstr(cell_name,"-"))){
                            k2 = atoi(cell_hyphen+1);
                            *cell_hyphen = '\0';
                            k1 = atoi(cell_name);
                            for(k=k1;k<=k2;k++){ cells[k] = 1; }    
                        } else {
                            k = atoi(cell_name);
                            cells[k] = 1;
                        }
                    }
                }
                break;
            case 'r': // print a rotated pdb file, but this is not used
                usage(EXIT_FAILURE);
                exit(EXIT_FAILURE);
            case 'm': // use the rotation matrices in another PDB file
                if(*f_option_m == 1) {
                    usage(EXIT_FAILURE);
                    exit(EXIT_FAILURE);
                }
                *f_option_m = 1;
                strncpy(fname_m, optarg, MAX_PDB_NAME_LENGTH);
                break;
            case 'h': // help
                usage(0);
                exit(EXIT_SUCCESS);                
            case 'v': // vervose
                if (*f_verbose != 0) {
                    usage(EXIT_FAILURE);
                    exit(EXIT_FAILURE);
                }
                *f_verbose = 1;
                break;
            case 't': // No BIOMT lines are used. In stead, you can select 2PLV or 4RHV matrix.
                printf("%c \n",optarg[0]);
                if (optarg[0]=='p') {
                    // use the 2PLV matrices
                    *f_option_t = 1;
                } else if(optarg[0]=='r') {
                    // use the 2PLV matrices, assuming the orientation of the 4RHV coordinate system
                    *f_option_t = 2;
                } else {
                    fprintf(stderr, "The -t option needs the p or r parameters either\n");
                    usage(EXIT_FAILURE);
                    exit(EXIT_FAILURE);
                }
                break;
        }
    }    
    if(if_cell_list == 0) {
        for(k=0;k<60;k++) {cells[k]=1;}
    }

    argc -= optind;
    argv += optind;
    //check for consistency of options
    if (argc == 0) {
        // When no filenames are written in the command line
        usage(EXIT_FAILURE);
        exit(EXIT_FAILURE);
    } else if (argc == 1) {
        // The PDB file name is read from the command line.
        strncpy(pdbname,argv[0],MAX_PDB_NAME_LENGTH);
        pdbname[MAX_PDB_NAME_LENGTH-1] = '\0';
    } else {
        // Two or more file names are written in the command line
        usage(EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }    
    return argc;
}
