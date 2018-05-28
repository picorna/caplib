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
#include <unistd.h>

#ifndef MAX_PDB_NAME_LENGTH
#define MAX_PDB_NAME_LENGTH 81
#endif

void usage(int status);

int call_getopt(int argc,char *argv[],int *argc_return,
    int *f_option_v,int *f_option_l,
    char fname_pdblist[MAX_PDB_NAME_LENGTH], char pdbname[MAX_PDB_NAME_LENGTH]) {
    
    char ch;
    extern char *optarg;
    extern int optind,opterr;

//    int f_option_l;        // A pdb file list is read from a separate file.
    *f_option_l = 0;
    *f_option_v = 0;
    
    // read parameters from argc and argv
    while ((ch = getopt(argc,argv,"vhl:")) != -1) {
        switch (ch){
            case 'l':
                *f_option_l = 1;
                strncpy(fname_pdblist, optarg, MAX_PDB_NAME_LENGTH);
                break;
            case 'h':
                usage(EXIT_SUCCESS);
                exit(EXIT_SUCCESS);                
            case 'v':
                *f_option_v = 1;
        }
    }

    *argc_return = argc;
    *argc_return -= optind;
    argv += optind;
    if(*f_option_l == 1) {
       if (*argc_return < 0) {
          // When no PDB filename list
          usage(EXIT_FAILURE);
          exit(EXIT_FAILURE);
       } else if (*argc_return > 0) {
          // More than one PDB filename lists
          strncpy(pdbname,argv[0],MAX_PDB_NAME_LENGTH);
          pdbname[MAX_PDB_NAME_LENGTH-1] = '\0';
       }
    } else {
       if (*argc_return == 0) {
          // When no filenames are written in the command line
          usage(EXIT_FAILURE);
          exit(EXIT_FAILURE);
       } else if (*argc_return >= 1) {
          // The PDB file name is read from the command line.
      }
    }
    return argc;
}
