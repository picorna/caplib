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
            char pdbname1[MAX_PDB_NAME_LENGTH],char pdbname2[MAX_PDB_NAME_LENGTH]){

    char ch;
    extern char *optarg;
    extern int  optind, opterr;
    
    // read parameters from argc and argv
    while ((ch = getopt(argc,argv,"h")) != -1) {
        switch (ch){
            case 'h': // help
                usage(0);
                exit(EXIT_SUCCESS);                
                break;
        }
    }

    argc -= optind;
    argv += optind;
    if (argc == 0) {
        // When no filenames are written in the command line
        usage(EXIT_FAILURE);
        exit(EXIT_FAILURE);
    } else if (argc == 2) {
        // The two PDB file names are read from the command line.
        strncpy(pdbname1,argv[0],MAX_PDB_NAME_LENGTH);
        pdbname1[MAX_PDB_NAME_LENGTH-1] = '\0';
        argc -= optind;
        argv += optind;
        strncpy(pdbname2,argv[0],MAX_PDB_NAME_LENGTH);
        pdbname2[MAX_PDB_NAME_LENGTH-1] = '\0';
    } else {
        // One or more than 3 file names are written in the command line
        usage(EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }    
    return argc;
}
