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
#ifndef capgen_iequiv_h
#define capgen_iequiv_h

struct _contact {
    const char iequiv[9][9] = {
        //data ((iequiv(i,j),j=0,8),i=0,8)/
        {0, 1, 2, 3, 4, 5, 6, 7, 8},
        {8, 0, 7,-1,-1,-1,-1, 1, 2},
        {7, 2, 0, 5, 6,-1,-1, 8, 1},
        {3,-1, 4, 0, 2,-1,-1,-1,-1},
        {5,-1, 6, 7, 0, 4,-1,-1,-1},
        {4,-1,-1,-1, 5, 0, 2, 3,-1},
        {6,-1,-1,-1,-1, 7, 0, 5,-1},
        {2, 8, 1,-1,-1, 3, 4, 0, 7},
        {1, 7, 8,-1,-1,-1,-1, 2, 0}
    }
    const unsigned char iclnea[8]={4,3,5,2,1,6,0,4}
    const unsigned char icelnk[8]={4,0,5,1,2,6,3,4}
    }
        
        data mawari/1,3,4,11,2/

    const char iequiv12[12][12] = {
        //data ((iequi5(i,j),j=0,11),i=0,11)/
        {0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11},
        {2,  0,  1, 11,  9, 10,  5,  3,  4,  8,  6,  7},
        {1,  2,  0,  7,  8,  6, 10, 11,  9,  4,  5,  3},
        {3,  4,  5,  0,  1,  2,  9, 10, 11,  6,  7,  8},
        {11,  9, 10,  2,  0,  1,  8,  6,  7,  5,  3,  4},
        {7,  8,  6,  1,  2,  0,  4,  5,  3, 10, 11,  9},
        {6,  7,  8,  9, 10, 11,  0,  1,  2,  3,  4,  5},
        {5,  3,  4,  8,  6,  7,  2,  0,  1, 11,  9, 10},
        {10, 11,  9,  4,  5,  3,  1,  2,  0,  7,  8,  6},
        {9, 10, 11,  6,  7,  8,  3,  4,  5,  0,  1,  2},
        {8,  6,  7,  5,  3,  4, 11,  9, 10,  2,  0,  1},
        {4,  5,  3, 10, 11,  9,  7,  8,  6,  1,  2,  0}
    }
        koko
        data iclnea5/4,3,5,2,1,6,0,4/
    data icelnk5/4,0,5,1,2,6,3,4/
        // Through which plane, the cells contact each other.
    const unsigned char mawari/1,3,4,11,2/
}
