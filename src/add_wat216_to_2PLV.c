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
#include <math.h>
#include "wat216.h"
#include "planes.h"

int cell_number(double *);
int cell_number12(double *);

//
// Attention: This function performs the "packing" on the atom base, whereas the packing is done on the residue base in the apricot program.

void add_wat216_to_2PLV(int natoms,double *pdb_xyz,char *pdb_line){
    double *x,*y,*z;
    double xyz_rotated[3];
    int nedges;
    double edge[5][3],edge_max[3],edge_min[3];
    double wmin[3];
    double d,d1,d2,d3;
    double water_oxygen[3],water_hydrogen[3],w216_min[3];
    int i,j,m,p;
    int i0,i1,i2;
    int cell,cell12,nneighbor_cells;
    double inner2,outer2,short_contact2;
    double *rmat;
    int dont_adopt;
    int nborder;
    unsigned short short_int;
    
    if (wat216.type_of_rsbc != 60 && wat216.type_of_rsbc != 12) {
        fprintf(stderr,"Bug in program: wat216.type_of_rscb must be 60 or 12 in add_wat216_to_2PLV\n");
        exit(EXIT_FAILURE);
    }

    // Assume the PDB coordinates are in the 2PLV coordinate system. If not, apply this program with the -r
    // option to move the PDB into the 2PLV system.
    //
    // 1. calculate the cell number for each atom
    // 2. rotate the atom into the cell 0 with the inverse matrix
    // Thus all the atomic coordinates of protein are packed within cell 0
    //
    x = pdb_xyz;
    for (i=0; i<natoms; i++) {
        y = x+1;
        z = x+2;
        if(wat216.type_of_rsbc == 60){
            cell = cell_number(x);
            xyz_rotated[0] = rsbc_planes.invmat_2plv[cell][0][0]*(*x) + rsbc_planes.invmat_2plv[cell][0][1]*(*y) + rsbc_planes.invmat_2plv[cell][0][2]*(*z);
            xyz_rotated[1] = rsbc_planes.invmat_2plv[cell][1][0]*(*x) + rsbc_planes.invmat_2plv[cell][1][1]*(*y) + rsbc_planes.invmat_2plv[cell][1][2]*(*z);
            xyz_rotated[2] = rsbc_planes.invmat_2plv[cell][2][0]*(*x) + rsbc_planes.invmat_2plv[cell][2][1]*(*y) + rsbc_planes.invmat_2plv[cell][2][2]*(*z);
        } else {
            cell = cell_number(x);
            cell12 = cell % 5;
            xyz_rotated[0] = rsbc_planes.invmat_2plv_12[cell][0][0]*(*x) + rsbc_planes.invmat_2plv_12[cell][0][1]*(*y) + rsbc_planes.invmat_2plv_12[cell][0][2]*(*z);
            xyz_rotated[1] = rsbc_planes.invmat_2plv_12[cell][1][0]*(*x) + rsbc_planes.invmat_2plv_12[cell][1][1]*(*y) + rsbc_planes.invmat_2plv_12[cell][1][2]*(*z);
            xyz_rotated[2] = rsbc_planes.invmat_2plv_12[cell][2][0]*(*x) + rsbc_planes.invmat_2plv_12[cell][2][1]*(*y) + rsbc_planes.invmat_2plv_12[cell][2][2]*(*z);
        }

        //    2. rotate each atom by the inverse matrix for the cell number
        *x = xyz_rotated[0];
        *y = xyz_rotated[1];
        *z = xyz_rotated[2];
        x += 3;
    }
    
    //
    // 3. Calculate the size of rectangular box of water that contains cell 0
    if(wat216.type_of_rsbc == 60) {
        // The four edges of cell 0 are
        //  rsbc_planes.axis5[0][]  (up),
        //  rsbc_planes.axis3[4][]  (up),
        //  rsbc_planes.axis2[6][]  (down),
        //  rsbc_planes.axis2[9][]  (up)
        // (These are axis 1, axis 1-2-3, Z axis (axis 1-2),and axis 1-3.)
        // They are named "edge[4][3]" in this function.
        nedges = 4;
        for(i=0;i<3;i++){
            edge[0][i] =  rsbc_planes.axis5[0][i];
            edge[1][i] =  rsbc_planes.axis3[4][i];
            edge[2][i] = -rsbc_planes.axis2[6][i];
            edge[3][i] =  rsbc_planes.axis2[9][i];
        }
    } else {
        // for 1/12 cells,
        // The five edges of cell 0 (1/12) are the 3fold axes, 1-2-3, 1-3-12, 1-5-12, 1-4-5, 1-2-4. They are as follows, respectively,
        //  rsbc_planes.axis3[4][]
        //  rsbc_planes.axis3[13][]
        //  rsbc_planes.axis3[0][]
        //  rsbc_planes.axis3[12][]
        //  rsbc_planes.axis3[5][]
        // They are named "edge[5][3]" in this function.
        nedges = 5;
        for(i=0;i<3;i++){
            edge[0][i] =  rsbc_planes.axis3[4][i];
            edge[1][i] =  rsbc_planes.axis3[13][i];
            edge[2][i] =  rsbc_planes.axis3[0][i];
            edge[3][i] =  rsbc_planes.axis3[12][i];
            edge[4][i] =  rsbc_planes.axis3[5][i];
        }
    }
    // normalize edge[][]
    // I know that the lengths of axes (For example, the lengths of rsbc_planes.axis5, rsbc_planes.axis3, and rsbc_planes.axis2 are (1/2)*sqrt((5+sqrt(5))/2), (3+sqrt(5))*sqrt(3)/12, and sqrt((3+sqrt(5))/8), respectively.
    // However, I do not use this knowledge. The normalization is done numerically in this function.
    for(j=0;j<nedges;j++){
        d = edge[j][0]*edge[j][0] + edge[j][1]*edge[j][1] + edge[j][2]*edge[j][2];
        d = 1.0/sqrt(d);
        edge[j][0] *= d;
        edge[j][1] *= d;
        edge[j][2] *= d;
    }
    
    // find the size of the rectangular box that contains the cell 0 (1/60 or 1/12) with the inner and outer water radii,
    // wat216.inner and wat216.outer.
    //
    // The water region is a part of spherical shell with inner and outer radii, wat216.inner and wat216.outer,
    // so that determination of the rectangular box size needs a complicated calculation about the outer curved
    // surface in principle.  However, the 2PLV coordinate system is adopted so that the 3fold axis of cell 0 is
    // equal to the Z axis (ordinary rsbc with 1/60 cell, or the center of cell is the Z axis (1/12 cell). Therefore, the box size is determined only from the intersections of the axes (edges) including the Z axis (center) and the inner and outer sphere surface.
    for(i=0;i<3;i++){
        edge_max[i] = -1.0e30;
        edge_min[i] =  1.0e30;
        for(j=0;j<nedges;j++){
            d = edge[j][i]*wat216.inner;
            if(d < edge_min[i]){
                edge_min[i] = d;
            } else if(d > edge_max[i]){
                edge_max[i] = d;
            }
            d = edge[j][i]*wat216.outer;
            if(d < edge_min[i]){
                edge_min[i] = d;
            } else if(d > edge_max[i]){
                edge_max[i] = d;
            }
            // for the Z axis
            if(i == 2){
                d = wat216.inner;
            } else {
                d = 0;
            }
            if(d < edge_min[i]){
                edge_min[i] = d;
            } else if(d > edge_max[i]){
                edge_max[i] = d;
            }
            if(i == 2){
                d = wat216.outer;
            } else {
                d = 0;
            }
            if(d < edge_min[i]){
                edge_min[i] = d;
            } else if(d > edge_max[i]){
                edge_max[i] = d;
            }
            // The above "if"s are superflourious and they can be more simplified if you see the figure of icosahedron and cells, as shown in the articles. However, I wrote the above for simplicity in the program where no speed-up is needed.
        }
    }
    
    // 4. calculate the left low corner of the cubic box of wat216.dat
    // The 1st atom of a water molecule in wat216.dat is assumed to be oxygen.
    // wat216.xyz[3 or 4 atoms][x,y,z][216 mols]
    for(i=0;i<3;i++){
        wmin[i] = 1.0e30;
        for(j=0;j<wat216.nmols;j++){
            d = wat216.xyz[0][i][j];
            if(d < wmin[i]){
                wmin[i] = d;
            }
        }
    }

    
    // 5. The number of wat216 cubics that are needed to fill the water region is,
    //  ((edge_max[i]-edge_min[i])/wat216.cubic_size)+1
    // 6. For each water molecules of wat216 cubics, whether the criteria of water region are checked.
    // The 1st atom of a water molecule in wat216.dat is assumed to be oxygen.
    inner2 = wat216.inner*wat216.inner;
    outer2 = wat216.outer*wat216.outer;
    short_contact2 = wat216.short_contact*wat216.short_contact;
    if(wat216.type_of_rsbc == 60){
        nneighbor_cells = 8;
    } else {
        nneighbor_cells = 5;
    }
    for (i0=0;i0<((edge_max[0]-edge_min[0])/wat216.cubic_size)+1;i0++){
    for (i1=0;i1<((edge_max[1]-edge_min[1])/wat216.cubic_size)+1;i1++){
    for (i2=0;i2<((edge_max[2]-edge_min[2])/wat216.cubic_size)+1;i2++){
        for (m=0;m<wat216.nmols;m++){
            // wat216.xyz[3 or 4 atoms][x,y,z][216 mols]
            water_oxygen[0] = wat216.xyz[0][0][m] + i0*wat216.cubic_size + edge_min[0] - w216_min[0];
            water_oxygen[1] = wat216.xyz[0][1][m] + i1*wat216.cubic_size + edge_min[1] - w216_min[1];
            water_oxygen[2] = wat216.xyz[0][2][m] + i2*wat216.cubic_size + edge_min[2] - w216_min[2];
            // Oxygen of water must be in cell 0.
            if(wat216.type_of_rsbc == 60) {
                cell = cell_number(water_oxygen);
                if(cell != 0){
                    continue;
                }
            } else {
                cell = cell_number(water_oxygen);
                cell12 = cell % 5;
                if(cell12 != 0){
                    continue;
                }
            }

            // The distance of oxygen from origin must be from wat216.inner to wat216.outer.
            d = water_oxygen[0]*water_oxygen[0] + water_oxygen[1]*water_oxygen[1] + water_oxygen[2]*water_oxygen[2];
            if(d < inner2 || d > outer2){
                continue;
            }
            
            // Short contact with proteins in cell 0 and its neighbor cells
            dont_adopt = 0;
            for(p=0;p<nneighbor_cells;p++){
                if(p != 0){
                    if(wat216.type_of_rsbc == 60){
                        rmat = &rsbc_planes.rotmat_2plv[p][0][0];
                        nborder = 4;
                    } else {
                        rmat = &rsbc_planes.rotmat_2plv_12[p][0][0];
                        nborder = 5;
                    }
                    x = pdb_xyz;
                    for(i=0;i<natoms;i++){
                        // The border atom flag (In contrast, the APRICOT program uses the border "residue" flags).
                        // if x is close to the plane p, front() is 1
                        short_int = *(rsbc_planes.front+i);
                        short_int = short_int >> p;
                        if((short_int & 0x0001) == 0){
                            x += 3;
                            continue;
                        }
                        // rotate each protein atom by the matrix for cell p
                        y = x+1;
                        z = x+2;
                        xyz_rotated[0] = (*rmat)*(*x)     + (*(rmat+1))*(*y) + (*(rmat+2))*(*z);
                        xyz_rotated[1] = (*(rmat+3))*(*x) + (*(rmat+4))*(*y) + (*(rmat+5))*(*z);
                        xyz_rotated[2] = (*(rmat+6))*(*x) + (*(rmat+7))*(*y) + (*(rmat+8))*(*z);
                        x += 3;
                        d1 = water_oxygen[0] - xyz_rotated[0];
                        d2 = water_oxygen[1] - xyz_rotated[1];
                        d3 = water_oxygen[2] - xyz_rotated[2];
                        d = d1*d1 + d2*d2 + d3*d3;
                        if(d < short_contact2){
                            dont_adopt = 1;
                            break;
                        }
                    }
                } else {
                    x = pdb_xyz;
                    for(i=0;i<natoms;i++){
                        y = x+1;
                        z = x+2;
                        d1 = water_oxygen[0] - *x;
                        d2 = water_oxygen[1] - *y;
                        d3 = water_oxygen[2] - *z;
                        x += 3;
                        d = d1*d1 + d2*d2 + d3*d3;
                        if(d < short_contact2){
                            dont_adopt = 1;
                            break;
                        }
                    }
                }
                if(dont_adopt == 1){
                    break;
                }
            }
            
            // Finally, the atoms of the water m are added to the simulation
            if(dont_adopt == 0) {
                printf("ATOM   33 O  HOH 333 %8.3f%8.3f%8.3f\n",water_oxygen[0],water_oxygen[1],water_oxygen[2]);
                water_hydrogen[0] = wat216.xyz[1][0][m] + i0*wat216.cubic_size + edge_min[0] - w216_min[0];
                water_hydrogen[1] = wat216.xyz[1][1][m] + i1*wat216.cubic_size + edge_min[1] - w216_min[1];
                water_hydrogen[2] = wat216.xyz[1][2][m] + i2*wat216.cubic_size + edge_min[2] - w216_min[2];
                printf("ATOM   33 H1  HOH 333 %8.3f%8.3f%8.3f\n",water_hydrogen[0],water_hydrogen[1],water_hydrogen[2]);
                water_hydrogen[0] = wat216.xyz[2][0][m] + i0*wat216.cubic_size + edge_min[0] - w216_min[0];
                water_hydrogen[1] = wat216.xyz[2][1][m] + i1*wat216.cubic_size + edge_min[1] - w216_min[1];
                water_hydrogen[2] = wat216.xyz[2][2][m] + i2*wat216.cubic_size + edge_min[2] - w216_min[2];
                printf("ATOM   33 H2  HOH 333 %8.3f%8.3f%8.3f\n",water_hydrogen[0],water_hydrogen[1],water_hydrogen[2]);
            }
        }
    }}}
    return;
}
