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
def read_biomt_file(filename,rot_array,xyzcent_t_array):
   fp = open(filename,'r')
   line = fp.readline()
   if (line[0:6] == "REMARK") && (line[10:18] == "   BIOMT"):
      if nbiomt_t < 180:
         nrot = nbiomt_t / 3
         i    = nbiomt_t % 3
         rot[nrot][i][0] = 
         rot[nrot][i][1] = 
         rot[nrot][i][2] = 
      nbiomt_t++
    fp.close()

    
    l = 0;
    for (i = 0; i < 60; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                rot_array[l] = rot[i][j][k];
                l++;
            }
        }
    }
    
    l = 0;
    for (i = 0; i < 60; i++) {
        for (j = 0; j < 3; j++) {
            xyzcent_t_array[l] = xyzcent_t[i][j];
            l++;
        }
    }
}



void gen_entire(double *rot_array,double *xyzcent_t_array,double *xyz_x,double *xyz_y,double *xyz_z,double *xyz_x2,double *xyz_y2,double *xyz_z2,int n,int sn) {
    
    int i,j,k,l;
    double rot[60][3][3];
    double xyzcent_t[60][3];
    double xyz_rotated[3];
    l = 0;
    
    for (i = 0; i < 60; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                rot[i][j][k] = rot_array[l];
                l++;
            }
        }
    }
    
    l = 0;
    
    for (i = 0; i < 60; i++) {
        for (j = 0; j < 3; j++) {
            xyzcent_t[i][j] = xyzcent_t_array[l];
            l++;
        }
    }
    
    for (i = 0; i < n; i++) {
        xyz_rotated[0] = rot[sn][0][0]*(xyz_x[i]) + rot[sn][0][1]*(xyz_y[i]) + rot[sn][0][2]*(xyz_z[i]) + xyzcent_t[sn][0];
        xyz_rotated[1] = rot[sn][1][0]*(xyz_x[i]) + rot[sn][1][1]*(xyz_y[i]) + rot[sn][1][2]*(xyz_z[i]) + xyzcent_t[sn][1];
        xyz_rotated[2] = rot[sn][2][0]*(xyz_x[i]) + rot[sn][2][1]*(xyz_y[i]) + rot[sn][2][2]*(xyz_z[i]) + xyzcent_t[sn][2];
        xyz_x2[i] = xyz_rotated[0];
        xyz_y2[i] = xyz_rotated[1];
        xyz_z2[i] = xyz_rotated[2];
    }
}
