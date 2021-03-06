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
/* prepare and initialize the rotation matrices used for 2PLV and 4RHV */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef DEBUG
#define DEBUG 0
#endif

double matrix_determinant(double [3][3]);
void matrix_product(double [3][3],double [3][3],double [3][3]);
void make_unit_matrix(double *);
void matrix_inverse(double [3][3],double [3][3]);

void set_2plv_rotmat(double rotmat_2plv[60][3][3],double rotmat_4rhv[60][3][3],double rotmat_2plv_12[12][3][3],double invmat_2plv_12[12][3][3],double invmat_2plv[60][3][3],double invmat_4rhv[60][3][3]){
	
	double zero    = 0.0;
	double one     = 1.0;
	double half    = 1.0/2.0;
	double quarter = 1.0/4.0;
	double root5   = sqrt(5.0);
	int i,j,m;
	double det;
	
	double mat_A_2plv[3][3],mat_B_2plv[3][3],mat_C_2plv[3][3],mat_D_2plv[3][3];
	double mat_t[3][3]     = {0.,-1.,0., 1.,0.,0.,0.,0., 1.};
	double mat_t_inv[3][3] = {0., 1.,0.,-1.,0.,0.,0.,0., 1.};
	double mat_temp[3][3];
	
	// The rotation matrices of the old 2PLV file (The present 2PLV file is revised).
	/*
	 REMARK   7 EACH OF THE SIXTY ICOSAHEDRALLY-RELATED PROTOMERS CAN BE     2PLV 158
	 REMARK   7 GENERATED BY MULTIPLYING THE (CARTESIAN) ATOMIC              2PLV 159
	 REMARK   7 COORDINATES BY ONE OF SIXTY PRODUCT MATRICES.  THE PRODUCT   2PLV 160
	 REMARK   7 MATRICES HAVE THE FORM -                                     2PLV 161
	 REMARK   7                                                              2PLV 162
	 REMARK   7   MATRIX( J, K, L, M ) = ( A )**J * ( B )**K * ( C )**L *    2PLV 163
	 REMARK   7                          ( D )**M                            2PLV 164
	 REMARK   7                                                              2PLV 165
	 REMARK   7 WHERE J=(0 OR 1), K=(0 OR 1), L=(0, 1, OR 2), AND            2PLV 166
	 REMARK   7 M=(0, 1, 2, 3, OR 4), AND THE MATRIX MULTIPLICATIONS ARE     2PLV 167
	 REMARK   7 PERFORMED IN THE ORDER SPECIFIED.                            2PLV 168
	 REMARK   7                                                              2PLV 169
	 REMARK   7 ( D ) = FIVEFOLD  =                                          2PLV 170
	 REMARK   7  0.309017     -0.809017      0.500000                        2PLV 171
	 REMARK   7  0.809017      0.500000      0.309017                        2PLV 172
	 REMARK   7 -0.500000      0.309017      0.809017                        2PLV 173
	 REMARK   7                                                              2PLV 174
	 REMARK   7 ( C ) = THREEFOLD =                                          2PLV 175
	 REMARK   7 -0.309017     -0.809017      0.500000                        2PLV 176
	 REMARK   7  0.809017     -0.500000     -0.309017                        2PLV 177
	 REMARK   7  0.500000      0.309017      0.809017                        2PLV 178
	 REMARK   7                                                              2PLV 179
	 REMARK   7 ( B ) = TWOFOLD   =                                          2PLV 180
	 REMARK   7 -0.500000     -0.309017     -0.809017                        2PLV 181
	 REMARK   7 -0.309017     -0.809017      0.500000                        2PLV 182
	 REMARK   7 -0.809017      0.500000      0.309017                        2PLV 183
	 REMARK   7                                                              2PLV 184
	 REMARK   7 ( A ) = TWOFOLD   =                                          2PLV 185
	 REMARK   7 -0.809017     -0.500000      0.309017                        2PLV 186
	 REMARK   7 -0.500000      0.309017     -0.809017                        2PLV 187
	 REMARK   7  0.309017     -0.809017     -0.500000   
	 */
	
	// Rotation matrix ( D ) of the old 2PLV file (5fold axis)
	mat_D_2plv[0][0] =  (root5 - one)*quarter;
	mat_D_2plv[0][1] = -(root5 + one)*quarter;
	mat_D_2plv[0][2] = half;
	mat_D_2plv[1][0] =  (root5 + one)*quarter;
	mat_D_2plv[1][1] = half;
	mat_D_2plv[1][2] =  (root5 - one)*quarter;
	mat_D_2plv[2][0] = -half;
	mat_D_2plv[2][1] =  (root5 - one)*quarter;
	mat_D_2plv[2][2] =  (root5 + one)*quarter;
	// Rotation matrix ( C ) of the old 2plv file (3fold axis)
	mat_C_2plv[0][0] = -(root5 - one)*quarter; 
	mat_C_2plv[0][1] = -(root5 + one)*quarter;
	mat_C_2plv[0][2] = half;
	mat_C_2plv[1][0] =  (root5 + one)*quarter;
	mat_C_2plv[1][1] = -half;
	mat_C_2plv[1][2] = -(root5 - one)*quarter;
	mat_C_2plv[2][0] = half;
	mat_C_2plv[2][1] =  (root5 - one)*quarter;
	mat_C_2plv[2][2] =  (root5 + one)*quarter;
	// Rotation matrix ( B ) of the old 2plv file (2fold axis)
	mat_B_2plv[0][0] = -half;
	mat_B_2plv[0][1] = -(root5 - one)*quarter;
	mat_B_2plv[0][2] = -(root5 + one)*quarter;
	mat_B_2plv[1][0] = -(root5 - one)*quarter;
	mat_B_2plv[1][1] = -(root5 + one)*quarter;
	mat_B_2plv[1][2] = half;
	mat_B_2plv[2][0] = -(root5 + one)*quarter;
	mat_B_2plv[2][1] = half;
	mat_B_2plv[2][2] =  (root5 - one)*quarter;
	// Rotation matrix ( A ) of the old 2plv file (2fold axis)
	mat_A_2plv[0][0] = -(root5 + one)*quarter;
	mat_A_2plv[0][1] = -half;
	mat_A_2plv[0][2] =  (root5 - one)*quarter;
	mat_A_2plv[1][0] = -half;
	mat_A_2plv[1][1] =  (root5 - one)*quarter;
	mat_A_2plv[1][2] = -(root5 + one)*quarter;
	mat_A_2plv[2][0] =  (root5 - one)*quarter;
	mat_A_2plv[2][1] = -(root5 + one)*quarter;
	mat_A_2plv[2][2] = -half;
	
	if (DEBUG >= 3) {
		for (i=0; i<3; i++) {
			printf("D of 2plv = ");
			for (j=0; j<3; j++) {
				printf("%10.6f ",mat_D_2plv[i][j]);
			}
			printf("\n");
		}
		for (i=0; i<3; i++) {
			printf("C of 2plv = ");
			for (j=0; j<3; j++) {
				printf("%10.6f ",mat_C_2plv[i][j]);
			}
			printf("\n");
		}	for (i=0; i<3; i++) {
			printf("B of 2plv = ");
			for (j=0; j<3; j++) {
				printf("%10.6f ",mat_B_2plv[i][j]);
			}
			printf("\n");
		}	for (i=0; i<3; i++) {
			printf("A of 2plv = ");
			for (j=0; j<3; j++) {
				printf("%10.6f ",mat_A_2plv[i][j]);
			}
			printf("\n");
		}		
	}
		
	// The first matrix is unit matrix.
	rotmat_2plv[0][0][0] = one;
	rotmat_2plv[0][0][1] = zero;
	rotmat_2plv[0][0][2] = zero;
	rotmat_2plv[0][1][0] = zero;
	rotmat_2plv[0][1][1] = one;
	rotmat_2plv[0][1][2] = zero;
	rotmat_2plv[0][2][0] = zero;
	rotmat_2plv[0][2][1] = zero;
	rotmat_2plv[0][2][2] = one;
	for (i=1; i<5; i++) {
		matrix_product(mat_D_2plv, (double (*)[3])&rotmat_2plv[i-1][0][0], (double (*)[3])&rotmat_2plv[i][0][0]);
	}
	for (i=0; i<5; i++) {
		matrix_product(mat_C_2plv, (double (*)[3])&rotmat_2plv[i][0][0], (double (*)[3])&rotmat_2plv[i+5][0][0]);
	}
	for (i=5; i<10; i++) {
		matrix_product(mat_C_2plv, (double (*)[3])&rotmat_2plv[i][0][0], (double (*)[3])&rotmat_2plv[i+5][0][0]);
	}
	for (i=0; i<15; i++) {
		matrix_product(mat_B_2plv, (double (*)[3])&rotmat_2plv[i][0][0], (double (*)[3])&rotmat_2plv[i+15][0][0]);
	}
	for (i=0; i<30; i++) {
		matrix_product(mat_A_2plv, (double (*)[3])&rotmat_2plv[i][0][0], (double (*)[3])&rotmat_2plv[i+30][0][0]);
	}
    
    // The same as above, but for 1/12 cells
	rotmat_2plv_12[0][0][0] = one;
	rotmat_2plv_12[0][0][1] = zero;
	rotmat_2plv_12[0][0][2] = zero;
	rotmat_2plv_12[0][1][0] = zero;
	rotmat_2plv_12[0][1][1] = one;
	rotmat_2plv_12[0][1][2] = zero;
	rotmat_2plv_12[0][2][0] = zero;
	rotmat_2plv_12[0][2][1] = zero;
	rotmat_2plv_12[0][2][2] = one;
	for (i=1; i<2; i++) {
		matrix_product(mat_C_2plv, (double (*)[3])&rotmat_2plv_12[i-1][0][0], (double (*)[3])&rotmat_2plv_12[i][0][0]);
	}
	for (i=0; i<3; i++) {
		matrix_product(mat_B_2plv, (double (*)[3])&rotmat_2plv_12[i][0][0], (double (*)[3])&rotmat_2plv_12[i+3][0][0]);
	}
	for (i=0; i<6; i++) {
		matrix_product(mat_A_2plv, (double (*)[3])&rotmat_2plv_12[i][0][0], (double (*)[3])&rotmat_2plv_12[i+6][0][0]);
	}

    for(m=0;m<60;m++){
        matrix_inverse((double (*)[3])&rotmat_2plv[m][0][0],(double (*)[3])&invmat_2plv[m][0][0]);
    }
    for(m=0;m<12;m++){
        matrix_inverse((double (*)[3])&rotmat_2plv_12[m][0][0],(double (*)[3])&invmat_2plv_12[m][0][0]);
    }

	// The rotation matrices of the old 4RHV file (The present 4RHV file is revised).
	/*
	 REMARK   5                                                              4RHV 132
	 REMARK   5    LET P1 = COORDINATES OF ANY ATOM AS LISTED IN ENTRY       4RHV 133
	 REMARK   5                                                              4RHV 134	 
	 REMARK   5    1.                                                        4RHV 135
	 REMARK   5       TRNSF1 1    .500000  -.809017   .309017                4RHV 136
	 REMARK   5       TRNSF2 1    .809017   .309017  -.500000                4RHV 137
	 REMARK   5       TRNSF3 1    .309017   .500000   .809017                4RHV 138
	 REMARK   5                                                              4RHV 139
	 REMARK   5       APPLY THE FIVE-FOLD GIVEN BY TRNSF 1                   4RHV 140
	 REMARK   5                                                              4RHV 141
	 REMARK   5       P2 = TRNSF 1 * P1                                      4RHV 142
	 REMARK   5       P3 = TRNSF 1 * P2                                      4RHV 143
	 REMARK   5       P4 = TRNSF 1 * P3                                      4RHV 144
	 REMARK   5       P5 = TRNSF 1 * P4                                      4RHV 145
	 REMARK   5                                                              4RHV 146
	 REMARK   5       P1 THROUGH P5 CONSTITUTE AN ENTIRE PENTAMER,           4RHV 147
	 REMARK   5       CONSISTING OF FIVE PROTOMERIC UNITS.  ONE COMPLETE     4RHV 148
	 REMARK   5       PENTAMER SET OF COORDINATES WILL CONTAIN FIVE CHAINS   4RHV 149
	 REMARK   5       EACH OF VP1, VP2, VP3, AND VP4.                        4RHV 150
	 REMARK   5                                                              4RHV 151
	 REMARK   5    2.                                                        4RHV 152
	 REMARK   5                                                              4RHV 153
	 REMARK   5       TRNSF1 2    .309017  -.500000   .809017                4RHV 154
	 REMARK   5       TRNSF2 2   -.500000  -.809017  -.309017                4RHV 155
	 REMARK   5       TRNSF3 2    .809017  -.309017  -.500000                4RHV 156
	 REMARK   5                                                              4RHV 157
	 REMARK   5       APPLY THE TWO-FOLD GIVEN BY TRNSF 2 TO P1 THROUGH P5   4RHV 158
	 REMARK   5       TO GENERATE P6 THROUGH P10 FOR A SECOND ENTIRE         4RHV 159
	 REMARK   5       PENTAMER.                                              4RHV 160
	 REMARK   5                                                              4RHV 161
	 REMARK   5       P6 = TRNSF 2 * P1                                      4RHV 162
	 REMARK   5       P7 = TRNSF 2 * P2                                      4RHV 163
	 REMARK   5       P8 = TRNSF 2 * P3                                      4RHV 164
	 REMARK   5       P9 = TRNSF 2 * P4                                      4RHV 165
	 REMARK   5       P10 = TRNSF 2 * P5                                     4RHV 166
	 REMARK   5                                                              4RHV 167	 
	 REMARK   5    3.                                                        4RHV 168
	 REMARK   5                                                              4RHV 169
	 REMARK   5       TRNSF1 3   -.500000  -.809017   .309017                4RHV 170
	 REMARK   5       TRNSF2 3   -.809017   .309017  -.500000                4RHV 171
	 REMARK   5       TRNSF3 3    .309017  -.500000  -.809017                4RHV 172
	 REMARK   5                                                              4RHV 173
	 REMARK   5       APPLY THE TWO-FOLD GIVEN BY TRNSF 3 TO EACH OF THE     4RHV 174
	 REMARK   5       FIRST TEN PROTOMERS (WILL YIELD TOTAL OF TWENTY        4RHV 175
	 REMARK   5       UNITS, FOUR COMPLETE PENTAMERS).  FOR EXAMPLE,         4RHV 176
	 REMARK   5                                                              4RHV 177
	 REMARK   5       P11 = TRNSF 3 * P1                                     4RHV 178
	 REMARK   5       P14 = TRNSF 3 * P4                                     4RHV 179
	 REMARK   5       P16 = TRNSF 3 * P6                                     4RHV 180
	 REMARK   5       P20 = TRNSF 3 * P10                                    4RHV 181
	 REMARK   5                                                              4RHV 182	 
	 REMARK   5    4.                                                        4RHV 183
	 REMARK   5                                                              4RHV 184
	 REMARK   5       TRNSF1 4    .000000   .000000  1.000000                4RHV 185
	 REMARK   5       TRNSF2 4   1.000000   .000000   .000000                4RHV 186
	 REMARK   5       TRNSF3 4    .000000  1.000000   .000000                4RHV 187
	 REMARK   5       FIRST, APPLY THE THREE-FOLD GIVEN BY TRNSF 4 TO EACH   4RHV 189
	 REMARK   5       OF THE TWENTY PROTOMERS P1 THROUGH P20 TO GENERATE     4RHV 190
	 REMARK   5       PROTOMERS P21 THROUGH P40.  FOR EXAMPLE,               4RHV 191
	 REMARK   5                                                              4RHV 192
	 REMARK   5       P21 = TRNSF 4 * P2                                     4RHV 193
	 REMARK   5       P22 = TRNSF 4 * P2                                     4RHV 194
	 REMARK   5       P31 = TRNSF 4 * P11                                    4RHV 195
	 REMARK   5       P40 = TRNSF 4 * P20                                    4RHV 196
	 REMARK   5                                                              4RHV 197
	 REMARK   5       FINALLY, PROTOMERS P41 THROUGH P60 MAY BE GENERATED    4RHV 198
	 REMARK   5       BY APPLYING TRNSF 4 TO PROTOMERS P21 THROUGH P40.      4RHV 199
	 REMARK   5       FOR EXAMPLE,                                           4RHV 200
	 REMARK   5                                                              4RHV 201
	 REMARK   5       P41 = TRNSF 4 * P21                                    4RHV 202
	 REMARK   5       P42 = TRNSF 4 * P22                                    4RHV 203
	 REMARK   5       P51 = TRNSF 4 * P31                                    4RHV 204
	 REMARK   5       P60 = TRNSF 4 * P40                                    4RHV 205
	 REMARK   5                                                              4RHV 206
	 REMARK   5                                                              4RHV 207
	 REMARK   5       THIS YIELDS A TOTAL OF 60 PROTOMERS IN AN ICOSAHEDRAL  4RHV 208
	 REMARK   5       VIRION.                                                4RHV 209
	*/
	
	/* It was very difficult to decide which method should be used to generate the 60
	   matrices of the 4RHV coordinate system: There can be two possible methods.
	   1) Use the rotation matrices of TRANSF 1-4 that were written in the old 4RHV
	      file. In this method, the rotation matrices are written as,
	                R = F4^3 * F3^2 * F2^2 * F1^5
	      so that the style of cell number is different from that in 2PLV,
	      because in 2PLV, the rotation matrices are,
	                R = A^2 * B^2 * C^3 * D^5.
	      However, this method is royal to the 4RHV cell number style.
	   2) Rotate the coordinate system to be equal to the 2PLV coordinate system, then
	      apply the rotation matrices of 2PLV, and inversely rotate the coordinate
	      system to be equal to the 4RHV coordinate system. In this method, the style
	      of cell number is changed and is made to be equal to that in 2PLV.
	   After some consideration, I adopted the second method. I believe the use of one
	   cell number style (2PLV style of cell number) for all PDB files are important.
	   In the following, the rotation matrices of 2PLV are transformed as,
	              (new matrix) = T^-1 (2PLV matrix) T
	   and the results are used for the 4RHV coordinates. T is written as mat_t[3][3]
	   and T^-1 as mat_t_inv[3][3] in this program.
	 */
	
	for (m=0; m<60; m++) {
		matrix_product((double (*)[3])&rotmat_2plv[m][0][0],mat_t,mat_temp);
		matrix_product(mat_t_inv,mat_temp,(double (*)[3])&rotmat_4rhv[m][0][0]);
	}
	
	// Debug lines
	if (DEBUG >= 3) {
		for (i=0; i<60; i++) {
		printf("2plv %d:\n",i);
		printf("%7.3f%7.3f%7.3f\n",rotmat_2plv[i][0][0],rotmat_2plv[i][0][1],rotmat_2plv[i][0][2]);
		printf("%7.3f%7.3f%7.3f\n",rotmat_2plv[i][1][0],rotmat_2plv[i][1][1],rotmat_2plv[i][1][2]);
		printf("%7.3f%7.3f%7.3f\n",rotmat_2plv[i][2][0],rotmat_2plv[i][2][1],rotmat_2plv[i][2][2]);
		printf("\n");
		}
	}
	if (DEBUG >= 1) {
		for (i=0; i<60; i++) {
			det = matrix_determinant((double (*)[3])&rotmat_2plv[i][0][0]);
			if (fabs(1.0 - det) > 0.0001) {
				printf("The determinant of matrix %d of 2PLV is not 1, but is %f\n",i,det);
			}
		}
		for (m=0; m<60; m++) {
			det = matrix_determinant((double (*)[3])&rotmat_4rhv[m][0][0]);
			if (fabs(1.0 - det) > 0.0001) {
				printf("The determinant of matrix %d of 4RHV is not 1, but is %f\n",m,det);
			}
		}
	}
    
    for(m=0;m<60;m++){
        matrix_inverse((double (*)[3])&rotmat_4rhv[m][0][0],(double (*)[3])&invmat_4rhv[m][0][0]);
    }
}
