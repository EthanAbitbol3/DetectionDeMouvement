/* -------------------------------- */
/* --- morpho_ouverture_swp64.c --- */
/* -------------------------------- */

/*
 * Copyright (c) 2004 - 2013, Lionel Lacassagne, All rights reserved
 * University of Paris Sud, Laboratoire de Recherche en Informatique
*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "nrtype.h"
#include "nrdef.h"
#include "nrutil.h"
//#include "nrutil_ext.h" // printfM8

//#include "sequence.h"

#include "swp.h"  // left right
#include "morpho.h"

#include "morpho_erosion_swp64.h"
#include "morpho_dilatation_swp64.h"
#include "morpho_ouverture_swp64.h"
// -----------------------------------------------------------------------------------------------
void line_ouverture3_ui64matrix_swp64_fusion_basic(uint64 **X, int i, int j0, int j1, uint64 **Y)
// -----------------------------------------------------------------------------------------------
{
    /*
        Version basique de fusion
    */
    uint64 t00, t01, t02; // colonnes temporaire
    uint64 t10, t11, t12;
    uint64 t20, t21, t22; 

    uint64 l0,r0,l1,r1,l2,r2;

    int j;
    for(j=j0; j<=j1; j++)
    {
        /*
            Calcule des cases temporaire 
        */
        t00 = MIN9( LEFT_64(X[i-2][j-2], X[i-2][j-1]),X[i-2][j-1],RIGHT_64(X[i-2][j-1], X[i-2][j]),
                    LEFT_64(X[i-1][j-2], X[i-1][j-1]),X[i-1][j-1],RIGHT_64(X[i-1][j-1], X[i-1][j]),
                    LEFT_64(X[i+0][j-2], X[i+0][j-1]),X[i+0][j-1],RIGHT_64(X[i+0][j-1], X[i+0][j]));

        t01 = MIN9( LEFT_64(X[i-2][j-1], X[i-2][j+0]),X[i-2][j+0],RIGHT_64(X[i-2][j-0], X[i-2][j+1]),
                    LEFT_64(X[i-1][j-1], X[i-1][j+0]),X[i-1][j+0],RIGHT_64(X[i-1][j-0], X[i-1][j+1]),
                    LEFT_64(X[i+0][j-1], X[i+0][j+0]),X[i+0][j+0],RIGHT_64(X[i+0][j-0], X[i+0][j+1]));
       
        t02 = MIN9( LEFT_64(X[i-2][j+0], X[i-2][j+1]),X[i-2][j+1],RIGHT_64(X[i-2][j+1], X[i-2][j+2]),
                    LEFT_64(X[i-1][j+0], X[i-1][j+1]),X[i-1][j+1],RIGHT_64(X[i-1][j+1], X[i-1][j+2]),
                    LEFT_64(X[i+0][j+0], X[i+0][j+1]),X[i+0][j+1],RIGHT_64(X[i+0][j+1], X[i+0][j+2])); 

        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t10 = MIN9( LEFT_64(X[i-1][j-2], X[i-1][j-1]),X[i-1][j-1],RIGHT_64(X[i-1][j-1], X[i-1][j+0]),
                    LEFT_64(X[i+0][j-2], X[i+0][j-1]),X[i+0][j-1],RIGHT_64(X[i+0][j-1], X[i+0][j+0]),
                    LEFT_64(X[i+1][j-2], X[i+1][j-1]),X[i+1][j-1],RIGHT_64(X[i+1][j-1], X[i+1][j+0]));
        
        t11 = MIN9( LEFT_64(X[i-1][j-1], X[i-1][j+0]),X[i-1][j+0],RIGHT_64(X[i-1][j+0], X[i-1][j+1]),
                    LEFT_64(X[i+0][j-1], X[i+0][j+0]),X[i+0][j+0],RIGHT_64(X[i+0][j+0], X[i+0][j+1]),
                    LEFT_64(X[i+1][j-1], X[i+1][j+0]),X[i+1][j+0],RIGHT_64(X[i+1][j+0], X[i+1][j+1]));
        
        t12 = MIN9( LEFT_64(X[i-1][j+0], X[i-1][j+1]),X[i-1][j+1],RIGHT_64(X[i-1][j+1], X[i-1][j+2]),
                    LEFT_64(X[i+0][j+0], X[i+0][j+1]),X[i+0][j+1],RIGHT_64(X[i+0][j+1], X[i+0][j+2]),
                    LEFT_64(X[i+1][j+0], X[i+1][j+1]),X[i+1][j+1],RIGHT_64(X[i+1][j+1], X[i+1][j+2])); 
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t20 = MIN9( LEFT_64(X[i+0][j-2], X[i+0][j-1]),X[i+0][j-1],RIGHT_64(X[i+0][j-1], X[i+0][j+0]),
                    LEFT_64(X[i+1][j-2], X[i+1][j-1]),X[i+1][j-1],RIGHT_64(X[i+1][j-1], X[i+1][j+0]),
                    LEFT_64(X[i+2][j-2], X[i+2][j-1]),X[i+2][j-1],RIGHT_64(X[i+2][j-1], X[i+2][j+0]));
        
        t21 = MIN9( LEFT_64(X[i+0][j-1], X[i+0][j+0]),X[i+0][j+0],RIGHT_64(X[i+0][j+0], X[i+0][j+1]),
                    LEFT_64(X[i+1][j-1], X[i+1][j+0]),X[i+1][j+0],RIGHT_64(X[i+1][j+0], X[i+1][j+1]),
                    LEFT_64(X[i+2][j-1], X[i+2][j+0]),X[i+2][j+0],RIGHT_64(X[i+2][j+0], X[i+2][j+1]));
        
        t22 = MIN9( LEFT_64(X[i+0][j+0], X[i+0][j+1]),X[i+0][j+1],RIGHT_64(X[i+0][j+1], X[i+0][j+2]),
                    LEFT_64(X[i+1][j+0], X[i+1][j+1]),X[i+1][j+1],RIGHT_64(X[i+1][j+1], X[i+1][j+2]),
                    LEFT_64(X[i+2][j+0], X[i+2][j+1]),X[i+2][j+1],RIGHT_64(X[i+2][j+1], X[i+2][j+2]));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22);            
        // ----------------------------------
        // calcule de la case finale 
        Y[i][j] = MAX9( l0, t01, r0,
                     l1, t11, r1,
                     l2, t21, r2
                     );
    }
}
// --------------------------------------------------------------------------------------------
void line_ouverture3_ui64matrix_swp64_fusion_red(uint64 **X, int i, int j0, int j1, uint64 **Y)
// --------------------------------------------------------------------------------------------
{
    /*
        Version reduction de colonne de fusion ouverture
    */
    // Var des 5 lignes 
    uint64 *a0 = X[i-2]; // colonnes X
    uint64 *a1 = X[i-1];
    uint64 *a2 = X[i+0];
    uint64 *a3 = X[i+1]; 
    uint64 *a4 = X[i+2];

    uint64 *y = Y[i];

    uint64 l0,r0,l1,r1,l2,r2;

    uint64 t00, t01, t02; // colonnes temporaire
    uint64 t10, t11, t12;
    uint64 t20, t21, t22; 

    uint64 c00, c01, c02, c03, c04;
    uint64 c10, c11, c12, c13, c14;
    uint64 c20, c21, c22, c23, c24;
    
    int j = j0;

    c00 = MIN3(a0[j-2], a1[j-2], a2[j-2]); c10 = MIN3(a1[j-2], a2[j-2], a3[j-2]); c20 = MIN3(a2[j-2], a3[j-2], a4[j-2]);
    c01 = MIN3(a0[j-1], a1[j-1], a2[j-1]); c11 = MIN3(a1[j-1], a2[j-1], a3[j-1]); c21 = MIN3(a2[j-1], a3[j-1], a4[j-1]);
    c02 = MIN3(a0[j+0], a1[j+0], a2[j+0]); c12 = MIN3(a1[j+0], a2[j+0], a3[j+0]); c22 = MIN3(a2[j+0], a3[j+0], a4[j+0]);
    c03 = MIN3(a0[j+1], a1[j+1], a2[j+1]); c13 = MIN3(a1[j+1], a2[j+1], a3[j+1]); c23 = MIN3(a2[j+1], a3[j+1], a4[j+1]);
    

    t00 = MIN3(LEFT_64(c00, c01),c01,RIGHT_64(c01, c02)); 
    t01 = MIN3(LEFT_64(c01, c02),c02,RIGHT_64(c02, c03));

    t10 = MIN3(LEFT_64(c10, c11),c11,RIGHT_64(c11, c12));
    t11 = MIN3(LEFT_64(c11, c12),c12,RIGHT_64(c12, c13));

    t20 = MIN3(LEFT_64(c20, c21),c21,RIGHT_64(c21, c22));
    t21 = MIN3(LEFT_64(c21, c22),c22,RIGHT_64(c22, c23));

    for(j=j0; j<=j1; j++)
    {
        /*
            Calcule des cases temporaire 
        */
        c04 = MIN3(a0[j+2], a1[j+2], a2[j+2]); c14 = MIN3(a1[j+2], a2[j+2], a3[j+2]); c24 = MIN3(a2[j+2], a3[j+2], a4[j+2]);
        
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        // calcule de la case finale 
        Y[i][j] = MAX9( 
                        l0, t01, r0,
                        l1, t11, r1,
                        l2, t21, r2
                    );
        
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;
    }
}
// -------------------------------------------------------------------------------------------------
void line_ouverture3_ui64matrix_swp64_fusion_ilu3_red(uint64 **X, int i, int j0, int j1, uint64 **Y)
// -------------------------------------------------------------------------------------------------
{
    /*
        Version avec deroulage de boucle interne de degre 3 et
        reduction de colonnes de fusion 
    */

    // Var des 5 lignes 
    uint64 *a0 = X[i-2]; // colonnes X
    uint64 *a1 = X[i-1];
    uint64 *a2 = X[i+0];
    uint64 *a3 = X[i+1]; 
    uint64 *a4 = X[i+2];

    uint64 *Yi = Y[i];
    uint64 y;

    uint64 l0,r0,l1,r1,l2,r2;

    uint64 t00, t01, t02; // colonnes temporaire
    uint64 t10, t11, t12;
    uint64 t20, t21, t22; 

    // cases de X 
    uint64 a00, a01, a02, a03, a04;
    uint64 a10, a11, a12, a13, a14;
    uint64 a20, a21, a22, a23, a24;
    uint64 a30, a31, a32, a33, a34;
    uint64 a40, a41, a42, a43, a44;

    uint64 c00, c01, c02, c03, c04;
    uint64 c10, c11, c12, c13, c14;
    uint64 c20, c21, c22, c23, c24;
    
    int j = j0;

    // load 
    a00 = load1(a0, j-2); a10 = load1(a1, j-2); a20 = load1(a2, j-2); a30 = load1(a3, j-2); a40 = load1(a4, j-2);
    a01 = load1(a0, j-1); a11 = load1(a1, j-1); a21 = load1(a2, j-1); a31 = load1(a3, j-1); a41 = load1(a4, j-1);
    a02 = load1(a0, j  ); a12 = load1(a1, j  ); a22 = load1(a2, j  ); a32 = load1(a3, j  ); a42 = load1(a4, j  );
    a03 = load1(a0, j+1); a13 = load1(a1, j+1); a23 = load1(a2, j+1); a33 = load1(a3, j+1); a43 = load1(a4, j+1);

    // calcule des colonnes 
    c00 = MIN3( a00, a10, a20);
    c01 = MIN3( a01, a11, a21);
    c02 = MIN3( a02, a12, a22);
    c03 = MIN3( a03, a13, a23);

    c10 = MIN3( a10, a20, a30);
    c11 = MIN3( a11, a21, a31);
    c12 = MIN3( a12, a22, a32);
    c13 = MIN3( a13, a23, a33);

    c20 = MIN3(a20, a30, a40);
    c21 = MIN3(a21, a31, a41);
    c22 = MIN3(a22, a32, a42);
    c23 = MIN3(a23, a33, a43);

    t00 = MIN3(LEFT_64(c00, c01),c01,RIGHT_64(c01, c02)); 
    t01 = MIN3(LEFT_64(c01, c02),c02,RIGHT_64(c02, c03));

    t10 = MIN3(LEFT_64(c10, c11),c11,RIGHT_64(c11, c12));
    t11 = MIN3(LEFT_64(c11, c12),c12,RIGHT_64(c12, c13));

    t20 = MIN3(LEFT_64(c20, c21),c21,RIGHT_64(c21, c22));
    t21 = MIN3(LEFT_64(c21, c22),c22,RIGHT_64(c22, c23));
    
    // determination de l'intervalle de j
    int r = ((j1 - 2) % 3);
    for (j = j0; j <= ((j1 - 2) - r); j += 3){
        // ----------------------------------LOOP 1-----------------------------------------------
        // ----------------------------------
        // load
        // ----------------------------------
        a04 = load1(a0, j+2); a14 = load1(a1, j+2); a24 = load1(a2, j+2); a34 = load1(a3, j+2); a44 = load1(a4, j+2);
        // ----------------------------------
        // calc
        // ----------------------------------
        c04 = MIN3( a04, a14, a24);
        c14 = MIN3( a14, a24, a34);
        c24 = MIN3(a24, a34, a44);
        //-----------------------------------
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        // calcule de la case finale 
        y = MAX9( 
                        l0, t01, r0,
                        l1, t11, r1,
                        l2, t21, r2
                    );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi, j, y);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;


        // ----------------------------------LOOP 2-----------------------------------------------
        // ----------------------------------
        // load
        // ----------------------------------
        a04 = load1(a0, j+3); a14 = load1(a1, j+3); a24 = load1(a2, j+3); a34 = load1(a3, j+3); a44 = load1(a4, j+3);
        // ----------------------------------
        // calc
        // ----------------------------------
        c04 = MIN3( a04, a14, a24);
        c14 = MIN3( a14, a24, a34);
        c24 = MIN3(a24, a34, a44);
        //-----------------------------------
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        // calcule de la case finale 
        y = MAX9( 
                        l0, t01, r0,
                        l1, t11, r1,
                        l2, t21, r2
                    );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi, j+1, y);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;

        // ----------------------------------LOOP 3-----------------------------------------------
        // ----------------------------------
        // load
        // ----------------------------------
        a04 = load1(a0, j+4); a14 = load1(a1, j+4); a24 = load1(a2, j+4); a34 = load1(a3, j+4); a44 = load1(a4, j+4);
        // ----------------------------------
        // calc
        // ----------------------------------
        c04 = MIN3( a04, a14, a24);
        c14 = MIN3( a14, a24, a34);
        c24 = MIN3(a24, a34, a44);
        //-----------------------------------
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        // calcule de la case finale 
        y = MAX9( 
                        l0, t01, r0,
                        l1, t11, r1,
                        l2, t21, r2
                    );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi, j+2, y);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;
    }

    // epilogue
    j = ((j1-2)-r);
    j++;
    // load 
    a00 = load1(a0, j-2); a10 = load1(a1, j-2); a20 = load1(a2, j-2); a30 = load1(a3, j-2); a40 = load1(a4, j-2);
    a01 = load1(a0, j-1); a11 = load1(a1, j-1); a21 = load1(a2, j-1); a31 = load1(a3, j-1); a41 = load1(a4, j-1);
    a02 = load1(a0, j  ); a12 = load1(a1, j  ); a22 = load1(a2, j  ); a32 = load1(a3, j  ); a42 = load1(a4, j  );
    a03 = load1(a0, j+1); a13 = load1(a1, j+1); a23 = load1(a2, j+1); a33 = load1(a3, j+1); a43 = load1(a4, j+1);

    // calcule des colonnes 
    c00 = MIN3( a00, a10, a20);
    c01 = MIN3( a01, a11, a21);
    c02 = MIN3( a02, a12, a22);
    c03 = MIN3( a03, a13, a23);

    c10 = MIN3( a10, a20, a30);
    c11 = MIN3( a11, a21, a31);
    c12 = MIN3( a12, a22, a32);
    c13 = MIN3( a13, a23, a33);

    c20 = MIN3(a20, a30, a40);
    c21 = MIN3(a21, a31, a41);
    c22 = MIN3(a22, a32, a42);
    c23 = MIN3(a23, a33, a43);

    t00 = MIN3(LEFT_64(c00, c01),c01,RIGHT_64(c01, c02)); 
    t01 = MIN3(LEFT_64(c01, c02),c02,RIGHT_64(c02, c03));

    t10 = MIN3(LEFT_64(c10, c11),c11,RIGHT_64(c11, c12));
    t11 = MIN3(LEFT_64(c11, c12),c12,RIGHT_64(c12, c13));

    t20 = MIN3(LEFT_64(c20, c21),c21,RIGHT_64(c21, c22));
    t21 = MIN3(LEFT_64(c21, c22),c22,RIGHT_64(c22, c23));

    switch (r+2)
    {
        case 4:
        // ----------------------------------
        // load
        // ----------------------------------
        a04 = load1(a0, j+2); a14 = load1(a1, j+2); a24 = load1(a2, j+2); a34 = load1(a3, j+2); a44 = load1(a4, j+2);
        // ----------------------------------
        // calc
        // ----------------------------------
        c04 = MIN3( a04, a14, a24);
        c14 = MIN3( a14, a24, a34);
        c24 = MIN3(a24, a34, a44);
        //-----------------------------------
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        // calcule de la case finale 
        y = MAX9( 
                        l0, t01, r0,
                        l1, t11, r1,
                        l2, t21, r2
                    );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi, j, y);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;
        j++;

        case 3:
        // ----------------------------------
        // load
        // ----------------------------------
        a04 = load1(a0, j+2); a14 = load1(a1, j+2); a24 = load1(a2, j+2); a34 = load1(a3, j+2); a44 = load1(a4, j+2);
        // ----------------------------------
        // calc
        // ----------------------------------
        c04 = MIN3( a04, a14, a24);
        c14 = MIN3( a14, a24, a34);
        c24 = MIN3(a24, a34, a44);
        //-----------------------------------
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        // calcule de la case finale 
        y = MAX9( 
                        l0, t01, r0,
                        l1, t11, r1,
                        l2, t21, r2
                    );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi, j, y);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;
        j++;

        case 2:
        // ----------------------------------
        // load
        // ----------------------------------
        a04 = load1(a0, j+2); a14 = load1(a1, j+2); a24 = load1(a2, j+2); a34 = load1(a3, j+2); a44 = load1(a4, j+2);
        // ----------------------------------
        // calc
        // ----------------------------------
        c04 = MIN3( a04, a14, a24);
        c14 = MIN3( a14, a24, a34);
        c24 = MIN3(a24, a34, a44);
        //-----------------------------------
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        // calcule de la case finale 
        y = MAX9( 
                        l0, t01, r0,
                        l1, t11, r1,
                        l2, t21, r2
                    );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi, j, y);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;
        j++;

        case 1:
        // ----------------------------------
        // load
        // ----------------------------------
        a04 = load1(a0, j+2); a14 = load1(a1, j+2); a24 = load1(a2, j+2); a34 = load1(a3, j+2); a44 = load1(a4, j+2);
        // ----------------------------------
        // calc
        // ----------------------------------
        c04 = MIN3( a04, a14, a24);
        c14 = MIN3( a14, a24, a34);
        c24 = MIN3(a24, a34, a44);
        //-----------------------------------
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        // calcule de la case finale 
        y = MAX9( 
                        l0, t01, r0,
                        l1, t11, r1,
                        l2, t21, r2
                    );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi, j, y);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;
        j++;

        default:
            break;
    }
}
// -------------------------------------------------------------------------------------------------
void line_ouverture3_ui64matrix_swp64_fusion_elu2_red(uint64 **X, int i, int j0, int j1, uint64 **Y)
// -------------------------------------------------------------------------------------------------
{
    /*
        Version avec deroulage de boucle externe de degree 2 
        et reduction de colonnes
    */
   // Var des 5 lignes 
    uint64 *a0 = X[i-2]; // colonnes X
    uint64 *a1 = X[i-1];
    uint64 *a2 = X[i+0];
    uint64 *a3 = X[i+1]; 
    uint64 *a4 = X[i+2];
    uint64 *a5 = X[i+3]; //elu 2

    uint64 *Yi0 = Y[i];
    uint64 *Yi1 = Y[i+1];
    uint64 y0,y1;

    uint64 l0,r0,l1,r1,l2,r2,l3,r3;

    uint64 t00, t01, t02; // colonnes temporaire
    uint64 t10, t11, t12;
    uint64 t20, t21, t22; 
    uint64 t30, t31, t32; // elu2

    // cases de X 
    uint64 a00, a01, a02, a03, a04;
    uint64 a10, a11, a12, a13, a14;
    uint64 a20, a21, a22, a23, a24;
    uint64 a30, a31, a32, a33, a34;
    uint64 a40, a41, a42, a43, a44;
    uint64 a50, a51, a52, a53, a54;

    uint64 c00, c01, c02, c03, c04;
    uint64 c10, c11, c12, c13, c14;
    uint64 c20, c21, c22, c23, c24;
    uint64 c30, c31, c32, c33, c34;

    int j = j0;

    // load 
    a00 = load1(a0, j-2); a10 = load1(a1, j-2); a20 = load1(a2, j-2); a30 = load1(a3, j-2); a40 = load1(a4, j-2);
    a01 = load1(a0, j-1); a11 = load1(a1, j-1); a21 = load1(a2, j-1); a31 = load1(a3, j-1); a41 = load1(a4, j-1);
    a02 = load1(a0, j  ); a12 = load1(a1, j  ); a22 = load1(a2, j  ); a32 = load1(a3, j  ); a42 = load1(a4, j  );
    a03 = load1(a0, j+1); a13 = load1(a1, j+1); a23 = load1(a2, j+1); a33 = load1(a3, j+1); a43 = load1(a4, j+1);
    // elu 2
    a50 = load1(a5, j-2);
    a51 = load1(a5, j-1);
    a52 = load1(a5, j  );
    a53 = load1(a5, j+1);

    // calcule des colonnes 
    c00 = MIN3( a00, a10, a20);
    c01 = MIN3( a01, a11, a21);
    c02 = MIN3( a02, a12, a22);
    c03 = MIN3( a03, a13, a23);

    c10 = MIN3( a10, a20, a30);
    c11 = MIN3( a11, a21, a31);
    c12 = MIN3( a12, a22, a32);
    c13 = MIN3( a13, a23, a33);

    c20 = MIN3(a20, a30, a40);
    c21 = MIN3(a21, a31, a41);
    c22 = MIN3(a22, a32, a42);
    c23 = MIN3(a23, a33, a43);

    //elu 2
    c30 = MIN3(a30, a40, a50);
    c31 = MIN3(a31, a41, a51);
    c32 = MIN3(a32, a42, a52);
    c33 = MIN3(a33, a43, a53);

    t00 = MIN3(LEFT_64(c00, c01),c01,RIGHT_64(c01, c02)); 
    t01 = MIN3(LEFT_64(c01, c02),c02,RIGHT_64(c02, c03));

    t10 = MIN3(LEFT_64(c10, c11),c11,RIGHT_64(c11, c12));
    t11 = MIN3(LEFT_64(c11, c12),c12,RIGHT_64(c12, c13));

    t20 = MIN3(LEFT_64(c20, c21),c21,RIGHT_64(c21, c22));
    t21 = MIN3(LEFT_64(c21, c22),c22,RIGHT_64(c22, c23));

    //elu 2
    t30 = MIN3(LEFT_64(c30, c31),c31,RIGHT_64(c31, c32));
    t31 = MIN3(LEFT_64(c31, c32),c32,RIGHT_64(c32, c33));

    for (j = j0; j <= j1; j++){
        // load
        // ----------------------------------
        a04 = load1(a0, j+2); a14 = load1(a1, j+2); a24 = load1(a2, j+2); a34 = load1(a3, j+2); a44 = load1(a4, j+2);
        a54 = load1(a5, j+2);
        // ----------------------------------
        // calc
        // ----------------------------------
        c04 = MIN3( a04, a14, a24);
        c14 = MIN3( a14, a24, a34);
        c24 = MIN3(a24, a34, a44);
        c34 = MIN3( a34, a44, a54);
        //-----------------------------------
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        //elu 2
        t32 = MIN3(LEFT_64(c32, c33),c33,RIGHT_64(c33, c34));
        l3 = LEFT_64(t30,t31);
        r3 = RIGHT_64(t31,t32); 

        // calcule de la case finale 
        y0 = MAX9( 
                        l0, t01, r0,
                        l1, t11, r1,
                        l2, t21, r2
                    );
        
        y1 = MAX9( 
                        l1, t11, r1,
                        l2, t21, r2,
                        l3, t31, r3
                    );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi0, j, y0);
        store1(Yi1, j, y1);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;

        c32 = c33;
        c33 = c34;

        t30 = t31;
        t31 = t32;

    }

}
// ------------------------------------------------------------------------------------------------------
void line_ouverture3_ui64matrix_swp64_fusion_ilu3_elu2_red(uint64 **X, int i, int j0, int j1, uint64 **Y)
// ------------------------------------------------------------------------------------------------------
{
    /*
        Version avec deroulage de boucle interne de degre 3,
        deroulage de boucle interne de degre 2 et reduction de colonnes de fusion 
    */
   // Var des 5 lignes 
    uint64 *a0 = X[i-2]; // colonnes X
    uint64 *a1 = X[i-1];
    uint64 *a2 = X[i+0];
    uint64 *a3 = X[i+1]; 
    uint64 *a4 = X[i+2];
    uint64 *a5 = X[i+3]; //elu 2

    uint64 *Yi0 = Y[i];
    uint64 *Yi1 = Y[i+1];
    uint64 y0,y1;

    uint64 l0,r0,l1,r1,l2,r2,l3,r3;

    uint64 t00, t01, t02; // colonnes temporaire
    uint64 t10, t11, t12;
    uint64 t20, t21, t22; 
    uint64 t30, t31, t32; // elu2

    // cases de X 
    uint64 a00, a01, a02, a03, a04;
    uint64 a10, a11, a12, a13, a14;
    uint64 a20, a21, a22, a23, a24;
    uint64 a30, a31, a32, a33, a34;
    uint64 a40, a41, a42, a43, a44;
    uint64 a50, a51, a52, a53, a54;

    uint64 c00, c01, c02, c03, c04;
    uint64 c10, c11, c12, c13, c14;
    uint64 c20, c21, c22, c23, c24;
    uint64 c30, c31, c32, c33, c34;

    int j = j0;

    // load 
    a00 = load1(a0, j-2); a10 = load1(a1, j-2); a20 = load1(a2, j-2); a30 = load1(a3, j-2); a40 = load1(a4, j-2);
    a01 = load1(a0, j-1); a11 = load1(a1, j-1); a21 = load1(a2, j-1); a31 = load1(a3, j-1); a41 = load1(a4, j-1);
    a02 = load1(a0, j  ); a12 = load1(a1, j  ); a22 = load1(a2, j  ); a32 = load1(a3, j  ); a42 = load1(a4, j  );
    a03 = load1(a0, j+1); a13 = load1(a1, j+1); a23 = load1(a2, j+1); a33 = load1(a3, j+1); a43 = load1(a4, j+1);
    // elu 2
    a50 = load1(a5, j-2);
    a51 = load1(a5, j-1);
    a52 = load1(a5, j  );
    a53 = load1(a5, j+1);

    // calcule des colonnes 
    c00 = MIN3( a00, a10, a20);
    c01 = MIN3( a01, a11, a21);
    c02 = MIN3( a02, a12, a22);
    c03 = MIN3( a03, a13, a23);

    c10 = MIN3( a10, a20, a30);
    c11 = MIN3( a11, a21, a31);
    c12 = MIN3( a12, a22, a32);
    c13 = MIN3( a13, a23, a33);

    c20 = MIN3(a20, a30, a40);
    c21 = MIN3(a21, a31, a41);
    c22 = MIN3(a22, a32, a42);
    c23 = MIN3(a23, a33, a43);

    //elu 2
    c30 = MIN3(a30, a40, a50);
    c31 = MIN3(a31, a41, a51);
    c32 = MIN3(a32, a42, a52);
    c33 = MIN3(a33, a43, a53);

    t00 = MIN3(LEFT_64(c00, c01),c01,RIGHT_64(c01, c02)); 
    t01 = MIN3(LEFT_64(c01, c02),c02,RIGHT_64(c02, c03));

    t10 = MIN3(LEFT_64(c10, c11),c11,RIGHT_64(c11, c12));
    t11 = MIN3(LEFT_64(c11, c12),c12,RIGHT_64(c12, c13));

    t20 = MIN3(LEFT_64(c20, c21),c21,RIGHT_64(c21, c22));
    t21 = MIN3(LEFT_64(c21, c22),c22,RIGHT_64(c22, c23));

    //elu 2
    t30 = MIN3(LEFT_64(c30, c31),c31,RIGHT_64(c31, c32));
    t31 = MIN3(LEFT_64(c31, c32),c32,RIGHT_64(c32, c33));

    // determination de l'intervalle de j
    int r = ((j1 - 2) % 3);
    for (j = j0; j <= ((j1 - 2) - r); j += 3){
        // ----------------------------------LOOP 1-----------------------------------------------
        // load
        // ----------------------------------
        a04 = load1(a0, j+2); a14 = load1(a1, j+2); a24 = load1(a2, j+2); a34 = load1(a3, j+2); a44 = load1(a4, j+2);
        a54 = load1(a5, j+2);
        // ----------------------------------
        // calc
        // ----------------------------------
        c04 = MIN3( a04, a14, a24);
        c14 = MIN3( a14, a24, a34);
        c24 = MIN3(a24, a34, a44);
        c34 = MIN3( a34, a44, a54);
        //-----------------------------------
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        //elu 2
        t32 = MIN3(LEFT_64(c32, c33),c33,RIGHT_64(c33, c34));
        l3 = LEFT_64(t30,t31);
        r3 = RIGHT_64(t31,t32); 

        // calcule de la case finale 
        y0 = MAX9( 
                        l0, t01, r0,
                        l1, t11, r1,
                        l2, t21, r2
                    );
        
        y1 = MAX9( 
                        l1, t11, r1,
                        l2, t21, r2,
                        l3, t31, r3
                    );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi0, j, y0);
        store1(Yi1, j, y1);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;

        c32 = c33;
        c33 = c34;

        t30 = t31;
        t31 = t32;


        // ----------------------------------LOOP 2-----------------------------------------------
        // ----------------------------------
        // load
        // ----------------------------------
        a04 = load1(a0, j+3); a14 = load1(a1, j+3); a24 = load1(a2, j+3); a34 = load1(a3, j+3); a44 = load1(a4, j+3);
        a54 = load1(a5, j+3);
        // ----------------------------------
        // calc
        // ----------------------------------
        c04 = MIN3( a04, a14, a24);
        c14 = MIN3( a14, a24, a34);
        c24 = MIN3(a24, a34, a44);
        c34 = MIN3( a34, a44, a54);
        //-----------------------------------
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        //elu 2
        t32 = MIN3(LEFT_64(c32, c33),c33,RIGHT_64(c33, c34));
        l3 = LEFT_64(t30,t31);
        r3 = RIGHT_64(t31,t32); 

        // calcule de la case finale 
        y0 = MAX9( 
                        l0, t01, r0,
                        l1, t11, r1,
                        l2, t21, r2
                    );
        
        y1 = MAX9( 
                        l1, t11, r1,
                        l2, t21, r2,
                        l3, t31, r3
                    );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi0, j+1, y0);
        store1(Yi1, j+1, y1);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;

        c32 = c33;
        c33 = c34;

        t30 = t31;
        t31 = t32;

        // ----------------------------------LOOP 3-----------------------------------------------
        // ----------------------------------
        // load
        // ----------------------------------
        a04 = load1(a0, j+4); a14 = load1(a1, j+4); a24 = load1(a2, j+4); a34 = load1(a3, j+4); a44 = load1(a4, j+4);
        a54 = load1(a5, j+4);
        // ----------------------------------
        // calc
        // ----------------------------------
        c04 = MIN3( a04, a14, a24);
        c14 = MIN3( a14, a24, a34);
        c24 = MIN3(a24, a34, a44);
        c34 = MIN3( a34, a44, a54);
        //-----------------------------------
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        //elu 2
        t32 = MIN3(LEFT_64(c32, c33),c33,RIGHT_64(c33, c34));
        l3 = LEFT_64(t30,t31);
        r3 = RIGHT_64(t31,t32); 

        // calcule de la case finale 
        y0 = MAX9( 
                        l0, t01, r0,
                        l1, t11, r1,
                        l2, t21, r2
                    );
        
        y1 = MAX9( 
                        l1, t11, r1,
                        l2, t21, r2,
                        l3, t31, r3
                    );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi0, j+2, y0);
        store1(Yi1, j+2, y1);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;

        c32 = c33;
        c33 = c34;

        t30 = t31;
        t31 = t32;
    }

    // epilogue
    j = ((j1-2)-r);
    j++;
    // load 
    a00 = load1(a0, j-2); a10 = load1(a1, j-2); a20 = load1(a2, j-2); a30 = load1(a3, j-2); a40 = load1(a4, j-2);
    a01 = load1(a0, j-1); a11 = load1(a1, j-1); a21 = load1(a2, j-1); a31 = load1(a3, j-1); a41 = load1(a4, j-1);
    a02 = load1(a0, j  ); a12 = load1(a1, j  ); a22 = load1(a2, j  ); a32 = load1(a3, j  ); a42 = load1(a4, j  );
    a03 = load1(a0, j+1); a13 = load1(a1, j+1); a23 = load1(a2, j+1); a33 = load1(a3, j+1); a43 = load1(a4, j+1);
    // elu 2
    a50 = load1(a5, j-2);
    a51 = load1(a5, j-1);
    a52 = load1(a5, j  );
    a53 = load1(a5, j+1);

    // calcule des colonnes 
    c00 = MIN3( a00, a10, a20);
    c01 = MIN3( a01, a11, a21);
    c02 = MIN3( a02, a12, a22);
    c03 = MIN3( a03, a13, a23);

    c10 = MIN3( a10, a20, a30);
    c11 = MIN3( a11, a21, a31);
    c12 = MIN3( a12, a22, a32);
    c13 = MIN3( a13, a23, a33);

    c20 = MIN3(a20, a30, a40);
    c21 = MIN3(a21, a31, a41);
    c22 = MIN3(a22, a32, a42);
    c23 = MIN3(a23, a33, a43);

    //elu 2
    c30 = MIN3(a30, a40, a50);
    c31 = MIN3(a31, a41, a51);
    c32 = MIN3(a32, a42, a52);
    c33 = MIN3(a33, a43, a53);

    t00 = MIN3(LEFT_64(c00, c01),c01,RIGHT_64(c01, c02)); 
    t01 = MIN3(LEFT_64(c01, c02),c02,RIGHT_64(c02, c03));

    t10 = MIN3(LEFT_64(c10, c11),c11,RIGHT_64(c11, c12));
    t11 = MIN3(LEFT_64(c11, c12),c12,RIGHT_64(c12, c13));

    t20 = MIN3(LEFT_64(c20, c21),c21,RIGHT_64(c21, c22));
    t21 = MIN3(LEFT_64(c21, c22),c22,RIGHT_64(c22, c23));

    //elu 2
    t30 = MIN3(LEFT_64(c30, c31),c31,RIGHT_64(c31, c32));
    t31 = MIN3(LEFT_64(c31, c32),c32,RIGHT_64(c32, c33));

    switch (r+2)
    {
        case 4:
        // ----------------------------------LOOP 1-----------------------------------------------
        // load
        // ----------------------------------
        a04 = load1(a0, j+2); a14 = load1(a1, j+2); a24 = load1(a2, j+2); a34 = load1(a3, j+2); a44 = load1(a4, j+2);
        a54 = load1(a5, j+2);
        // ----------------------------------
        // calc
        // ----------------------------------
        c04 = MIN3( a04, a14, a24);
        c14 = MIN3( a14, a24, a34);
        c24 = MIN3(a24, a34, a44);
        c34 = MIN3( a34, a44, a54);
        //-----------------------------------
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        //elu 2
        t32 = MIN3(LEFT_64(c32, c33),c33,RIGHT_64(c33, c34));
        l3 = LEFT_64(t30,t31);
        r3 = RIGHT_64(t31,t32); 

        // calcule de la case finale 
        y0 = MAX9( 
                        l0, t01, r0,
                        l1, t11, r1,
                        l2, t21, r2
                    );
        
        y1 = MAX9( 
                        l1, t11, r1,
                        l2, t21, r2,
                        l3, t31, r3
                    );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi0, j, y0);
        store1(Yi1, j, y1);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;

        c32 = c33;
        c33 = c34;

        t30 = t31;
        t31 = t32;
        j++;

        case 3:
        // ----------------------------------LOOP 1-----------------------------------------------
        // load
        // ----------------------------------
        a04 = load1(a0, j+2); a14 = load1(a1, j+2); a24 = load1(a2, j+2); a34 = load1(a3, j+2); a44 = load1(a4, j+2);
        a54 = load1(a5, j+2);
        // ----------------------------------
        // calc
        // ----------------------------------
        c04 = MIN3( a04, a14, a24);
        c14 = MIN3( a14, a24, a34);
        c24 = MIN3(a24, a34, a44);
        c34 = MIN3( a34, a44, a54);
        //-----------------------------------
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        //elu 2
        t32 = MIN3(LEFT_64(c32, c33),c33,RIGHT_64(c33, c34));
        l3 = LEFT_64(t30,t31);
        r3 = RIGHT_64(t31,t32); 

        // calcule de la case finale 
        y0 = MAX9( 
                        l0, t01, r0,
                        l1, t11, r1,
                        l2, t21, r2
                    );
        
        y1 = MAX9( 
                        l1, t11, r1,
                        l2, t21, r2,
                        l3, t31, r3
                    );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi0, j, y0);
        store1(Yi1, j, y1);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;

        c32 = c33;
        c33 = c34;

        t30 = t31;
        t31 = t32;
        j++;

        case 2:
        // ----------------------------------LOOP 1-----------------------------------------------
        // load
        // ----------------------------------
        a04 = load1(a0, j+2); a14 = load1(a1, j+2); a24 = load1(a2, j+2); a34 = load1(a3, j+2); a44 = load1(a4, j+2);
        a54 = load1(a5, j+2);
        // ----------------------------------
        // calc
        // ----------------------------------
        c04 = MIN3( a04, a14, a24);
        c14 = MIN3( a14, a24, a34);
        c24 = MIN3(a24, a34, a44);
        c34 = MIN3( a34, a44, a54);
        //-----------------------------------
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        //elu 2
        t32 = MIN3(LEFT_64(c32, c33),c33,RIGHT_64(c33, c34));
        l3 = LEFT_64(t30,t31);
        r3 = RIGHT_64(t31,t32); 

        // calcule de la case finale 
        y0 = MAX9( 
                        l0, t01, r0,
                        l1, t11, r1,
                        l2, t21, r2
                    );
        
        y1 = MAX9( 
                        l1, t11, r1,
                        l2, t21, r2,
                        l3, t31, r3
                    );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi0, j, y0);
        store1(Yi1, j, y1);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;

        c32 = c33;
        c33 = c34;

        t30 = t31;
        t31 = t32;
        j++;

        case 1:
        // ----------------------------------LOOP 1-----------------------------------------------
        // load
        // ----------------------------------
        a04 = load1(a0, j+2); a14 = load1(a1, j+2); a24 = load1(a2, j+2); a34 = load1(a3, j+2); a44 = load1(a4, j+2);
        a54 = load1(a5, j+2);
        // ----------------------------------
        // calc
        // ----------------------------------
        c04 = MIN3( a04, a14, a24);
        c14 = MIN3( a14, a24, a34);
        c24 = MIN3(a24, a34, a44);
        c34 = MIN3( a34, a44, a54);
        //-----------------------------------
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        //elu 2
        t32 = MIN3(LEFT_64(c32, c33),c33,RIGHT_64(c33, c34));
        l3 = LEFT_64(t30,t31);
        r3 = RIGHT_64(t31,t32); 

        // calcule de la case finale 
        y0 = MAX9( 
                        l0, t01, r0,
                        l1, t11, r1,
                        l2, t21, r2
                    );
        
        y1 = MAX9( 
                        l1, t11, r1,
                        l2, t21, r2,
                        l3, t31, r3
                    );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi0, j, y0);
        store1(Yi1, j, y1);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;

        c32 = c33;
        c33 = c34;

        t30 = t31;
        t31 = t32;
        j++;

        default:
            break;
    }
}
// -------------------------------------------------------------------------------------------------------------
void line_ouverture3_ui64matrix_swp64_fusion_ilu3_elu2_red_factor(uint64 **X, int i, int j0, int j1, uint64 **Y)
// -------------------------------------------------------------------------------------------------------------
{
    /*
        Version avec deroulage de boucle interne de degre 3,
        deroulage de boucle interne de degre 2 et reduction de colonnes de fusion 
    */
   // Var des 5 lignes 
    uint64 *a0 = X[i-2]; // colonnes X
    uint64 *a1 = X[i-1];
    uint64 *a2 = X[i+0];
    uint64 *a3 = X[i+1]; 
    uint64 *a4 = X[i+2];
    uint64 *a5 = X[i+3]; //elu 2

    uint64 *Yi0 = Y[i];
    uint64 *Yi1 = Y[i+1];
    uint64 y0,y1;

    uint64 l0,r0,l1,r1,l2,r2,l3,r3;

    uint64 f; // factorisation

    uint64 f0, f1, f2, f3;

    uint64 t00, t01, t02; // colonnes temporaire
    uint64 t10, t11, t12;
    uint64 t20, t21, t22; 
    uint64 t30, t31, t32; // elu2

    // cases de X 
    uint64 a00, a01, a02, a03, a04;
    uint64 a10, a11, a12, a13, a14;
    uint64 a20, a21, a22, a23, a24;
    uint64 a30, a31, a32, a33, a34;
    uint64 a40, a41, a42, a43, a44;
    uint64 a50, a51, a52, a53, a54;

    uint64 c00, c01, c02, c03, c04;
    uint64 c10, c11, c12, c13, c14;
    uint64 c20, c21, c22, c23, c24;
    uint64 c30, c31, c32, c33, c34;

    int j = j0;

    // load 
    a00 = load1(a0, j-2); a10 = load1(a1, j-2); a20 = load1(a2, j-2); a30 = load1(a3, j-2); a40 = load1(a4, j-2);
    a01 = load1(a0, j-1); a11 = load1(a1, j-1); a21 = load1(a2, j-1); a31 = load1(a3, j-1); a41 = load1(a4, j-1);
    a02 = load1(a0, j  ); a12 = load1(a1, j  ); a22 = load1(a2, j  ); a32 = load1(a3, j  ); a42 = load1(a4, j  );
    a03 = load1(a0, j+1); a13 = load1(a1, j+1); a23 = load1(a2, j+1); a33 = load1(a3, j+1); a43 = load1(a4, j+1);
    // elu 2
    a50 = load1(a5, j-2);
    a51 = load1(a5, j-1);
    a52 = load1(a5, j  );
    a53 = load1(a5, j+1);

    // calcule des colonnes 
    f   = MIN(  a10, a20);
    c00 = MIN(  a00, f);
    c10 = MIN(  f, a30);

    f   = MIN(  a11, a21);
    c01 = MIN(  a01, f);
    c11 = MIN(  f, a31);

    f   = MIN(  a12, a22);
    c02 = MIN(  a02, f);
    c12 = MIN(  f, a32);
    
    f   = MIN(  a13, a23);
    c03 = MIN(  a03, f);
    c13 = MIN(  f, a33);

    f   = MIN(  a30, a40);
    c20 = MIN(  a20, f);
    c30 = MIN(  f, a50);

    f   = MIN(  a31, a41);
    c21 = MIN(  a21, f);
    c31 = MIN(  f, a51);

    f   = MIN(  a32, a42);
    c22 = MIN(  a22, f);
    c32 = MIN(  f, a52);

    f   = MIN(  a33, a43);
    c23 = MIN( a23, f);
    c33 = MIN( f, a53);

    
    t00 = MIN3(LEFT_64(c00, c01),c01,RIGHT_64(c01, c02)); 
    t01 = MIN3(LEFT_64(c01, c02),c02,RIGHT_64(c02, c03));

    t10 = MIN3(LEFT_64(c10, c11),c11,RIGHT_64(c11, c12));
    t11 = MIN3(LEFT_64(c11, c12),c12,RIGHT_64(c12, c13));

    t20 = MIN3(LEFT_64(c20, c21),c21,RIGHT_64(c21, c22));
    t21 = MIN3(LEFT_64(c21, c22),c22,RIGHT_64(c22, c23));

    //elu 2
    t30 = MIN3(LEFT_64(c30, c31),c31,RIGHT_64(c31, c32));
    t31 = MIN3(LEFT_64(c31, c32),c32,RIGHT_64(c32, c33));

    // determination de l'intervalle de j
    int r = ((j1 - 2) % 3);
    for (j = j0; j <= ((j1 - 2) - r); j += 3){
        // ----------------------------------LOOP 1-----------------------------------------------
        // load
        // ----------------------------------
        a04 = load1(a0, j+2); a14 = load1(a1, j+2); a24 = load1(a2, j+2); a34 = load1(a3, j+2); a44 = load1(a4, j+2);
        a54 = load1(a5, j+2);
        // ----------------------------------
        // calc
        // ----------------------------------
        f   = MIN(a14, a24);
        c04 = MIN( a04, f);
        c14 = MIN( f, a34);

        f   = MIN(a34, a44);
        c24 = MIN(a24, f);
        c34 = MIN( f, a54);
        //-----------------------------------
        
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        //elu 2
        t32 = MIN3(LEFT_64(c32, c33),c33,RIGHT_64(c33, c34));
        l3 = LEFT_64(t30,t31);
        r3 = RIGHT_64(t31,t32); 

        // calcule de la case finale
        f = MAX6(l1, t11, r1,
                l2, t21, r2);

        y0 = MAX4( l0, t01, r0, f );
        
        y1 = MAX4( f, l3, t31, r3 );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi0, j, y0);
        store1(Yi1, j, y1);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;

        c32 = c33;
        c33 = c34;

        t30 = t31;
        t31 = t32;


        // ----------------------------------LOOP 2-----------------------------------------------
        // ----------------------------------
        // load
        // ----------------------------------
        a04 = load1(a0, j+3); a14 = load1(a1, j+3); a24 = load1(a2, j+3); a34 = load1(a3, j+3); a44 = load1(a4, j+3);
        a54 = load1(a5, j+3);
        // ----------------------------------
        // calc
        // ----------------------------------
        f = MIN(a14, a24);
        c04 = MIN( a04, f);
        c14 = MIN( f, a34);

        f = MIN(a34, a44);
        c24 = MIN(a24, f);
        c34 = MIN( f, a54);
        //-----------------------------------
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        //elu 2
        t32 = MIN3(LEFT_64(c32, c33),c33,RIGHT_64(c33, c34));
        l3 = LEFT_64(t30,t31);
        r3 = RIGHT_64(t31,t32); 

        // calcule de la case finale
        f = MAX6(l1, t11, r1,
                l2, t21, r2);

        y0 = MAX4( l0, t01, r0, f );
        
        y1 = MAX4( f, l3, t31, r3 );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi0, j+1, y0);
        store1(Yi1, j+1, y1);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;

        c32 = c33;
        c33 = c34;

        t30 = t31;
        t31 = t32;

        // ----------------------------------LOOP 3-----------------------------------------------
        // ----------------------------------
        // load
        // ----------------------------------
        a04 = load1(a0, j+4); a14 = load1(a1, j+4); a24 = load1(a2, j+4); a34 = load1(a3, j+4); a44 = load1(a4, j+4);
        a54 = load1(a5, j+4);
        // ----------------------------------
        // calc
        // ----------------------------------
        f = MIN(a14, a24);
        c04 = MIN( a04, f);
        c14 = MIN( f, a34);

        f = MIN(a34, a44);
        c24 = MIN(a24, f);
        c34 = MIN( f, a54);
        //-----------------------------------
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        //elu 2
        t32 = MIN3(LEFT_64(c32, c33),c33,RIGHT_64(c33, c34));
        l3 = LEFT_64(t30,t31);
        r3 = RIGHT_64(t31,t32); 

        // calcule de la case finale
        f = MAX6(l1, t11, r1,
                l2, t21, r2);

        y0 = MAX4( l0, t01, r0, f );
        
        y1 = MAX4( f, l3, t31, r3 );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi0, j+2, y0);
        store1(Yi1, j+2, y1);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;

        c32 = c33;
        c33 = c34;

        t30 = t31;
        t31 = t32;
    }

    // epilogue
    if(r+2 == 0) return; // ne fait pas l'épilogue si inutile 
    j = ((j1-2)-r);
    j++;
    // load 
    a00 = load1(a0, j-2); a10 = load1(a1, j-2); a20 = load1(a2, j-2); a30 = load1(a3, j-2); a40 = load1(a4, j-2);
    a01 = load1(a0, j-1); a11 = load1(a1, j-1); a21 = load1(a2, j-1); a31 = load1(a3, j-1); a41 = load1(a4, j-1);
    a02 = load1(a0, j  ); a12 = load1(a1, j  ); a22 = load1(a2, j  ); a32 = load1(a3, j  ); a42 = load1(a4, j  );
    a03 = load1(a0, j+1); a13 = load1(a1, j+1); a23 = load1(a2, j+1); a33 = load1(a3, j+1); a43 = load1(a4, j+1);
    // elu 2
    a50 = load1(a5, j-2);
    a51 = load1(a5, j-1);
    a52 = load1(a5, j  );
    a53 = load1(a5, j+1);

    // calcule des colonnes 
    c00 = MIN3( a00, a10, a20); 
    c01 = MIN3( a01, a11, a21);
    c02 = MIN3( a02, a12, a22);
    c03 = MIN3( a03, a13, a23);

    c10 = MIN3( a10, a20, a30);
    c11 = MIN3( a11, a21, a31);
    c12 = MIN3( a12, a22, a32);
    c13 = MIN3( a13, a23, a33);

    c20 = MIN3(a20, a30, a40);
    c21 = MIN3(a21, a31, a41);
    c22 = MIN3(a22, a32, a42);
    c23 = MIN3(a23, a33, a43);

    //elu 2
    c30 = MIN3(a30, a40, a50);
    c31 = MIN3(a31, a41, a51);
    c32 = MIN3(a32, a42, a52);
    c33 = MIN3(a33, a43, a53);

    f   = MIN(  a10, a20);
    c00 = MIN(  a00, f);
    c10 = MIN(  f, a30);

    f   = MIN(  a11, a21);
    c01 = MIN(  a01, f);
    c11 = MIN(  f, a31);

    f   = MIN(  a12, a22);
    c02 = MIN(  a02, f);
    c12 = MIN(  f, a32);
    
    f   = MIN(  a13, a23);
    c03 = MIN(  a03, f);
    c13 = MIN(  f, a33);

    f   = MIN(  a30, a40);
    c20 = MIN(  a20, f);
    c30 = MIN(  f, a50);

    f   = MIN(  a31, a41);
    c21 = MIN(  a21, f);
    c31 = MIN(  f, a51);

    f   = MIN(  a32, a42);
    c22 = MIN(  a22, f);
    c32 = MIN(  f, a52);

    f   = MIN(  a33, a43);
    c23 = MIN( a23, f);
    c33 = MIN( f, a53);

    t00 = MIN3(LEFT_64(c00, c01),c01,RIGHT_64(c01, c02)); 
    t01 = MIN3(LEFT_64(c01, c02),c02,RIGHT_64(c02, c03));

    t10 = MIN3(LEFT_64(c10, c11),c11,RIGHT_64(c11, c12));
    t11 = MIN3(LEFT_64(c11, c12),c12,RIGHT_64(c12, c13));

    t20 = MIN3(LEFT_64(c20, c21),c21,RIGHT_64(c21, c22));
    t21 = MIN3(LEFT_64(c21, c22),c22,RIGHT_64(c22, c23));

    //elu 2
    t30 = MIN3(LEFT_64(c30, c31),c31,RIGHT_64(c31, c32));
    t31 = MIN3(LEFT_64(c31, c32),c32,RIGHT_64(c32, c33));

    switch (r+2)
    {
        case 4:
        // load
        // ----------------------------------
        a04 = load1(a0, j+2); a14 = load1(a1, j+2); a24 = load1(a2, j+2); a34 = load1(a3, j+2); a44 = load1(a4, j+2);
        a54 = load1(a5, j+2);
        // ----------------------------------
        // calc
        // ----------------------------------
        f   = MIN(a14, a24);
        c04 = MIN( a04, f);
        c14 = MIN( f, a34);

        f   = MIN(a34, a44);
        c24 = MIN(a24, f);
        c34 = MIN( f, a54);
        //-----------------------------------
        
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        //elu 2
        t32 = MIN3(LEFT_64(c32, c33),c33,RIGHT_64(c33, c34));
        l3 = LEFT_64(t30,t31);
        r3 = RIGHT_64(t31,t32); 

        // calcule de la case finale
        f = MAX6(l1, t11, r1,
                l2, t21, r2);

        y0 = MAX4( l0, t01, r0, f );
        
        y1 = MAX4( f, l3, t31, r3 );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi0, j, y0);
        store1(Yi1, j, y1);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;

        c32 = c33;
        c33 = c34;

        t30 = t31;
        t31 = t32;
        j++;

        case 3:
        // load
        // ----------------------------------
        a04 = load1(a0, j+2); a14 = load1(a1, j+2); a24 = load1(a2, j+2); a34 = load1(a3, j+2); a44 = load1(a4, j+2);
        a54 = load1(a5, j+2);
        // ----------------------------------
        // calc
        // ----------------------------------
        f   = MIN(a14, a24);
        c04 = MIN( a04, f);
        c14 = MIN( f, a34);

        f   = MIN(a34, a44);
        c24 = MIN(a24, f);
        c34 = MIN( f, a54);
        //-----------------------------------
        
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        //elu 2
        t32 = MIN3(LEFT_64(c32, c33),c33,RIGHT_64(c33, c34));
        l3 = LEFT_64(t30,t31);
        r3 = RIGHT_64(t31,t32); 

        // calcule de la case finale
        f = MAX6(l1, t11, r1,
                l2, t21, r2);

        y0 = MAX4( l0, t01, r0, f );
        
        y1 = MAX4( f, l3, t31, r3 );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi0, j, y0);
        store1(Yi1, j, y1);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;

        c32 = c33;
        c33 = c34;

        t30 = t31;
        t31 = t32;
        j++;

        case 2:
        // load
        // ----------------------------------
        a04 = load1(a0, j+2); a14 = load1(a1, j+2); a24 = load1(a2, j+2); a34 = load1(a3, j+2); a44 = load1(a4, j+2);
        a54 = load1(a5, j+2);
        // ----------------------------------
        // calc
        // ----------------------------------
        f   = MIN(a14, a24);
        c04 = MIN( a04, f);
        c14 = MIN( f, a34);

        f   = MIN(a34, a44);
        c24 = MIN(a24, f);
        c34 = MIN( f, a54);
        //-----------------------------------
        
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        //elu 2
        t32 = MIN3(LEFT_64(c32, c33),c33,RIGHT_64(c33, c34));
        l3 = LEFT_64(t30,t31);
        r3 = RIGHT_64(t31,t32); 

        // calcule de la case finale
        f = MAX6(l1, t11, r1,
                l2, t21, r2);

        y0 = MAX4( l0, t01, r0, f );
        
        y1 = MAX4( f, l3, t31, r3 );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi0, j, y0);
        store1(Yi1, j, y1);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;

        c32 = c33;
        c33 = c34;

        t30 = t31;
        t31 = t32;
        j++;

        case 1:
        // load
        // ----------------------------------
        a04 = load1(a0, j+2); a14 = load1(a1, j+2); a24 = load1(a2, j+2); a34 = load1(a3, j+2); a44 = load1(a4, j+2);
        a54 = load1(a5, j+2);
        // ----------------------------------
        // calc
        // ----------------------------------
        f   = MIN(a14, a24);
        c04 = MIN( a04, f);
        c14 = MIN( f, a34);

        f   = MIN(a34, a44);
        c24 = MIN(a24, f);
        c34 = MIN( f, a54);
        //-----------------------------------
        
        t02 = MIN3(LEFT_64(c02, c03),c03,RIGHT_64(c03, c04));
        l0 = LEFT_64(t00,t01);
        r0 = RIGHT_64(t01,t02);
        // ----------------------------------
        t12 = MIN3(LEFT_64(c12, c13),c13,RIGHT_64(c13, c14));
        l1 = LEFT_64(t10,t11);
        r1 = RIGHT_64(t11,t12);
        // ----------------------------------
        t22 = MIN3(LEFT_64(c22, c23),c23,RIGHT_64(c23, c24));
        l2 = LEFT_64(t20,t21);
        r2 = RIGHT_64(t21,t22); 
        // ----------------------------------
        //elu 2
        t32 = MIN3(LEFT_64(c32, c33),c33,RIGHT_64(c33, c34));
        l3 = LEFT_64(t30,t31);
        r3 = RIGHT_64(t31,t32); 

        // calcule de la case finale
        f = MAX6(l1, t11, r1,
                l2, t21, r2);

        y0 = MAX4( l0, t01, r0, f );
        
        y1 = MAX4( f, l3, t31, r3 );
        
        // ----------------------------------
        // store
        // ----------------------------------
        store1(Yi0, j, y0);
        store1(Yi1, j, y1);
        
        // ----------------------------------
        // Rotation de colonne et des var intermediaire
        // ----------------------------------
        c02 = c03;
        c03 = c04;

        t00 = t01;
        t01 = t02;

        c12 = c13;
        c13 = c14;

        t10 = t11;
        t11 = t12;

        c22 = c23;
        c23 = c24;

        t20 = t21;
        t21 = t22;

        c32 = c33;
        c33 = c34;

        t30 = t31;
        t31 = t32;
        j++;

        default:
            break;
    }
}
// -------------------------------------------------------------------------------------------------------
void ouverture3_ui64matrix_swp64_basic(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y, uint64 **Z)
// -------------------------------------------------------------------------------------------------------
{
    erosion3_ui64matrix_swp64_basic(X, i0-1, i1+1, j0-1, j1+1, Y);
    dilatation3_ui64matrix_swp64_basic(Y, i0,   i1,   j0,   j1,   Z);
}
// --------------------------------------------------------------------------------------------------
void ouverture3_ui64matrix_swp64_fusion_basic(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// --------------------------------------------------------------------------------------------------
{
    int i;
    for(i = i0; i <=i1; i++)
    {
        line_ouverture3_ui64matrix_swp64_fusion_basic(X, i, j0, j1, Y);
    }
}
// ------------------------------------------------------------------------------------------------
void ouverture3_ui64matrix_swp64_fusion_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// ------------------------------------------------------------------------------------------------
{
    int i;
    for(i = i0; i <=i1; i++)
    {
        line_ouverture3_ui64matrix_swp64_fusion_red(X, i, j0, j1, Y);
    }
}
// -----------------------------------------------------------------------------------------------------
void ouverture3_ui64matrix_swp64_fusion_ilu3_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// -----------------------------------------------------------------------------------------------------
{
    int i;
    for(i = i0; i <=i1; i++)
    {
        line_ouverture3_ui64matrix_swp64_fusion_ilu3_red(X, i, j0, j1, Y);
    }
}
// -----------------------------------------------------------------------------------------------------
void ouverture3_ui64matrix_swp64_fusion_elu2_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// -----------------------------------------------------------------------------------------------------
{
    int i;
    int r = ((i1-1)%2);
    for(i = i0; i <=((i1-1)-r); i+=2)
    {
        line_ouverture3_ui64matrix_swp64_fusion_elu2_red(X, i, j0, j1, Y);
    }
    i = ((i1-1)-r);
    switch (r+1)
    {
    case 2:
        i++;
        line_ouverture3_ui64matrix_swp64_fusion_red(X, i, j0, j1, Y);
    case 1:
        i++;
        line_ouverture3_ui64matrix_swp64_fusion_red(X, i, j0, j1, Y);
    default:
        break;
    }
}
// ----------------------------------------------------------------------------------------------------------
void ouverture3_ui64matrix_swp64_fusion_ilu3_elu2_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// ----------------------------------------------------------------------------------------------------------
{
    int i;
    int r = ((i1-1)%2);
    for(i = i0; i <=((i1-1)-r); i+=2)
    {
        line_ouverture3_ui64matrix_swp64_fusion_ilu3_elu2_red(X, i, j0, j1, Y);
    }
    i = ((i1-1)-r);
    switch (r+1)
    {
    case 2:
        i++;
        line_ouverture3_ui64matrix_swp64_fusion_ilu3_red(X, i, j0, j1, Y);
    case 1:
        i++;
        line_ouverture3_ui64matrix_swp64_fusion_ilu3_red(X, i, j0, j1, Y);
    default:
        break;
    }
}
// -----------------------------------------------------------------------------------------------------------------
void ouverture3_ui64matrix_swp64_fusion_ilu3_elu2_red_factor(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// -----------------------------------------------------------------------------------------------------------------
{
    int i;
    int r = ((i1-1)%2);
    for(i = i0; i <=((i1-1)-r); i+=2)
    {
        line_ouverture3_ui64matrix_swp64_fusion_ilu3_elu2_red_factor(X, i, j0, j1, Y);
    }
    i = ((i1-1)-r);
    switch (r+1)
    {
    case 2:
        i++;
        line_ouverture3_ui64matrix_swp64_fusion_ilu3_red(X, i, j0, j1, Y);
    case 1:
        i++;
        line_ouverture3_ui64matrix_swp64_fusion_ilu3_red(X, i, j0, j1, Y);
    default:
        break;
    }
}
// ----------------------------------------------------------------------------------------------------------------
void ouverture3_ui64matrix_swp64_pipeline_basic(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y, uint64 **Z)
// ----------------------------------------------------------------------------------------------------------------
{
    int i;
    int b = 1;
    i =i0-b;
    line_erosion3_ui64matrix_swp64_basic(X, i  , j0-b, j1+b, Y);
    line_erosion3_ui64matrix_swp64_basic(X, i+1, j0-b, j1+b, Y);
    for(i=i0; i<=i1; i++)
    {
        line_erosion3_ui64matrix_swp64_basic(X, i+1, j0-b  , j1+b  , Y);
        line_dilatation3_ui64matrix_swp64_basic(Y, i  , j0-b+1, j1+b-1, Z);
    }
}
// --------------------------------------------------------------------------------------------------------------
void ouverture3_ui64matrix_swp64_pipeline_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y, uint64 **Z)
// --------------------------------------------------------------------------------------------------------------
{
    int i;
    int b = 1;
    i =i0-b;
    line_erosion3_ui64matrix_swp64_red(X, i  , j0-b, j1+b, Y);
    line_erosion3_ui64matrix_swp64_red(X, i+1, j0-b, j1+b, Y);
    for(i=i0; i<=i1; i++)
    {
        line_erosion3_ui64matrix_swp64_red(X, i+1, j0-b  , j1+b  , Y);
        line_dilatation3_ui64matrix_swp64_red(Y, i  , j0-b+1, j1+b-1, Z);
    }
}
// -------------------------------------------------------------------------------------------------------------------
void ouverture3_ui64matrix_swp64_pipeline_ilu3_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y, uint64 **Z)
// -------------------------------------------------------------------------------------------------------------------
{
    int i;
    int b = 1;
    i =i0-b;
    line_erosion3_ui64matrix_swp64_ilu3_red(X, i  , j0-b, j1+b, Y);
    line_erosion3_ui64matrix_swp64_ilu3_red(X, i+1, j0-b, j1+b, Y);
    for(i=i0; i<=i1; i++)
    {
        line_erosion3_ui64matrix_swp64_ilu3_red(X, i+1, j0-b  , j1+b  , Y);
        line_dilatation3_ui64matrix_swp64_ilu3_red(Y, i  , j0-b+1, j1+b-1, Z);
    }
}
// -------------------------------------------------------------------------------------------------------------------
void ouverture3_ui64matrix_swp64_pipeline_elu2_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y, uint64 **Z)
// -------------------------------------------------------------------------------------------------------------------
{
    int i;
    int b = 1;
    i =i0-b;
    line_erosion3_ui64matrix_swp64_red(X, i  , j0-b, j1+b, Y);
    line_erosion3_ui64matrix_swp64_red(X, i+1, j0-b, j1+b, Y);
    
    int r = ((i1-1)%2);

    for(i=i0; i<=((i1-1)-r); i+=2)
    {
        line_erosion3_ui64matrix_swp64_elu2_red(X, i+1, j0-b  , j1+b  , Y);
        line_dilatation3_ui64matrix_swp64_elu2_red(Y, i  , j0-b+1, j1+b-1, Z);
    }

    i = ((i1-1)-r);
    switch(r+1)
    {
    case 2:
        i++;
        line_erosion3_ui64matrix_swp64_red(X, i+1, j0-b  , j1+b  , Y);
        line_dilatation3_ui64matrix_swp64_red(Y, i  , j0-b+1, j1+b-1, Z);
    case 1:
        i++; 
        line_erosion3_ui64matrix_swp64_red(X, i+1, j0-b  , j1+b  , Y);
        line_dilatation3_ui64matrix_swp64_red(Y, i  , j0-b+1, j1+b-1, Z);
    default:
        break;
    }
}
// --------------------------------------------------------------------------------------------------------------------------
void ouverture3_ui64matrix_swp64_pipeline_elu2_red_factor(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y, uint64 **Z)
// --------------------------------------------------------------------------------------------------------------------------
{
    int i;
    int b = 1;
    i =i0-b;
    line_erosion3_ui64matrix_swp64_red(X, i  , j0-b, j1+b, Y);
    line_erosion3_ui64matrix_swp64_red(X, i+1, j0-b, j1+b, Y);
    
    int r = ((i1-1)%2);

    for(i=i0; i<=((i1-1)-r); i+=2)
    {
        line_erosion3_ui64matrix_swp64_elu2_red_factor(X, i+1, j0-b  , j1+b  , Y);
        line_dilatation3_ui64matrix_swp64_elu2_red_factor(Y, i  , j0-b+1, j1+b-1, Z);
    }

    i = ((i1-1)-r);
    switch(r+1)
    {
    case 2:
        i++;
        line_erosion3_ui64matrix_swp64_red(X, i+1, j0-b  , j1+b  , Y);
        line_dilatation3_ui64matrix_swp64_red(Y, i  , j0-b+1, j1+b-1, Z);
    case 1:
        i++; 
        line_erosion3_ui64matrix_swp64_red(X, i+1, j0-b  , j1+b  , Y);
        line_dilatation3_ui64matrix_swp64_red(Y, i  , j0-b+1, j1+b-1, Z);
    default:
        break;
    }
}
// ------------------------------------------------------------------------------------------------------------------------
void ouverture3_ui64matrix_swp64_pipeline_ilu3_elu2_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y, uint64 **Z)
// ------------------------------------------------------------------------------------------------------------------------
{
    int i;
    int b = 1;
    i =i0-b;
    line_erosion3_ui64matrix_swp64_ilu3_red(X, i  , j0-b, j1+b, Y);
    line_erosion3_ui64matrix_swp64_ilu3_red(X, i+1, j0-b, j1+b, Y);
    
    int r = ((i1-1)%2);

    for(i=i0; i<=((i1-1)-r); i+=2)
    {
        line_erosion3_ui64matrix_swp64_ilu3_elu2_red(X, i+1, j0-b  , j1+b  , Y);
        line_dilatation3_ui64matrix_swp64_ilu3_elu2_red(Y, i  , j0-b+1, j1+b-1, Z);
    }

    i = ((i1-1)-r);
    switch(r+1)
    {
    case 2:
        i++;
        line_erosion3_ui64matrix_swp64_ilu3_red(X, i+1, j0-b  , j1+b  , Y);
        line_dilatation3_ui64matrix_swp64_ilu3_red(Y, i  , j0-b+1, j1+b-1, Z);
    case 1:
        i++; 
        line_erosion3_ui64matrix_swp64_ilu3_red(X, i+1, j0-b  , j1+b  , Y);
        line_dilatation3_ui64matrix_swp64_ilu3_red(Y, i  , j0-b+1, j1+b-1, Z);
    default:
        break;
    }
}
// -------------------------------------------------------------------------------------------------------------------------------
void ouverture3_ui64matrix_swp64_pipeline_ilu3_elu2_red_factor(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y, uint64 **Z)
// -------------------------------------------------------------------------------------------------------------------------------
{
    int i;
    int b = 1;
    i =i0-b;
    line_erosion3_ui64matrix_swp64_ilu3_red(X, i  , j0-b, j1+b, Y);
    line_erosion3_ui64matrix_swp64_ilu3_red(X, i+1, j0-b, j1+b, Y);
    
    int r = ((i1-1)%2);

    for(i=i0; i<=((i1-1)-r); i+=2)
    {
        line_erosion3_ui64matrix_swp64_ilu3_elu2_red_factor(X, i+1, j0-b  , j1+b  , Y);
        line_dilatation3_ui64matrix_swp64_ilu3_elu2_red_factor(Y, i  , j0-b+1, j1+b-1, Z);
    }

    i = ((i1-1)-r);
    switch(r+1)
    {
    case 2:
        i++;
        line_erosion3_ui64matrix_swp64_ilu3_red(X, i+1, j0-b  , j1+b  , Y);
        line_dilatation3_ui64matrix_swp64_ilu3_red(Y, i  , j0-b+1, j1+b-1, Z);
    case 1:
        i++; 
        line_erosion3_ui64matrix_swp64_ilu3_red(X, i+1, j0-b  , j1+b  , Y);
        line_dilatation3_ui64matrix_swp64_ilu3_red(Y, i  , j0-b+1, j1+b-1, Z);
    default:
        break;
    }
}
