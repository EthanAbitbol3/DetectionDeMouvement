/* ------------------------------ */
/* --- morpho_erosion_swp64.c --- */
/* ------------------------------ */

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
//#include "sequence.h"

#include "swp.h"
#include "morpho.h"

// -------------------------------------------------------------------------------------
void line_erosion3_ui64matrix_swp64_basic(uint64 **X, int i, int j0, int j1, uint64 **Y)
// -------------------------------------------------------------------------------------
{
    int j;
    for (j = j0; j <= j1; j++){
        Y[i][j]=MIN9(
                    LEFT_64(X[i-1][j-1], X[i-1][j]), X[i-1][j], RIGHT_64(X[i-1][j],X[i-1][j+1]),
                    LEFT_64(X[i  ][j-1], X[i  ][j]), X[i  ][j], RIGHT_64(X[i  ][j],X[i  ][j+1]),
                    LEFT_64(X[i+1][j-1], X[i+1][j]), X[i+1][j], RIGHT_64(X[i+1][j],X[i+1][j+1])
                );
    }
}
// -----------------------------------------------------------------------------------
void line_erosion3_ui64matrix_swp64_reg(uint64 **X, int i, int j0, int j1, uint64 **Y)
// -----------------------------------------------------------------------------------
{
    // debug
    setborder1(j0, j1);

    uint64 *a0, *a1, *a2;
    a0 = X[i-1];
    a1 = X[i  ];
    a2 = X[i+1];

    uint64 y;
    uint64 *y0 = Y[i];
    
    uint64 l0,b0,r0,l1,b1,r1,l2,b2,r2; 
    
    uint64 a00, a01, a02;
    uint64 a10, a11, a12;
    uint64 a20, a21, a22;

    int j;
    for (j = j0; j <= j1; j++){
        
        // load
        a00 = load1(a0, j - 1); a10 = load1(a1, j - 1); a20 = load1(a2, j - 1);
        a01 = load1(a0, j + 0); a11 = load1(a1, j + 0); a21 = load1(a2, j + 0);
        a02 = load1(a0, j + 1); a12 = load1(a1, j + 1); a22 = load1(a2, j + 1);

        l0 = LEFT_64(a00, a01);
        b0 = a01;
        r0 = RIGHT_64(a01, a02);

        l1 = LEFT_64(a10, a11);
        b1 = a11;
        r1 = RIGHT_64(a11, a12);
        
        l2 = LEFT_64(a20, a21);
        b2 = a21;
        r2 = RIGHT_64(a21,a22);


        y = MIN9(
                l0, b0, r0,
                l1, b1, r1,
                l2, b2, r2
            );

        store1(y0, j, y);

    }
}
// -----------------------------------------------------------------------------------
void line_erosion3_ui64matrix_swp64_rot(uint64 **X, int i, int j0, int j1, uint64 **Y)
// -----------------------------------------------------------------------------------
{
    // debug
    setborder1(j0, j1);

    uint64 *a0, *a1, *a2;
    a0 = X[i-1];
    a1 = X[i  ];
    a2 = X[i+1];

    uint64 y;
    uint64 *y0 = Y[i];
    
    uint64 l0,b0,r0,l1,b1,r1,l2,b2,r2; 
    
    uint64 a00, a01, a02;
    uint64 a10, a11, a12;
    uint64 a20, a21, a22;

    // prologue 
    int j = j0 ;
    a00 = load1(a0, j - 1); a01 = load1(a0, j + 0);  
    a10 = load1(a1, j - 1); a11 = load1(a1, j + 0); 
    a20 = load1(a2, j - 1); a21 = load1(a2, j + 0);

    for (j = j0; j <= j1; j++){
        
        // load
        a02 = load1(a0, j + 1); a12 = load1(a1, j + 1); a22 = load1(a2, j + 1);

        l0 = LEFT_64(a00, a01);
        b0 = a01;
        r0 = RIGHT_64(a01, a02);

        l1 = LEFT_64(a10, a11);
        b1 = a11;
        r1 = RIGHT_64(a11, a12);
        
        l2 = LEFT_64(a20, a21);
        b2 = a21;
        r2 = RIGHT_64(a21,a22);


        y = MIN9(
                l0, b0, r0,
                l1, b1, r1,
                l2, b2, r2
            );

        // store
        store1(y0, j, y);

        // reg rot
        a00 = a01; a10 = a11; a20 = a21;
        a01 = a02; a11 = a12; a21 = a22;
    }
}
// -----------------------------------------------------------------------------------
void line_erosion3_ui64matrix_swp64_red(uint64 **X, int i, int j0, int j1, uint64 **Y)
// -----------------------------------------------------------------------------------
{
    // debug
    setborder1(j0-1, j1+1);

    uint64 *a0, *a1, *a2;
    a0 = X[i-1];
    a1 = X[i  ];
    a2 = X[i+1];

    uint64 y;
    uint64 *y0 = Y[i];
    
    uint64 l0,r0; 
    
    uint64 a00, a01, a02;
    uint64 a10, a11, a12;
    uint64 a20, a21, a22;

    uint64 c0 , c1 , c2; 
    int j = j0 ;
    a00 = load1(a0, j - 1); a01 = load1(a0, j + 0);  
    a10 = load1(a1, j - 1); a11 = load1(a1, j + 0); 
    a20 = load1(a2, j - 1); a21 = load1(a2, j + 0);
    
    c0 = MIN3(a00, a10, a20); c1 = MIN3(a01, a11, a21);

    for (j = j0; j <= j1; j++){
        
        // load
        a02 = load1(a0, j + 1); 
        a12 = load1(a1, j + 1); 
        a22 = load1(a2, j + 1);
        
        c2 = MIN3(a02, a12, a22);

        // proc
        l0 = LEFT_64(c0,c1);
        r0 = RIGHT_64 (c1, c2);

        y = MIN3(l0,c1,r0);

        // store
        store1(y0, j, y);

        // rotation de colonne
        c0 = c1;
        c1 = c2;
    }
}
// ------------------------------------------------------------------------------------
void line_erosion3_ui64matrix_swp64_ilu3(uint64 **X, int i, int j0, int j1, uint64 **Y)
// ------------------------------------------------------------------------------------
{
    // debug
    setborder1(j0-1, j1+1);

    uint64 *a0, *a1, *a2;
    a0 = X[i-1];
    a1 = X[i  ];
    a2 = X[i+1];
    uint64 y;
    uint64 *y0 = Y[i];

    uint64 a00, a01, a02;
    uint64 a10, a11, a12;
    uint64 a20, a21, a22;

    uint64 l0,b0,r0,l1,b1,r1,l2,b2,r2;

// determination de l'intervalle de j
    int r = ((j1 - 2) % 3);
    int j; 
    for (j = j0; j <= ((j1 - 2) - r); j += 3){
        // loop 1 -------------------------------------------------------------------------------
        // load
        a00 = load1(a0, j - 1);  a01 = load1(a0, j + 0); a02 = load1(a0, j + 1);
        a10 = load1(a1, j - 1);  a11 = load1(a1, j + 0); a12 = load1(a1, j + 1);
        a20 = load1(a2, j - 1);  a21 = load1(a2, j + 0); a22 = load1(a2, j + 1);

        l0 = LEFT_64(a00, a01); b0 = a01;r0 = RIGHT_64(a01, a02);

        l1 = LEFT_64(a10, a11); b1 = a11;r1 = RIGHT_64(a11, a12);
        
        l2 = LEFT_64(a20, a21); b2 = a21;r2 = RIGHT_64(a21,a22);

        y = MIN9(
            l0, b0, r0,
            l1, b1, r1,
            l2, b2, r2
            );

        store1(y0, j, y);

        // loop 2 -------------------------------------------------------------------------------
        // load
        a00 = load1(a0, j + 0);  a01 = load1(a0, j + 1); a02 = load1(a0, j + 2);
        a10 = load1(a1, j + 0);  a11 = load1(a1, j + 1); a12 = load1(a1, j + 2);
        a20 = load1(a2, j + 0);  a21 = load1(a2, j + 1); a22 = load1(a2, j + 2);

        l0 = LEFT_64(a00, a01); b0 = a01;r0 = RIGHT_64(a01, a02);

        l1 = LEFT_64(a10, a11);b1 = a11;r1 = RIGHT_64(a11, a12);
        
        l2 = LEFT_64(a20, a21);b2 = a21;r2 = RIGHT_64(a21,a22);

        y =MIN9(
            l0, b0, r0,
            l1, b1, r1,
            l2, b2, r2
            );

        store1(y0, j + 1, y);

        // loop 3 -------------------------------------------------------------------------------
        // load
        a00 = load1(a0, j + 1);  a01 = load1(a0, j + 2); a02 = load1(a0, j + 3);
        a10 = load1(a1, j + 1);  a11 = load1(a1, j + 2); a12 = load1(a1, j + 3);
        a20 = load1(a2, j + 1);  a21 = load1(a2, j + 2); a22 = load1(a2, j + 3);

        l0 = LEFT_64(a00, a01); b0 = a01;r0 = RIGHT_64(a01, a02);

        l1 = LEFT_64(a10, a11);b1 = a11;r1 = RIGHT_64(a11, a12);
        
        l2 = LEFT_64(a20, a21);b2 = a21;r2 = RIGHT_64(a21,a22);

        y =MIN9(
            l0, b0, r0,
            l1, b1, r1,
            l2, b2, r2
            );

        store1(y0, j+2, y);
    }
        // epilogue
    j = ((j1 - 2) - r);
    switch (r + 2)
    {
    case 4:
        j++;
        // load
        a00 = load1(a0, j - 1);  a01 = load1(a0, j + 0); a02 = load1(a0, j + 1);
        a10 = load1(a1, j - 1);  a11 = load1(a1, j + 0); a12 = load1(a1, j + 1);
        a20 = load1(a2, j - 1);  a21 = load1(a2, j + 0); a22 = load1(a2, j + 1);

        l0 = LEFT_64(a00, a01); b0 = a01;r0 = RIGHT_64(a01, a02);

        l1 = LEFT_64(a10, a11);b1 = a11;r1 = RIGHT_64(a11, a12);
        
        l2 = LEFT_64(a20, a21);b2 = a21;r2 = RIGHT_64(a21,a22);

        y =MIN9(
            l0, b0, r0,
            l1, b1, r1,
            l2, b2, r2
            );

        store1(y0, j, y);
        
    case 3:
        j++;
        // load
        a00 = load1(a0, j - 1);  a01 = load1(a0, j + 0); a02 = load1(a0, j + 1);
        a10 = load1(a1, j - 1);  a11 = load1(a1, j + 0); a12 = load1(a1, j + 1);
        a20 = load1(a2, j - 1);  a21 = load1(a2, j + 0); a22 = load1(a2, j + 1);

        l0 = LEFT_64(a00, a01); b0 = a01;r0 = RIGHT_64(a01, a02);

        l1 = LEFT_64(a10, a11);b1 = a11;r1 = RIGHT_64(a11, a12);
        
        l2 = LEFT_64(a20, a21);b2 = a21;r2 = RIGHT_64(a21,a22);
        y =MIN9(
            l0, b0, r0,
            l1, b1, r1,
            l2, b2, r2
            );

        store1(y0, j, y);

    case 2:
        j++;
        // load
        a00 = load1(a0, j - 1);  a01 = load1(a0, j + 0); a02 = load1(a0, j + 1);
        a10 = load1(a1, j - 1);  a11 = load1(a1, j + 0); a12 = load1(a1, j + 1);
        a20 = load1(a2, j - 1);  a21 = load1(a2, j + 0); a22 = load1(a2, j + 1);

        l0 = LEFT_64(a00, a01); b0 = a01;r0 = RIGHT_64(a01, a02);

        l1 = LEFT_64(a10, a11);b1 = a11;r1 = RIGHT_64(a11, a12);
        
        l2 = LEFT_64(a20, a21);b2 = a21;r2 = RIGHT_64(a21,a22);

        y =MIN9(
            l0, b0, r0,
            l1, b1, r1,
            l2, b2, r2
            );

        store1(y0, j, y);

    case 1:
        j++;
        // load
        a00 = load1(a0, j - 1);  a01 = load1(a0, j + 0); a02 = load1(a0, j + 1);
        a10 = load1(a1, j - 1);  a11 = load1(a1, j + 0); a12 = load1(a1, j + 1);
        a20 = load1(a2, j - 1);  a21 = load1(a2, j + 0); a22 = load1(a2, j + 1);

        l0 = LEFT_64(a00, a01); b0 = a01;r0 = RIGHT_64(a01, a02);

        l1 = LEFT_64(a10, a11);b1 = a11;r1 = RIGHT_64(a11, a12);
        
        l2 = LEFT_64(a20, a21);b2 = a21;r2 = RIGHT_64(a21,a22);

        y =MIN9(
            l0, b0, r0,
            l1, b1, r1,
            l2, b2, r2
            );

        store1(y0, j, y);
    default:
        break;
    }
}
// ----------------------------------------------------------------------------------------
void line_erosion3_ui64matrix_swp64_ilu3_red(uint64 **X, int i, int j0, int j1, uint64 **Y)
// ----------------------------------------------------------------------------------------
{
    // debug
    setborder1(j0-1, j1+1);

    uint64 *a0, *a1, *a2;
    a0 = X[i-1];
    a1 = X[i  ];
    a2 = X[i+1];

    uint64 y;
    uint64 *y0 = Y[i];
    
    uint64 l0,r0; 
    
    uint64 a00, a01, a02;
    uint64 a10, a11, a12;
    uint64 a20, a21, a22;

    uint64 c0 , c1 , c2; 
    int j = j0 ;
    a00 = load1(a0, j - 1); a01 = load1(a0, j + 0);  
    a10 = load1(a1, j - 1); a11 = load1(a1, j + 0); 
    a20 = load1(a2, j - 1); a21 = load1(a2, j + 0);
    
    c0 = MIN3(a00, a10, a20); c1 = MIN3(a01, a11, a21);

    int r = ((j1 - 2) % 3);
    for (j = j0; j <= ((j1 - 2) - r); j += 3){
        // loop 1 -------------------------------------------------------------------------------
        // load
        a02 = load1(a0, j + 1); 
        a12 = load1(a1, j + 1); 
        a22 = load1(a2, j + 1);
        
        c2 = MIN3(a02, a12, a22);

        // proc
        l0 = LEFT_64(c0,c1);
        r0 = RIGHT_64(c1, c2);

        y = MIN3(l0,c1,r0);

        // store
        store1(y0, j, y);

        // loop 2 -------------------------------------------------------------------------------
        // load
        a02 = load1(a0, j + 2); 
        a12 = load1(a1, j + 2); 
        a22 = load1(a2, j + 2);
        
        c0 = MIN3(a02, a12, a22);

        // proc
        l0 = LEFT_64(c1,c2);
        r0 = RIGHT_64 (c2, c0);

        y = MIN3(l0,c2,r0);

        // store
        store1(y0, j+1, y);

        // loop 3 -------------------------------------------------------------------------------
        // load
        a02 = load1(a0, j + 3); 
        a12 = load1(a1, j + 3); 
        a22 = load1(a2, j + 3);
        
        c1 = MIN3(a02, a12, a22);

        // proc
        l0 = LEFT_64(c2,c0);
        r0 = RIGHT_64 (c0, c1);

        y = MIN3(l0,c0,r0);

        // store
        store1(y0, j+2, y);

    }
    // epilogue
    j = ((j1 - 2) - r);
    j++;

    a00 = load1(a0, j - 1); a01 = load1(a0, j + 0);  
    a10 = load1(a1, j - 1); a11 = load1(a1, j + 0); 
    a20 = load1(a2, j - 1); a21 = load1(a2, j + 0);
    c0 = MIN3(a00, a10, a20); c1 = MIN3(a01, a11, a21);
    switch (r + 2)
    {
    case 4:
        // load
        a02 = load1(a0, j + 1); 
        a12 = load1(a1, j + 1); 
        a22 = load1(a2, j + 1);
        
        c2 = MIN3(a02, a12, a22);

        // proc
        l0 = LEFT_64(c0,c1);
        r0 = RIGHT_64 (c1, c2);

        y = MIN3(l0,c1,r0);

        // store
        store1(y0, j, y);

        // rotation de colonne
        c0 = c1;
        c1 = c2;
        j++;

    case 3:
        // load
        a02 = load1(a0, j + 1); 
        a12 = load1(a1, j + 1); 
        a22 = load1(a2, j + 1);
        
        c2 = MIN3(a02, a12, a22);

        // proc
        l0 = LEFT_64(c0,c1);
        r0 = RIGHT_64 (c1, c2);

        y = MIN3(l0,c1,r0);

        // store
        store1(y0, j, y);

        // rotation de colonne
        c0 = c1;
        c1 = c2;
        j++;

    case 2:
        // load
        a02 = load1(a0, j + 1); 
        a12 = load1(a1, j + 1); 
        a22 = load1(a2, j + 1);
        
        c2 = MIN3(a02, a12, a22);

        // proc
        l0 = LEFT_64(c0,c1);
        r0 = RIGHT_64 (c1, c2);

        y = MIN3(l0,c1,r0);

        // store
        store1(y0, j, y);

        // rotation de colonne
        c0 = c1;
        c1 = c2;
        j++;

    case 1:
        // load
        a02 = load1(a0, j + 1); 
        a12 = load1(a1, j + 1); 
        a22 = load1(a2, j + 1);
        
        c2 = MIN3(a02, a12, a22);

        // proc
        l0 = LEFT_64(c0,c1);
        r0 = RIGHT_64 (c1, c2);

        y = MIN3(l0,c1,r0);

        // store
        store1(y0, j, y);

    default:
        break;
    }
}
// ----------------------------------------------------------------------------------------
void line_erosion3_ui64matrix_swp64_elu2_red(uint64 **X, int i, int j0, int j1, uint64 **Y)
// ----------------------------------------------------------------------------------------
{
    /*
        Version avec deroulage de boucle externe de degree 2 
        et reduction de colonnes
    */

    // debug
    setborder1(j0-1, j1+1);

    uint64 *a0, *a1, *a2, *a3;
    a0 = X[i-1];
    a1 = X[i  ];
    a2 = X[i+1];
    a3 = X[i+2];

    uint64 y00, y01;
    uint64 *y0 = Y[i], *y1 = Y[i+1];

    // var des 12 cases de X
    uint64 a00, a01, a02;
    uint64 a10, a11, a12;
    uint64 a20, a21, a22;
    uint64 a30, a31, a32;

    // var des colonnes 
    uint64 c00, c10, c20;
    uint64 c01, c11, c21;

    uint64 l0,r0,l1,r1;

    // var temporaire pour la factorisation
    uint64 f0;

    int j = j0;

    // prologue
    a00 = load1(a0, j - 1); a01 = load1(a0, j + 0);  
    a10 = load1(a1, j - 1); a11 = load1(a1, j + 0); 
    a20 = load1(a2, j - 1); a21 = load1(a2, j + 0);
    a30 = load1(a3, j - 1); a31 = load1(a3, j + 0);

    // premieres colonnes
    c00 = MIN3(a00, a10, a20); c01 = MIN3(a10, a20, a30);
    c10 = MIN3(a01, a11, a21); c11 = MIN3(a11, a21, a31);

    for (j = j0; j <= j1; j++){
        // load
        a02 = load1(a0, j + 1);
        a12 = load1(a1, j + 1);
        a22 = load1(a2, j + 1);

        a32 = load1(a3, j + 1); // elu2

        // nouvelles colonnes
        c20 = MIN3(a02, a12, a22);

        l0 = LEFT_64(c00,c10);
        r0 = RIGHT_64 (c10, c20);

        c21 = MIN3(a12, a22, a32); // elu2

        l1 = LEFT_64(c01,c11);
        r1 = RIGHT_64 (c11, c21);

         // calc
        y00 = MIN3(l0, c10, r0);

        y01 = MIN3(l1, c11, r1); // elu2

        // store
        store1(y0, j, y00);

        store1(y1, j, y01); // elu2

        // rotation des colonnes
        c00 = c10;
        c10 = c20;

        c01 = c11;
        c11 = c21;
    }
}
// -----------------------------------------------------------------------------------------------
void line_erosion3_ui64matrix_swp64_elu2_red_factor(uint64 **X, int i, int j0, int j1, uint64 **Y)
// -----------------------------------------------------------------------------------------------
{
    /*
        Version avec deroulage de boucle externe de degree 2 
        et reduction de colonnes
    */

       // debug
        setborder1(j0-1, j1+1);

    uint64 *a0, *a1, *a2, *a3;
    a0 = X[i-1];
    a1 = X[i  ];
    a2 = X[i+1];
    a3 = X[i+2];

    uint64 y00, y01;
    uint64 *y0 = Y[i], *y1 = Y[i+1];

    // var des 12 cases de X
    uint64 a00, a01, a02;
    uint64 a10, a11, a12;
    uint64 a20, a21, a22;
    uint64 a30, a31, a32;

    // var des colonnes 
    uint64 c00, c10, c20;
    uint64 c01, c11, c21;

    uint64 l0,r0,l1,r1;

    // var temporaire pour la factorisation
    uint64 f0;

    int j = j0;

    // prologue
    a00 = load1(a0, j - 1); a01 = load1(a0, j + 0);  
    a10 = load1(a1, j - 1); a11 = load1(a1, j + 0); 
    a20 = load1(a2, j - 1); a21 = load1(a2, j + 0);
    a30 = load1(a3, j - 1); a31 = load1(a3, j + 0);

    // premieres colonnes
    f0  = MIN(a10, a20);
    c00 = MIN(a00, f0); c01 = MIN(f0, a30);
    f0  = MIN(a11, a21);
    c10 = MIN(a01, f0); c11 = MIN(f0, a31);

    for (j = j0; j <= j1; j++){
        // load
        a02 = load1(a0, j + 1);
        a12 = load1(a1, j + 1);
        a22 = load1(a2, j + 1);

        a32 = load1(a3, j + 1); // elu2

        // nouvelles colonnes
        f0 = MIN(a12, a22);
        c20 = MIN(a02, f0);

        l0 = LEFT_64(c00,c10);
        r0 = RIGHT_64 (c10, c20);

        c21 = MIN(f0, a32); // elu2

        l1 = LEFT_64(c01,c11);
        r1 = RIGHT_64 (c11, c21);

         // calc
        y00 = MIN3(l0, c10, r0);

        y01 = MIN3(l1, c11, r1); // elu2

        // store
        store1(y0, j, y00);

        store1(y1, j, y01); // elu2

        // rotation des colonnes
        c00 = c10;
        c10 = c20;

        c01 = c11;
        c11 = c21;
    }
}
// ---------------------------------------------------------------------------------------------
void line_erosion3_ui64matrix_swp64_ilu3_elu2_red(uint64 **X, int i, int j0, int j1, uint64 **Y)
// ---------------------------------------------------------------------------------------------
{
    /*
        Version avec deroulage de boucle externe de degre 3, 
        deroulage de boucle interne de degre 2 et reduction de colonnes 
    */
    // debug
    setborder1(j0-1, j1+1);

    uint64 *a0, *a1, *a2, *a3;
    a0 = X[i-1];
    a1 = X[i  ];
    a2 = X[i+1];
    a3 = X[i+2];

    uint64 y;
    uint64 *y0 = Y[i];
    uint64 *y1 = Y[i+1]; 
    
    uint64 l0,r0;
    uint64 l1, r1;
    
    uint64 a00, a01, a02;
    uint64 a10, a11, a12;
    uint64 a20, a21, a22;
    uint64 a30, a31, a32;

    uint64 c00, c10, c20;
    uint64 c01, c11, c21;

    int j = j0 ;
    a00 = load1(a0, j - 1); a01 = load1(a0, j + 0);  
    a10 = load1(a1, j - 1); a11 = load1(a1, j + 0); 
    a20 = load1(a2, j - 1); a21 = load1(a2, j + 0);
    a30 = load1(a3, j - 1); a31 = load1(a3, j + 0);
    
    c00 = MIN3(a00, a10, a20); c10 = MIN3(a01, a11, a21);
    c01 = MIN3(a10, a20, a30); c11 = MIN3(a11, a21, a31);

    int r = (((j1-j0)+1) % 3);
    for (j = j0; j <= (j1 - r); j += 3){
        // loop 1 -------------------------------------------------------------------------------
        // load
        a02 = load1(a0, j + 1); 
        a12 = load1(a1, j + 1); 
        a22 = load1(a2, j + 1);

        a32 = load1(a3, j + 1);
        
        c20 = MIN3(a02, a12, a22);
        c21 = MIN3(a12, a22, a32);

        // proc
        l0 = LEFT_64(c00,c10);
        r0 = RIGHT_64(c10, c20);

        l1 = LEFT_64(c01, c11);
        r1 = RIGHT_64(c11, c21);


        y = MIN3(l0,c10,r0);
        store1(y0, j, y);

        y = MIN3(l1, c11, r1);
        store1(y1, j, y);

        // loop 2 -------------------------------------------------------------------------------
        // load
        a02 = load1(a0, j + 2); 
        a12 = load1(a1, j + 2); 
        a22 = load1(a2, j + 2);

        a32 = load1(a3, j + 2);

        c00 = MIN3(a02, a12, a22);
        c01 = MIN3(a12, a22, a32);
        
        // proc
        l0 = LEFT_64(c10, c20);
        r0 = RIGHT_64 (c20, c00);
        l1 = LEFT_64(c11, c21);
        r1 = RIGHT_64(c21, c01);


        y = MIN3(l0,c20,r0);
        store1(y0, j+1, y);

        y = MIN3(l1, c21, r1);
        store1(y1, j+1, y);
        

        // loop 3 -------------------------------------------------------------------------------
        // load
        a02 = load1(a0, j + 3); 
        a12 = load1(a1, j + 3); 
        a22 = load1(a2, j + 3);

        a32 = load1(a3, j + 3);
        
        c10 = MIN3(a02, a12, a22);
        c11 = MIN3(a12, a22, a32);

        l0 = LEFT_64(c20,c00);
        r0 = RIGHT_64(c00, c10);

        l1 = LEFT_64(c21, c01);
        r1 = RIGHT_64(c01, c11);

        y = MIN3(l0,c00,r0);
        store1(y0, j+2, y);

        y = MIN3(l1, c01, r1);
        store1(y1, j+2, y);

    }
    // epilogue
    j = (j1- r);
    j++ ;
    a00 = load1(a0, j - 1); a01 = load1(a0, j + 0);  
    a10 = load1(a1, j - 1); a11 = load1(a1, j + 0); 
    a20 = load1(a2, j - 1); a21 = load1(a2, j + 0);
    a30 = load1(a3, j - 1); a31 = load1(a3, j + 0);
    
    c00 = MIN3(a00, a10, a20); c10 = MIN3(a01, a11, a21);
    c01 = MIN3(a10, a20, a30); c11 = MIN3(a11, a21, a31);
    
    switch (r)
    {
    case 4:
        // load
        a02 = load1(a0, j + 1); 
        a12 = load1(a1, j + 1); 
        a22 = load1(a2, j + 1);

        a32 = load1(a3, j + 1);
        
        c20 = MIN3(a02, a12, a22);
        c21 = MIN3(a12, a22, a32);

        // proc
        l0 = LEFT_64(c00,c10);
        r0 = RIGHT_64(c10, c20);
        l1 = LEFT_64(c01, c11);
        r1 = RIGHT_64(c11, c21);


        y = MIN3(l0,c10,r0);
        store1(y0, j, y);

        y = MIN3(l1, c11, r1);
        store1(y1, j, y);

        // rotation de colonne
        c00 = c10;
        c10 = c20;

        c01 = c11;
        c11 = c21;
        j++;

    case 3:
        // load
        a02 = load1(a0, j + 1); 
        a12 = load1(a1, j + 1); 
        a22 = load1(a2, j + 1);

        a32 = load1(a3, j + 1);
        
        c20 = MIN3(a02, a12, a22);
        c21 = MIN3(a12, a22, a32);

        // proc
        l0 = LEFT_64(c00,c10);
        r0 = RIGHT_64(c10, c20);
        l1 = LEFT_64(c01, c11);
        r1 = RIGHT_64(c11, c21);


        y = MIN3(l0,c10,r0);
        store1(y0, j, y);

        y = MIN3(l1, c11, r1);
        store1(y1, j, y);

        // rotation de colonne
        c00 = c10;
        c10 = c20;
        c01 = c11;
        c11 = c21;
        j++;

    case 2:
        // load
        a02 = load1(a0, j + 1); 
        a12 = load1(a1, j + 1); 
        a22 = load1(a2, j + 1);

        a32 = load1(a3, j + 1);
        
        c20 = MIN3(a02, a12, a22);
        c21 = MIN3(a12, a22, a32);

        // proc
        l0 = LEFT_64(c00,c10);
        r0 = RIGHT_64(c10, c20);
        l1 = LEFT_64(c01, c11);
        r1 = RIGHT_64(c11, c21);


        y = MIN3(l0,c10,r0);
        store1(y0, j, y);

        y = MIN3(l1, c11, r1);
        store1(y1, j, y);

        // rotation de colonne
        c00 = c10;
        c10 = c20;
        c01 = c11;
        c11 = c21;
        j++;

    case 1:
        // load
        a02 = load1(a0, j + 1); 
        a12 = load1(a1, j + 1); 
        a22 = load1(a2, j + 1);

        a32 = load1(a3, j + 1);
        
        c20 = MIN3(a02, a12, a22);
        c21 = MIN3(a12, a22, a32);

        // proc
        l0 = LEFT_64(c00,c10);
        r0 = RIGHT_64(c10, c20);
        l1 = LEFT_64(c01, c11);
        r1 = RIGHT_64(c11, c21);


        y = MIN3(l0,c10,r0);
        store1(y0, j, y);

        y = MIN3(l1, c11, r1);
        store1(y1, j, y);

        // rotation de colonne
        c00 = c10;
        c10 = c20;
        c01 = c11;
        c11 = c21;
        j++;
    default:
        break;
    }
}
// ----------------------------------------------------------------------------------------------------
void line_erosion3_ui64matrix_swp64_ilu3_elu2_red_factor(uint64 **X, int i, int j0, int j1, uint64 **Y)
// ----------------------------------------------------------------------------------------------------
{
    /*
        Version avec deroulage de boucle externe de degre 3, 
        deroulage de boucle interne de degre 2 et reduction de colonnes 
    */
    // debug
    setborder1(j0-1, j1+1);

    uint64 *a0, *a1, *a2, *a3;
    a0 = X[i-1];
    a1 = X[i  ];
    a2 = X[i+1];
    a3 = X[i+2];

    uint64 y;
    uint64 *y0 = Y[i];
    uint64 *y1 = Y[i+1]; 
    
    uint64 l0,r0;
    uint64 l1, r1;
    
    uint64 f; // factorisation

    uint64 a00, a01, a02;
    uint64 a10, a11, a12;
    uint64 a20, a21, a22;
    uint64 a30, a31, a32;

    uint64 c00, c10, c20;
    uint64 c01, c11, c21;

    int j = j0 ;
    a00 = load1(a0, j - 1); a01 = load1(a0, j + 0);  
    a10 = load1(a1, j - 1); a11 = load1(a1, j + 0); 
    a20 = load1(a2, j - 1); a21 = load1(a2, j + 0);
    a30 = load1(a3, j - 1); a31 = load1(a3, j + 0);
    
    f = MIN(a10,a20);
    c00 = MIN(a00, f); 
    c01 = MIN(f, a30); 
    
    f = MIN (a11, a21);
    c10 = MIN(a01, f);
    c11 = MIN(f, a31);

    int r = (((j1-j0)+1) % 3);
    for (j = j0; j <= (j1 - r); j += 3){
        // loop 1 -------------------------------------------------------------------------------
        // load
        a02 = load1(a0, j + 1); 
        a12 = load1(a1, j + 1); 
        a22 = load1(a2, j + 1);

        a32 = load1(a3, j + 1);
        
        f = MIN(a12, a22);
        c20 = MIN(a02, f);
        c21 = MIN(f, a32);

        // proc
        l0 = LEFT_64(c00,c10);
        r0 = RIGHT_64(c10, c20);

        l1 = LEFT_64(c01, c11);
        r1 = RIGHT_64(c11, c21);


        y = MIN3(l0,c10,r0);
        store1(y0, j, y);

        y = MIN3(l1, c11, r1);
        store1(y1, j, y);

        // loop 2 -------------------------------------------------------------------------------
        // load
        a02 = load1(a0, j + 2); 
        a12 = load1(a1, j + 2); 
        a22 = load1(a2, j + 2);

        a32 = load1(a3, j + 2);
        
        f= MIN(a12,a22);
        c00 = MIN(a02, f);
        c01 = MIN(f, a32);
        
        // proc
        l0 = LEFT_64(c10, c20);
        r0 = RIGHT_64 (c20, c00);
        l1 = LEFT_64(c11, c21);
        r1 = RIGHT_64(c21, c01);


        y = MIN3(l0,c20,r0);
        store1(y0, j+1, y);

        y = MIN3(l1, c21, r1);
        store1(y1, j+1, y);
        

        // loop 3 -------------------------------------------------------------------------------
        // load
        a02 = load1(a0, j + 3); 
        a12 = load1(a1, j + 3); 
        a22 = load1(a2, j + 3);

        a32 = load1(a3, j + 3);
        
        f = MIN(a12,a22);
        c10 = MIN(a02, f);
        c11 = MIN(f, a32);

        l0 = LEFT_64(c20,c00);
        r0 = RIGHT_64(c00, c10);

        l1 = LEFT_64(c21, c01);
        r1 = RIGHT_64(c01, c11);

        y = MIN3(l0,c00,r0);
        store1(y0, j+2, y);

        y = MIN3(l1, c01, r1);
        store1(y1, j+2, y);

    }
    // epilogue
    j = (j1- r);
    j++ ;
    a00 = load1(a0, j - 1); a01 = load1(a0, j + 0);  
    a10 = load1(a1, j - 1); a11 = load1(a1, j + 0); 
    a20 = load1(a2, j - 1); a21 = load1(a2, j + 0);
    a30 = load1(a3, j - 1); a31 = load1(a3, j + 0);
    
    f = MIN(a10,a20);
    c00 = MIN(a00, f); 
    c01 = MIN(f, a30); 
    
    f = MIN (a11, a21);
    c10 = MIN(a01, f);
    c11 = MIN(f, a31);
    
    switch (r)
    {
    case 4:
        // load
        a02 = load1(a0, j + 1); 
        a12 = load1(a1, j + 1); 
        a22 = load1(a2, j + 1);

        a32 = load1(a3, j + 1);
        
        f = MIN(a12,a22);
        c20 = MIN(a02, f);
        c21 = MIN(f, a32);

        // proc
        l0 = LEFT_64(c00,c10);
        r0 = RIGHT_64(c10, c20);
        l1 = LEFT_64(c01, c11);
        r1 = RIGHT_64(c11, c21);


        y = MIN3(l0,c10,r0);
        store1(y0, j, y);

        y = MIN3(l1, c11, r1);
        store1(y1, j, y);

        // rotation de colonne
        c00 = c10;
        c10 = c20;

        c01 = c11;
        c11 = c21;
        j++;

    case 3:
        // load
        a02 = load1(a0, j + 1); 
        a12 = load1(a1, j + 1); 
        a22 = load1(a2, j + 1);

        a32 = load1(a3, j + 1);
        
        f = MIN(a12,a22);
        c20 = MIN(a02, f);
        c21 = MIN(f, a32);

        // proc
        l0 = LEFT_64(c00,c10);
        r0 = RIGHT_64(c10, c20);
        l1 = LEFT_64(c01, c11);
        r1 = RIGHT_64(c11, c21);


        y = MIN3(l0,c10,r0);
        store1(y0, j, y);

        y = MIN3(l1, c11, r1);
        store1(y1, j, y);

        // rotation de colonne
        c00 = c10;
        c10 = c20;
        c01 = c11;
        c11 = c21;
        j++;

    case 2:
        // load
        a02 = load1(a0, j + 1); 
        a12 = load1(a1, j + 1); 
        a22 = load1(a2, j + 1);

        a32 = load1(a3, j + 1);
        
        f = MIN(a12,a22);
        c20 = MIN(a02, f);
        c21 = MIN(f, a32);

        // proc
        l0 = LEFT_64(c00,c10);
        r0 = RIGHT_64(c10, c20);
        l1 = LEFT_64(c01, c11);
        r1 = RIGHT_64(c11, c21);


        y = MIN3(l0,c10,r0);
        store1(y0, j, y);

        y = MIN3(l1, c11, r1);
        store1(y1, j, y);

        // rotation de colonne
        c00 = c10;
        c10 = c20;
        c01 = c11;
        c11 = c21;
        j++;

    case 1:
        // load
        a02 = load1(a0, j + 1); 
        a12 = load1(a1, j + 1); 
        a22 = load1(a2, j + 1);

        a32 = load1(a3, j + 1);
        
        f = MIN(a12,a22);
        c20 = MIN(a02, f);
        c21 = MIN(f, a32);

        // proc
        l0 = LEFT_64(c00,c10);
        r0 = RIGHT_64(c10, c20);
        l1 = LEFT_64(c01, c11);
        r1 = RIGHT_64(c11, c21);


        y = MIN3(l0,c10,r0);
        store1(y0, j, y);

        y = MIN3(l1, c11, r1);
        store1(y1, j, y);
    default:
        break;
    }
}
// -----------------------------------------------------------------------------------------
void erosion3_ui64matrix_swp64_basic(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// -----------------------------------------------------------------------------------------
{
    int i;
    for (i = i0; i <= i1; i++)
    {
        line_erosion3_ui64matrix_swp64_basic(X, i, j0, j1, Y);
    }
}
// ---------------------------------------------------------------------------------------
void erosion3_ui64matrix_swp64_reg(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// ---------------------------------------------------------------------------------------
{
    int i;
    for (i = i0; i <= i1; i++)
    {
        line_erosion3_ui64matrix_swp64_reg(X, i, j0, j1, Y);
    }
}
// ---------------------------------------------------------------------------------------
void erosion3_ui64matrix_swp64_rot(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// ---------------------------------------------------------------------------------------
{
    int i;
    for (i = i0; i <= i1; i++)
    {
        line_erosion3_ui64matrix_swp64_rot(X, i, j0, j1, Y);
    }
}
// ---------------------------------------------------------------------------------------
void erosion3_ui64matrix_swp64_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// ---------------------------------------------------------------------------------------
{
    int i;
    for (i = i0; i <= i1; i++)
    {
        line_erosion3_ui64matrix_swp64_red(X, i, j0, j1, Y);
    }
}
// ----------------------------------------------------------------------------------------
void erosion3_ui64matrix_swp64_ilu3(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// ----------------------------------------------------------------------------------------
{
    int i;
    for (i = i0; i <= i1; i++)
    {
        line_erosion3_ui64matrix_swp64_ilu3(X, i, j0, j1, Y);
    }
}
// --------------------------------------------------------------------------------------------
void erosion3_ui64matrix_swp64_ilu3_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// --------------------------------------------------------------------------------------------
{
    int i;
    for (i = i0; i <= i1; i++)
    {
        line_erosion3_ui64matrix_swp64_ilu3_red(X, i, j0, j1, Y);
    }
}
// --------------------------------------------------------------------------------------------
void erosion3_ui64matrix_swp64_elu2_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// --------------------------------------------------------------------------------------------
{
    int i;
    int r = ((i1-1)%2);
    for (i = i0; i <= ((i1-1)-r); i+=2)
    {
        line_erosion3_ui64matrix_swp64_elu2_red(X, i, j0, j1, Y);
    }
    i = ((i1-1)-r);
    switch(r+1)
    {
    case 2:
        i++;
        line_erosion3_ui64matrix_swp64_red(X, i, j0, j1, Y);
    case 1:
        i++;
        line_erosion3_ui64matrix_swp64_red(X, i, j0, j1, Y);    
    default:
        break;
    }
}
// -------------------------------------------------------------------------------------------------
void erosion3_ui64matrix_swp64_ilu3_elu2_red(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// -------------------------------------------------------------------------------------------------
{
    int i;
    int r = ((i1-1)%2);
    for (i = i0; i <= ((i1-1)-r); i+=2)
    {
        line_erosion3_ui64matrix_swp64_ilu3_elu2_red(X, i, j0, j1, Y);
    }

    i = ((i1-1)-r);
    switch (r+1)
    {
    case 2:
        i++;
        line_erosion3_ui64matrix_swp64_ilu3_red(X, i, j0, j1, Y);
    case 1:
        i++;
        line_erosion3_ui64matrix_swp64_ilu3_red(X, i, j0, j1, Y);
    
    default:
        break;
    }
}
// --------------------------------------------------------------------------------------------------------
void erosion3_ui64matrix_swp64_ilu3_elu2_red_factor(uint64 **X, int i0, int i1, int j0, int j1, uint64 **Y)
// --------------------------------------------------------------------------------------------------------
{
    int i;
    int r = ((i1-1)%2);
    for (i = i0; i <= ((i1-1)-r); i+=2)
    {
        line_erosion3_ui64matrix_swp64_ilu3_elu2_red_factor(X, i, j0, j1, Y);
    }

    i = ((i1-1)-r);
    switch (r+1)
    {
    case 2:
        i++;
        line_erosion3_ui64matrix_swp64_ilu3_red(X, i, j0, j1, Y);
    case 1:
        i++;
        line_erosion3_ui64matrix_swp64_ilu3_red(X, i, j0, j1, Y);

    default:
        break;
    }
}