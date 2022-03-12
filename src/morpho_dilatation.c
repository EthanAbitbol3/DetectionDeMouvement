/* --------------------------- */
/* --- morpho_dilatation.c --- */
/* --------------------------- */

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
#include "morpho_dilatation.h"

// -------------------------------------------------------------------------------
void line_dilatation3_ui8matrix_basic(uint8 **X, int i, int j0, int j1, uint8 **Y)
// -------------------------------------------------------------------------------
{
    /*
        Version basique
    */
    int j;
    for (j = j0; j <= j1; j++)
    {

        Y[i][j] = MAX9(
            X[i - 1][j - 1], X[i - 1][j + 0], X[i - 1][j + 1],
            X[i + 0][j - 1], X[i + 0][j + 0], X[i + 0][j + 1],
            X[i + 1][j - 1], X[i + 1][j + 0], X[i + 1][j + 1]);
    }
}
// -----------------------------------------------------------------------------
void line_dilatation3_ui8matrix_reg(uint8 **X, int i, int j0, int j1, uint8 **Y)
// -----------------------------------------------------------------------------
{
    /*
        Version avec optimitation en  registres 
    */

    // var des 3 lignes
    uint8 *a0 = X[i - 1];
    uint8 *a1 = X[i + 0];
    uint8 *a2 = X[i + 1];

    uint8 *Yi = Y[i];

    // var des 9 cases de X
    uint8 a00, a01, a02;
    uint8 a10, a11, a12;
    uint8 a20, a21, a22;

    uint8 y;

    // debug
    setborder1(j0, j1);

    int j;
    for (j = j0; j <= j1; j++)
    {
        // load
        a00 = load1(a0, j - 1);
        a10 = load1(a1, j - 1);
        a20 = load1(a2, j - 1);
        a01 = load1(a0, j + 0);
        a11 = load1(a1, j + 0);
        a21 = load1(a2, j + 0);
        a02 = load1(a0, j + 1);
        a12 = load1(a1, j + 1);
        a22 = load1(a2, j + 1);

        // calc
        y = MAX9(
            a00, a01, a02,
            a10, a11, a12,
            a20, a21, a22);

        // store
        store1(Yi, j, y);
    }
}
// -----------------------------------------------------------------------------
void line_dilatation3_ui8matrix_rot(uint8 **X, int i, int j0, int j1, uint8 **Y)
// -----------------------------------------------------------------------------
{
    /*
        Version avec rotation de registre
    */

    // var des 3 lignes
    uint8 *a0 = X[i - 1];
    uint8 *a1 = X[i + 0];
    uint8 *a2 = X[i + 1];

    uint8 *Yi = Y[i];

    // var des 9 cases de X
    uint8 a00, a01, a02;
    uint8 a10, a11, a12;
    uint8 a20, a21, a22;

    uint8 y;

    // debug
    setborder1(j0, j1);

    int j = j0;

    // load des 2 premieres cases de chaques lignes
    a00 = load1(a0, j - 1);
    a10 = load1(a1, j - 1);
    a20 = load1(a2, j - 1);
    a01 = load1(a0, j + 0);
    a11 = load1(a1, j + 0);
    a21 = load1(a2, j + 0);

    for (j = j0; j <= j1; j++)
    {
        // load de la nouvelle case de chaques lignes
        a02 = load1(a0, j + 1);
        a12 = load1(a1, j + 1);
        a22 = load1(a2, j + 1);

        // calc
        y = MAX9(
            a00, a01, a02,
            a10, a11, a12,
            a20, a21, a22);

        // store
        store1(Yi, j, y);

        // rotation de registre
        a00 = a01;
        a10 = a11;
        a20 = a21;

        a01 = a02;
        a11 = a12;
        a21 = a22;
    }
}
// -----------------------------------------------------------------------------
void line_dilatation3_ui8matrix_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// -----------------------------------------------------------------------------
{
    /*
        Version avec Reduction de colonnes 
    */

    // var des 3 lignes
    uint8 *a0 = X[i - 1];
    uint8 *a1 = X[i + 0];
    uint8 *a2 = X[i + 1];

    uint8 *Yi = Y[i];

    // var des 9 cases de X
    uint8 a00, a01, a02;
    uint8 a10, a11, a12;
    uint8 a20, a21, a22;

    uint8 c0, c1, c2;
    uint8 y;

    // debug
    setborder1(j0, j1);

    int j = j0;

    // load des 2 premieres cases de chaques lignes
    a00 = load1(a0, j - 1);
    a10 = load1(a1, j - 1);
    a20 = load1(a2, j - 1);
    a01 = load1(a0, j + 0);
    a11 = load1(a1, j + 0);
    a21 = load1(a2, j + 0);

    // 2 premieres colonnes
    c0 = MAX3(a00, a10, a20);
    c1 = MAX3(a01, a11, a21);

    for (j = j0; j <= j1; j++)
    {
        // load
        a02 = load1(a0, j + 1);
        a12 = load1(a1, j + 1);
        a22 = load1(a2, j + 1);

        // nouvelle colonne
        c2 = MAX3(a02, a12, a22);

        // calc
        y = MAX3(c0, c1, c2);

        // store
        store1(Yi, j, y);

        // rotation de colonne
        c0 = c1;
        c1 = c2;
    }
}
// ------------------------------------------------------------------------------
void line_dilatation3_ui8matrix_ilu3(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ------------------------------------------------------------------------------
{
    /*
        Version avec deroulage de boucles interne de degree 3
    */

    // var des 3 lignes
    uint8 *a0 = X[i - 1];
    uint8 *a1 = X[i + 0];
    uint8 *a2 = X[i + 1];

    uint8 *Yi = Y[i];

    // var des 9 cases de X
    uint8 a00, a01, a02;
    uint8 a10, a11, a12;
    uint8 a20, a21, a22;

    uint8 y;

    // debug
    setborder1(j0, j1);
    int j;

    // determination de l'intervalle de j
    int r = ((j1 - 2) % 3);

    for (j = j0; j <= ((j1 - 2) - r); j += 3)
    {
        // loop 1 -------------------------------------------------------------------------------
        // load
        a00 = load1(a0, j - 1);
        a10 = load1(a1, j - 1);
        a20 = load1(a2, j - 1);
        a01 = load1(a0, j + 0);
        a11 = load1(a1, j + 0);
        a21 = load1(a2, j + 0);
        a02 = load1(a0, j + 1);
        a12 = load1(a1, j + 1);
        a22 = load1(a2, j + 1);

        // calc
        y = MAX9(
            a00, a01, a02,
            a10, a11, a12,
            a20, a21, a22);

        // store
        store1(Yi, j + 0, y);

        // loop 2 -------------------------------------------------------------------------------
        // load

        a00 = load1(a0, j + 0);
        a10 = load1(a1, j + 0);
        a20 = load1(a2, j + 0);
        a01 = load1(a0, j + 1);
        a11 = load1(a1, j + 1);
        a21 = load1(a2, j + 1);
        a02 = load1(a0, j + 2);
        a12 = load1(a1, j + 2);
        a22 = load1(a2, j + 2);

        // calc
        y = MAX9(
            a00, a01, a02,
            a10, a11, a12,
            a20, a21, a22);

        // store
        store1(Yi, j + 1, y);

        // loop 3 -------------------------------------------------------------------------------
        // load
        a00 = load1(a0, j + 1);
        a10 = load1(a1, j + 1);
        a20 = load1(a2, j + 1);
        a01 = load1(a0, j + 2);
        a11 = load1(a1, j + 2);
        a21 = load1(a2, j + 2);
        a02 = load1(a0, j + 3);
        a12 = load1(a1, j + 3);
        a22 = load1(a2, j + 3);

        // calc
        y = MAX9(
            a00, a01, a02,
            a10, a11, a12,
            a20, a21, a22);

        // store
        store1(Yi, j + 2, y);
    }

    // epilogue
    j = ((j1 - 2) - r);
    switch (r + 2)
    {
    case 4:
        j++;
        // load
        a00 = load1(a0, j - 1);
        a10 = load1(a1, j - 1);
        a20 = load1(a2, j - 1);
        a01 = load1(a0, j + 0);
        a11 = load1(a1, j + 0);
        a21 = load1(a2, j + 0);
        a02 = load1(a0, j + 1);
        a12 = load1(a1, j + 1);
        a22 = load1(a2, j + 1);

        // calc
        y = MAX9(
            a00, a01, a02,
            a10, a11, a12,
            a20, a21, a22);

        // store
        store1(Yi, j + 0, y);
    case 3:
        j++;
        // load
        a00 = load1(a0, j - 1);
        a10 = load1(a1, j - 1);
        a20 = load1(a2, j - 1);
        a01 = load1(a0, j + 0);
        a11 = load1(a1, j + 0);
        a21 = load1(a2, j + 0);
        a02 = load1(a0, j + 1);
        a12 = load1(a1, j + 1);
        a22 = load1(a2, j + 1);

        // calc
        y = MAX9(
            a00, a01, a02,
            a10, a11, a12,
            a20, a21, a22);

        // store
        store1(Yi, j + 0, y);

    case 2:
        j++;
        // load
        a00 = load1(a0, j - 1);
        a10 = load1(a1, j - 1);
        a20 = load1(a2, j - 1);
        a01 = load1(a0, j + 0);
        a11 = load1(a1, j + 0);
        a21 = load1(a2, j + 0);
        a02 = load1(a0, j + 1);
        a12 = load1(a1, j + 1);
        a22 = load1(a2, j + 1);

        // calc
        y = MAX9(
            a00, a01, a02,
            a10, a11, a12,
            a20, a21, a22);

        // store
        store1(Yi, j + 0, y);

    case 1:
        j++;
        // load
        a00 = load1(a0, j - 1);
        a10 = load1(a1, j - 1);
        a20 = load1(a2, j - 1);
        a01 = load1(a0, j + 0);
        a11 = load1(a1, j + 0);
        a21 = load1(a2, j + 0);
        a02 = load1(a0, j + 1);
        a12 = load1(a1, j + 1);
        a22 = load1(a2, j + 1);

        // calc
        y = MAX9(
            a00, a01, a02,
            a10, a11, a12,
            a20, a21, a22);

        // store
        store1(Yi, j + 0, y);
    default:
        break;
    }
}
// ----------------------------------------------------------------------------------
void line_dilatation3_ui8matrix_ilu3_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------------------
{
    /*
        Version avec deroulage de boucle interne de degree 3 
        et reduction de colonnes 
    */

    // var des 3 lignes
    uint8 *a0 = X[i - 1];
    uint8 *a1 = X[i + 0];
    uint8 *a2 = X[i + 1];

    uint8 *Yi = Y[i];

    // var des 9 cases de X
    uint8 a00, a01, a02;
    uint8 a10, a11, a12;
    uint8 a20, a21, a22;

    // var des colonnes pour la reduction
    uint8 c1, c2, c3;

    uint8 y;

    // debug
    setborder1(j0, j1);

    // prologue rot
    int j = j0;

    // load
    a00 = load1(a0, j - 1); a10 = load1(a1, j - 1); a20 = load1(a2, j - 1);
    a01 = load1(a0, j + 0); a11 = load1(a1, j + 0); a21 = load1(a2, j + 0);

    // 2 premieres colonnes
    c1 = MAX3(a00, a10, a20);
    c2 = MAX3(a01, a11, a21);

    int r = ((j1 - 2) % 3);

    for (j = j0; j <= ((j1 - 2) - r); j += 3)
    {
        // loop 1 -------------------------------------------------------------------------------
        // load
        a02 = load1(a0, j + 1);
        a12 = load1(a1, j + 1);
        a22 = load1(a2, j + 1);

        // nouvelle colonne
        c3 = MAX3(a02, a12, a22);

        // calc
        y = MAX3(c1, c2, c3);

        // store
        store1(Yi, j, y);

        // loop 2 -------------------------------------------------------------------------------
        // load
        a02 = load1(a0, j + 2);
        a12 = load1(a1, j + 2);
        a22 = load1(a2, j + 2);

        // nouvelle colonne
        c1 = MAX3(a02, a12, a22);
        // calc
        y = MAX3(c2, c3, c1);
        // store
        store1(Yi, j + 1, y);

        // loop 3 -------------------------------------------------------------------------------
        // load
        a02 = load1(a0, j + 3);
        a12 = load1(a1, j + 3);
        a22 = load1(a2, j + 3);

        // nouvelle colonne
        c2 = MAX3(a02, a12, a22);
        // calc
        y = MAX3(c3, c1, c2);
        // store
        store1(Yi, j + 2, y);


    }

    // epilogue
    
    j = ((j1 - 2) - r);
    j++;

    //load
    a00 = load1(a0, j - 1);a01 = load1(a0, j + 0);
    a10 = load1(a1, j - 1);a11 = load1(a1, j + 0);
    a20 = load1(a2, j - 1);a21 = load1(a2, j + 0);

    // 2 premieres colonnes
    c1 = MAX3(a00, a10, a20);
    c2 = MAX3(a01, a11, a21);

    switch (r + 2)
    {
    case 4:
        // load
        a02 = load1(a0, j + 1);
        a12 = load1(a1, j + 1);
        a22 = load1(a2, j + 1);

        // nouvelle colonne
        c3 = MAX3(a02, a12, a22);

        // calc
        y = MAX3(c1, c2, c3);

        // store
        store1(Yi, j, y);

        // colonne rot
        c1 = c2;
        c2 = c3;
        j++;

    case 3:
        // load
        a02 = load1(a0, j + 1);
        a12 = load1(a1, j + 1);
        a22 = load1(a2, j + 1);

        // nouvelle colonne
        c3 = MAX3(a02, a12, a22);

        // calc
        y = MAX3(c1, c2, c3);

        // store
        store1(Yi, j, y);

        // colonne rot
        c1 = c2;
        c2 = c3;
        j++;

    case 2:
        // load
        a02 = load1(a0, j + 1);
        a12 = load1(a1, j + 1);
        a22 = load1(a2, j + 1);

        // nouvelle colonne
        c3 = MAX3(a02, a12, a22);

        // calc
        y = MAX3(c1, c2, c3);

        // store
        store1(Yi, j, y);

        // colonne rot
        c1 = c2;
        c2 = c3;
        j++;

    case 1:
        // load
        a02 = load1(a0, j + 1);
        a12 = load1(a1, j + 1);
        a22 = load1(a2, j + 1);

        // nouvelle colonne
        c3 = MAX3(a02, a12, a22);

        // calc
        y = MAX3(c1, c2, c3);

        // store
        store1(Yi, j, y);
    default:
        break;
    }
}
// ----------------------------------------------------------------------------------
void line_dilatation3_ui8matrix_elu2_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------------------
{
    /*
        Version avec deroulage de boucle externe de degree 2 
        et reduction de colonnes
    */

    // var des 4 lignes
    uint8 *a0 = X[i - 1];
    uint8 *a1 = X[i + 0];
    uint8 *a2 = X[i + 1];
    uint8 *a3 = X[i + 2];

    uint8 *Yi0 = Y[i + 0];
    uint8 *Yi1 = Y[i + 1];

    // var des 12 cases de X
    uint8 a00, a01, a02;
    uint8 a10, a11, a12;
    uint8 a20, a21, a22;
    uint8 a30, a31, a32;

    // var des colonnes
    uint8 c00, c10, c20;
    uint8 c01, c11, c21;

    uint8 y0, y1;

    // debug
    setborder1(j0, j1);

    int j = j0;

    // prologue
    // Load
    a00 = load1(a0, j - 1);
    a10 = load1(a1, j - 1);
    a20 = load1(a2, j - 1);
    a30 = load1(a3, j - 1);
    a01 = load1(a0, j + 0);
    a11 = load1(a1, j + 0);
    a21 = load1(a2, j + 0);
    a31 = load1(a3, j + 0);

    // colonnes
    c00 = MAX3(a00, a10, a20);
    c01 = MAX3(a10, a20, a30);
    c10 = MAX3(a01, a11, a21);
    c11 = MAX3(a11, a21, a31);

    for (j = j0; j <= j1; j++)
    {
        // load
        a02 = load1(a0, j + 1);
        a12 = load1(a1, j + 1);
        a22 = load1(a2, j + 1);

        a32 = load1(a3, j + 1); // elu2

        // nouvelles colonnes
        c20 = MAX3(a02, a12, a22);

        c21 = MAX3(a12, a22, a32); // elu2

        // calc
        y0 = MAX3(c00, c10, c20);

        y1 = MAX3(c01, c11, c21); // elu2

        // store
        store1(Yi0, j, y0);

        store1(Yi1, j, y1); // elu2

        // rotation des colonnes
        c00 = c10;
        c10 = c20;

        c01 = c11;
        c11 = c21;
    }
}
// -----------------------------------------------------------------------------------------
void line_dilatation3_ui8matrix_elu2_red_factor(uint8 **X, int i, int j0, int j1, uint8 **Y)
// -----------------------------------------------------------------------------------------
{
    /*
        Version avec deroulage de boucle externe de degree 2 
        et reduction de colonnes en version factorise
    */

    // Var des 4 lignes
    uint8 *a0 = X[i - 1];
    uint8 *a1 = X[i + 0];
    uint8 *a2 = X[i + 1];
    uint8 *a3 = X[i + 2];

    uint8 *Yi0 = Y[i + 0];
    uint8 *Yi1 = Y[i + 1];

    // Var des 12 cases de X
    uint8 a00, a01, a02;
    uint8 a10, a11, a12;
    uint8 a20, a21, a22;
    uint8 a30, a31, a32;

    // var temporaire pour la factorisation
    uint8 f0;

    uint8 c00, c10, c20;
    uint8 c01, c11, c21;

    uint8 y0, y1;

    // debug
    setborder1(j0, j1);

    int j = j0;

    // load
    a00 = load1(a0, j - 1);
    a10 = load1(a1, j - 1);
    a20 = load1(a2, j - 1);
    a30 = load1(a3, j - 1);
    a01 = load1(a0, j + 0);
    a11 = load1(a1, j + 0);
    a21 = load1(a2, j + 0);
    a31 = load1(a3, j + 0);

    // premieres colonnes
    f0 = MAX(a10, a20);
    c00 = MAX(a00, f0);
    c01 = MAX(f0, a30);
    f0 = MAX(a11, a21);
    c10 = MAX(a01, f0);
    c11 = MAX(f0, a31);

    for (j = j0; j <= j1; j++)
    {
        // load
        a02 = load1(a0, j + 1);
        a12 = load1(a1, j + 1);
        a22 = load1(a2, j + 1);

        a32 = load1(a3, j + 1);
        // nouvelle colonne

        f0 = MAX(a12, a22);

        c20 = MAX(a02, f0);
        c21 = MAX(f0, a32);

        // calc
        y0 = MAX3(c00, c10, c20);

        y1 = MAX3(c01, c11, c21);

        // store
        store1(Yi0, j, y0);

        store1(Yi1, j, y1);

        // rotation de colonne
        c00 = c10;
        c10 = c20;

        c01 = c11;
        c11 = c21;
    }
}
// ---------------------------------------------------------------------------------------
void line_dilatation3_ui8matrix_ilu3_elu2_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------------------
{
    /*
        Version avec deroulage de boucle externe de degre 3, 
        deroulage de boucle interne de degre 2 et reduction de colonnes 
    */

    // Var des 4 lignes 
    uint8 *a0 = X[i - 1];
    uint8 *a1 = X[i + 0];
    uint8 *a2 = X[i + 1];
    uint8 *a3 = X[i + 2];

    uint8 *Yi0 = Y[i + 0];
    uint8 *Yi1 = Y[i + 1];

    // Var des 12 cases de X
    uint8 a00, a01, a02;
    uint8 a10, a11, a12;
    uint8 a20, a21, a22;
    uint8 a30, a31, a32;

    // colonnes 
    uint8 c00, c10, c20;
    uint8 c01, c11, c21;

    uint8 y0, y1;

    // debug
    setborder1(j0, j1);

    // prologue
    int j = j0;

    a00 = load1(a0, j - 1); a10 = load1(a1, j - 1); a20 = load1(a2, j - 1); a30 = load1(a3, j - 1);
    a01 = load1(a0, j + 0); a11 = load1(a1, j + 0); a21 = load1(a2, j + 0); a31 = load1(a3, j + 0);

    // premieres colonnes
    c00 = MAX3(a00, a10, a20);c01 = MAX3(a10, a20, a30);
    c10 = MAX3(a01, a11, a21); c11 = MAX3(a11, a21, a31);

    int r = ((j1 - 2) % 3);

    for (j = j0; j <= ((j1 - 2) - r); j += 3)
    {
        // loop 1 -------------------------------------------------------------------------------
        // load
        // load
        a02 = load1(a0, j + 1);
        a12 = load1(a1, j + 1);
        a22 = load1(a2, j + 1);

        a32 = load1(a3, j + 1);
        // nouvelle colonne
        c20 = MAX3(a02, a12, a22);

        c21 = MAX3(a12, a22, a32);

        // calc
        y0 = MAX3(c00, c10, c20);

        y1 = MAX3(c01, c11, c21);

        // store
        store1(Yi0, j, y0);

        store1(Yi1, j, y1);
        
        // loop 2 -------------------------------------------------------------------------------
        // load
        a02 = load1(a0, j + 2);
        a12 = load1(a1, j + 2);
        a22 = load1(a2, j + 2);

        a32 = load1(a3, j + 2);
        // nouvelle colonne
        c00 = MAX3(a02, a12, a22);

        c01 = MAX3(a12, a22, a32);

        // calc
        // y0 = MAX3(c00, c10, c20);
        y0 = MAX3(c10, c20, c00);

        // y1 = MAX3(c01, c11, c21);
        y1 = MAX3(c11, c21, c01);


        // store
        store1(Yi0, j + 1, y0);

        store1(Yi1, j + 1, y1);

        // loop 3 -------------------------------------------------------------------------------
        // load
        // load
        a02 = load1(a0, j + 3);
        a12 = load1(a1, j + 3);
        a22 = load1(a2, j + 3);

        a32 = load1(a3, j + 3);
        // nouvelle colonne
        c10 = MAX3(a02, a12, a22);

        c11 = MAX3(a12, a22, a32);

        // calc
        y0 = MAX3(c20, c00, c10);
        y1 = MAX3(c21, c01, c11);
        // store
        store1(Yi0, j + 2, y0);

        store1(Yi1, j + 2, y1);
    }

    // epilogue
    j = ((j1 - 2) - r);
    j++;

    a00 = load1(a0, j - 1); a10 = load1(a1, j - 1); a20 = load1(a2, j - 1); a30 = load1(a3, j - 1);
    a01 = load1(a0, j + 0); a11 = load1(a1, j + 0); a21 = load1(a2, j + 0); a31 = load1(a3, j + 0);

    // 2 premieres colonnes
    c00 = MAX3(a00, a10, a20); c01 = MAX3(a10, a20, a30);
    c10 = MAX3(a01, a11, a21); c11 = MAX3(a11, a21, a31);

    switch (r + 2)
    {
    case 4:

        // load
        a02 = load1(a0, j + 1); a12 = load1(a1, j + 1); a22 = load1(a2, j + 1);

        a32 = load1(a3, j + 1);
        // nouvelle colonne
        c20 = MAX3(a02, a12, a22);

        c21 = MAX3(a12, a22, a32);

        // calc
        y0 = MAX3(c00, c10, c20);

        y1 = MAX3(c01, c11, c21);

        // store
        store1(Yi0, j, y0);

        store1(Yi1, j, y1);

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
        // nouvelle colonne
        c20 = MAX3(a02, a12, a22);

        c21 = MAX3(a12, a22, a32);

        // calc
        y0 = MAX3(c00, c10, c20);

        y1 = MAX3(c01, c11, c21);

        // store
        store1(Yi0, j, y0);

        store1(Yi1, j, y1);

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
        // nouvelle colonne
        c20 = MAX3(a02, a12, a22);

        c21 = MAX3(a12, a22, a32);

        // calc
        y0 = MAX3(c00, c10, c20);

        y1 = MAX3(c01, c11, c21);

        // store
        store1(Yi0, j, y0);

        store1(Yi1, j, y1);

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
        // nouvelle colonne
        c20 = MAX3(a02, a12, a22);

        c21 = MAX3(a12, a22, a32);

        // calc
        y0 = MAX3(c00, c10, c20);

        y1 = MAX3(c01, c11, c21);

        // store
        store1(Yi0, j, y0);

        store1(Yi1, j, y1);

    default:
        break;
    }
}
// ----------------------------------------------------------------------------------------------
void line_dilatation3_ui8matrix_ilu3_elu2_red_factor(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------------------------------
{
    /*
        Version avec deroulage de boucle interne de degre 3,
        deroulage de boucle externe de degre 2 et reduction de colonnes en version factorise.
    */

    // Var des 4 lignes 
    uint8 *a0 = X[i - 1];
    uint8 *a1 = X[i + 0];
    uint8 *a2 = X[i + 1];
    uint8 *a3 = X[i + 2];

    uint8 *Yi0 = Y[i + 0];
    uint8 *Yi1 = Y[i + 1];

    // Var des 12 cases de X 
    uint8 a00, a01, a02;
    uint8 a10, a11, a12;
    uint8 a20, a21, a22;
    uint8 a30, a31, a32;

    
    uint8 f0; // var temporaire pour la factorisation

    //  colonnes 
    uint8 c00, c10, c20;
    uint8 c01, c11, c21;

    uint8 y0, y1;

    // debug
    setborder1(j0, j1);

    // prologue
    int j = j0;

    a00 = load1(a0, j - 1); a10 = load1(a1, j - 1); a20 = load1(a2, j - 1); a30 = load1(a3, j - 1);
    a01 = load1(a0, j + 0); a11 = load1(a1, j + 0); a21 = load1(a2, j + 0); a31 = load1(a3, j + 0);

    // 2 premieres colonnes
    f0  = MAX(a10, a20);
    c00 = MAX(a00, f0); c01 = MAX(f0, a30);
    f0  = MAX(a11, a21);
    c10 = MAX(a01, f0); c11 = MAX(f0, a31);

    int r = ((j1 - 2) % 3);

    for (j = j0; j <= ((j1 - 2) - r); j += 3)
    {
        // loop 1 -------------------------------------------------------------------------------
        // load
        // load
        a02 = load1(a0, j + 1);
        a12 = load1(a1, j + 1);
        a22 = load1(a2, j + 1);

        a32 = load1(a3, j + 1);

        // nouvelle colonne
        f0  = MAX(a12, a22);
        c20 = MAX(a02, f0);
        c21 = MAX(f0, a32);

        // calc
        y0 = MAX3(c00, c10, c20);

        y1 = MAX3(c01, c11, c21);

        // store
        store1(Yi0, j, y0);

        store1(Yi1, j, y1);

        // loop 2 -------------------------------------------------------------------------------
        // load
        a02 = load1(a0, j + 2);
        a12 = load1(a1, j + 2);
        a22 = load1(a2, j + 2);

        a32 = load1(a3, j + 2);
        // nouvelle colonne
        f0  = MAX(a12, a22);
        c00 = MAX(a02, f0);
        c01 = MAX(f0, a32);

        // calc
        y0 = MAX3(c10, c20, c00);
        y1 = MAX3(c11, c21, c01);

        // store
        store1(Yi0, j + 1, y0);
        store1(Yi1, j + 1, y1);

        // loop 3 -------------------------------------------------------------------------------
        // load
        // load
        a02 = load1(a0, j + 3);
        a12 = load1(a1, j + 3);
        a22 = load1(a2, j + 3);

        a32 = load1(a3, j + 3);
        // nouvelle colonne
        f0  = MAX(a12, a22);
        c10 = MAX(a02, f0);
        c11 = MAX(f0, a32);

        // calc
        y0 = MAX3(c20, c00, c10);
        y1 = MAX3(c21, c01, c11);

        // store
        store1(Yi0, j + 2, y0);
        store1(Yi1, j + 2, y1);
    }

    // epilogue
    j = ((j1 - 2) - r);
    j++;
    a00 = load1(a0, j - 1);
    a01 = load1(a0, j + 0);
    a10 = load1(a1, j - 1);
    a11 = load1(a1, j + 0);
    a20 = load1(a2, j - 1);
    a21 = load1(a2, j + 0);

    a30 = load1(a3, j - 1);
    a31 = load1(a3, j + 0);

    // 2 premieres colonnes
    f0  = MAX(a10, a20);
    c00 = MAX(a00, f0); c01 = MAX(f0, a30);
    f0  = MAX(a11, a21);
    c10 = MAX(a01, f0); c11 = MAX(f0, a31);

    switch (r + 2)
    {
    case 4:

        // load
        a02 = load1(a0, j + 1);
        a12 = load1(a1, j + 1);
        a22 = load1(a2, j + 1);

        a32 = load1(a3, j + 1);
        // nouvelle colonne
        f0  = MAX(a12, a22);
        c20 = MAX(a02, f0);
        c21 = MAX(f0, a32);

        // calc
        y0 = MAX3(c00, c10, c20);

        y1 = MAX3(c01, c11, c21);

        // store
        store1(Yi0, j, y0);

        store1(Yi1, j, y1);

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
        // nouvelle colonne
        f0  = MAX(a12, a22);
        c20 = MAX(a02, f0);
        c21 = MAX(f0, a32);

        // calc
        y0 = MAX3(c00, c10, c20);

        y1 = MAX3(c01, c11, c21);

        // store
        store1(Yi0, j, y0);

        store1(Yi1, j, y1);

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
        // nouvelle colonne
        f0  = MAX(a12, a22);
        c20 = MAX(a02, f0);
        c21 = MAX(f0, a32);

        // calc
        y0 = MAX3(c00, c10, c20);

        y1 = MAX3(c01, c11, c21);

        // store
        store1(Yi0, j, y0);

        store1(Yi1, j, y1);

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
        // nouvelle colonne
        f0  = MAX(a12, a22);
        c20 = MAX(a02, f0);
        c21 = MAX(f0, a32);

        // calc
        y0 = MAX3(c00, c10, c20);

        y1 = MAX3(c01, c11, c21);

        // store
        store1(Yi0, j, y0);

        store1(Yi1, j, y1);

    default:
        break;
    }
}
// -----------------------------------------------------------------------------------
void dilatation3_ui8matrix_basic(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// -----------------------------------------------------------------------------------
{
    int i;
    for (i = i0; i <= i1; i++)
    {
        line_dilatation3_ui8matrix_basic(X, i, j0, j1, Y);
    }
}
// ---------------------------------------------------------------------------------
void dilatation3_ui8matrix_reg(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------------
{
    int i;
    for (i = i0; i <= i1; i++)
    {
        line_dilatation3_ui8matrix_reg(X, i, j0, j1, Y);
    }
}
// ---------------------------------------------------------------------------------
void dilatation3_ui8matrix_rot(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------------
{
    int i;
    for (i = i0; i <= i1; i++)
    {
        line_dilatation3_ui8matrix_rot(X, i, j0, j1, Y);
    }
}
// ---------------------------------------------------------------------------------
void dilatation3_ui8matrix_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------------
{
    int i;
    for (i = i0; i <= i1; i++)
    {
        line_dilatation3_ui8matrix_red(X, i, j0, j1, Y);
    }
}
// --------------------------------------------------------------------------------------
void dilatation3_ui8matrix_ilu3                (uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------------------
{
    int i;
    for (i = i0; i <= i1; i++)
    {
        line_dilatation3_ui8matrix_ilu3(X, i, j0, j1, Y);
    }
}

// --------------------------------------------------------------------------------------
void dilatation3_ui8matrix_ilu3_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------------------
{
    int i;
    for (i = i0; i <= i1; i++)
    {
        line_dilatation3_ui8matrix_ilu3_red(X, i, j0, j1, Y);
    }
}
// --------------------------------------------------------------------------------------
void dilatation3_ui8matrix_elu2_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------------------
{
    int i;
    int r = ((i1 - 1) % 2);
    for (i = i0; i <= ((i1 - 1) - r); i += 2)
    {
        line_dilatation3_ui8matrix_elu2_red(X, i, j0, j1, Y);
    }
    i = ((i1 - 1) - r);
    switch (r + 1)
    {
    case 2:
        i++;
        line_dilatation3_ui8matrix_red(X, i, j0, j1, Y);
    case 1:
        i++;
        line_dilatation3_ui8matrix_red(X, i, j0, j1, Y);
    default:
        break;
    }
}
// ---------------------------------------------------------------------------------------------
void dilatation3_ui8matrix_elu2_red_factor(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------------------------
{
    int i;
    int r = ((i1 - 1) % 2);
    for (i = i0; i <= ((i1 - 1) - r); i += 2)
    {
        line_dilatation3_ui8matrix_elu2_red_factor(X, i, j0, j1, Y);
    }

    i = ((i1 - 1) - r);
    switch (r + 1)
    {
    case 2:
        i++;
        line_dilatation3_ui8matrix_red(X, i, j0, j1, Y);
    case 1:
        i++;
        line_dilatation3_ui8matrix_red(X, i, j0, j1, Y);
    default:
        break;
    }
}
// -------------------------------------------------------------------------------------------
void dilatation3_ui8matrix_ilu3_elu2_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// -------------------------------------------------------------------------------------------
{
    int i;
    int r = ((i1 - 1) % 2);
    for (i = i0; i <= ((i1 - 1) - r); i += 2)
    {
        line_dilatation3_ui8matrix_ilu3_elu2_red(X, i, j0, j1, Y);
    }

    i = ((i1 - 1) - r);
    switch (r + 1)
    {
    case 2:
        i++;
        line_dilatation3_ui8matrix_ilu3_red(X, i, j0, j1, Y);
    case 1:
        i++;
        line_dilatation3_ui8matrix_ilu3_red(X, i, j0, j1, Y);

    default:
        break;
    }
}
// --------------------------------------------------------------------------------------------------
void dilatation3_ui8matrix_ilu3_elu2_red_factor(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------------------------------
{
    int i;
    int r = ((i1 - 1) % 2);
    for (i = i0; i <= ((i1 - 1) - r); i += 2)
    {
        line_dilatation3_ui8matrix_ilu3_elu2_red_factor(X, i + 0, j0, j1, Y);
    }

    i = ((i1 - 1) - r);
    switch (r + 1)
    {
    case 2:
        i++;
        line_dilatation3_ui8matrix_ilu3_red(X, i, j0, j1, Y);
    case 1:
        i++;
        line_dilatation3_ui8matrix_ilu3_red(X, i, j0, j1, Y);

    default:
        break;
    }
}
