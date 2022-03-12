/* -------------------- */
/* --- sigmadelta.c --- */
/* -------------------- */

/*
 * Copyright (c) 2004 - 2013, Lionel Lacassagne, All rights reserved
 * University of Paris Sud, Laboratoire de Recherche en Informatique 
 * Creation: 2004-05-18 :
 * Creation: 2021-01-06 : version line pour pipeline
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "nrtype.h"
#include "nrdef.h"
#include "nrutil.h"

#include "sigmadelta.h"
#include "morpho.h"

/* --------------- */
int findPower2(int x)
/* --------------- */
{
    int p = 0;
    if(!x) return 0;
    
    while(!(x & 1)) {
        x = x / 2;
        p = p + 1;
    }
    return p;
}
// -----------------------------------------------------------------------------------------
void SigmaDelta_Step0_line(uint8 *I, uint8 *M, uint8 *O, uint8 *V, uint8 *E, int j0, int j1)
// -----------------------------------------------------------------------------------------
{
    /*

        Initialisation de SD pour t = 0 version ligne

    */

    // debug
    setborder1(j0, j1);

    int j;
    uint8 m, i, v;
    for (j = j0; j < j1; j++)
    {
        // Load
        i = load1(I, j);

        m = i;
        v = SD_VMIN;

        // Store
        store1(M, j, m);
        store1(V, j, v);
    }
}
// ------------------------------------------------------------------------------------------------
void SigmaDelta_1Step_line(uint8 *I, uint8 *M, uint8 *O, uint8 *V, uint8 *E, int k, int j0, int j1)
// ------------------------------------------------------------------------------------------------
{
    /*

        Debut de SD à partir de t = 1

    */

    // debug
    setborder1(j0, j1);
    
    int j;
    uint8 i, m, o, v, e;
    
    for (j = j0; j < j1; j++)
    {
        // load
        i = load1(I, j);
        m = load1(M, j);
        o = load1(O, j);
        v = load1(V, j);
        e = load1(E, j);

        // STEP 1 : Estimation Mt 
        if (m < i)
        {
            m = m + 1;
        }
        else if (m > i)
        {
            m = m - 1;
        }
        
        // STEP 2 : Ot computation 
        o = abs(m - i);

        // STEP 3 : Vt update and clamping 
        if (v < (k * o))
        {
            v = v + 1;
        }
        else if (v > (k * o))
        {
            v = v - 1;
        }

        v = ((v < SD_VMAX)? v: SD_VMAX);
        v = ((v > SD_VMIN)? v: SD_VMIN);

        // STEP 4 : Et estimation 
        if (o < v)
        {
            e = 0;
        }
        else
        {
            e = 1;
        }

        // store
        store1(M, j, m);
        store1(O, j, o);
        store1(V, j, v);
        store1(E, j, e);
    }
}
// ---------------------------------------------------------------------------------------------------------
void SigmaDelta_Step0(uint8 **I, uint8 **M, uint8 **O, uint8 **V, uint8 **E, int i0, int i1, int j0, int j1)
// ---------------------------------------------------------------------------------------------------------
{
    for(int i=i0; i<=i1; i++) {
        SigmaDelta_Step0_line(I[i], M[i], O[i], V[i], E[i], j0, j1);
    }
}
// ----------------------------------------------------------------------------------------------------------------
void SigmaDelta_1Step(uint8 **I, uint8 **M, uint8 **O, uint8 **V, uint8 **E, int k, int i0, int i1, int j0, int j1)
// ----------------------------------------------------------------------------------------------------------------
{
    for(int i=i0; i<=i1; i++) {
        SigmaDelta_1Step_line(I[i], M[i], O[i], V[i], E[i], k, j0, j1);
    }
}
