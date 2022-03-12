/* ---------------- */
/* --- morpho.h --- */
/* ---------------- */

/*
* Copyright (c) 2020 - 2021, Lionel Lacassagne, All rights reserved
* Sorbonne University, LIP6, CNRS
*/

#ifndef __MORPHO_H__
#define __MORPHO_H__

#ifdef __cplusplus
#ifdef PRAGMA_VERBOSE
#pragma message ("C++")
#endif
extern "C" {
#endif
    // #define DEBUG
    // #define ENABLE_CONTROL
    #ifdef DEBUG
        #ifndef ENABLE_CONTROL
            // display
            #define idisp(x) printf("%s = %d\n", #x, x)
            #define idisp9(a0, a1, a2, b0, b1, b2, c0, c1, c2) printf("-----\n%s = %d %s = %d %s = %d\n%s = %d %s = %d %s = %d\n%s = %d %s = %d %s = %d\n-------\n", #a0, a0, #a1, a1, #a2, a2, #b0, b0, #b1, b1, #b2, b2, #c0, c0, #c1, c1, #c2, c2)
            #define disp(x) printf("%s = %d\n", #x, x)
            #define VERBOSE(X) X
            #define PUTS(str) puts(str)
            #define CR putchar('\n');

            // LOAD STORE
            #define setborder1(i0, i1)
            #define load1(X, i) X[i]; //printf("Load : %s[%d] = %d\n", #X, i, X[i])
            #define store1(Y, i, y) Y[i]=y; //printf("Store : %s[%d] = %d\n", #Y, i, y)

        #endif // !ENABLE_CONTROL

        #ifdef ENABLE_CONTROL
            #define idisp(x) printf("%s = %3d\n", #x, x)
            #define idisp9(a0, a1, a2, b0, b1, b2, c0, c1, c2) printf("%s = %d %s = %d %s = %d\n%s = %d %s = %d %s = %d\n%s = %d %s = %d %s = %d\n", #a0, a0, #a1, a1, #a2, a2, #b0, b0, #b1, b1, #b2, b2, #c0, c0, #c1, c1, #c2, c2)
            #define disp(x) printf("%s = %5.0f\n", #x, x)
            #define VERBOSE(X) X
            #define PUTS(str) puts(str)
            #define CR putchar('\n');

            
            #define setborder1(i0, i1) int dbg_i0 = i0-1; int dbg_i1 = i1+1
            #define load1(X, i) X[i];\
                if((i < dbg_i0) || (i > dbg_i1)){\
                    printf("\n\e[31mFATAL ERROR\e[0m line %d in file %s :\nLoad %s[%s] => interval is [%d, %d] but accessing index %d\n", __LINE__, __FILE__, #X, #i, dbg_i0, dbg_i1, i);\
                    exit(EXIT_FAILURE);\
                }//else{printf("Load : %s[%d] = %5.0f\n", #X, i, X[i]);}
            #define store1(Y, i, y) Y[i]=y;\
                if((i < dbg_i0) || (i > dbg_i1)){\
                    printf("\n\e[31mFATAL ERROR\e[0m line %d in file %s :\nStore %s[%s] => interval is [%d, %d] but accessing index %d\n", __LINE__, __FILE__, #Y, #i, dbg_i0, dbg_i1, i);\
                    exit(EXIT_FAILURE);\
                }//else{printf("Store : %s[%d] = %5.0f\n", #Y, i, y);}

        #endif // ENABLE_CONTROL
    #else
        // release
        #define idisp(x)
        #define idisp9(a0, a1, a2, b0, b1, b2, c0, c1, c2)
        #define disp(x)

        #define VERBOSE(X)
        #define PUTS(str)
        #define CR


        #define setborder1(i0, i1)

        #define load1(X, i) X[i]
        #define store1(Y, i, y) Y[i]=y

        #define load2(X, i, j) X[i][j]
        #define store2(Y, i, j, y) Y[i][j]=y
    #endif

#ifdef DEBUG
    #define AND +
    #define OR +
#else
    #define AND &
    #define OR |
#endif

// MAX

#define OP1MAX(x, y) ((x) > (y) ? (x) : (y))
#define OP2MAX(x, y) ((x) OR (y))


/*
http://graphics.stanford.edu/~seander/bithacks.html#IntegerMinOrMax

Operation sur bit pour min et max si branchement couteux
Operationt ternaire => branchement
(<, >) peut dans certaines machines utilisé un branchement
*/
#define OP3MAX(x, y) (x ^ ((x ^ y) & -(x < y))) // clairement plus lent pour nous

#define MAX(x, y) OP2MAX(x, y)
#define MAX3(x, y, z) MAX(MAX(x, y), z)
#define MAX4(x0, x1, x2, x3) MAX(MAX(x0, x1), MAX(x2, x3))
#define MAX6(x0, x1, x2, x3, x4, x5) MAX(MAX3(x0, x1, x2), MAX3(x3, x4,x5))
#define MAX9(x0, y0, z0, x1, y1, z1, x2, y2, z2) MAX3(MAX3(x0, y0, z0), MAX3(x1, y1, z1), MAX3(x2, y2, z2))

// MIN
#define OP1MIN(x, y) ((x) < (y) ? (x) : (y))
#define OP2MIN(x, y) ((x) AND (y))

#define MIN(x, y) OP2MIN(x, y)
#define MIN3(x, y, z) MIN(MIN(x, y), z)
#define MIN4(x0, x1, x2, X3) MIN(MIN(x0, x1), MIN(x2, x3))
#define MIN6(x0, x1, x2, x3, x4, x5) MIN(MIN3(x0, x1, x2), MIN3(x3, x4, x5))
#define MIN9(x0, y0, z0, x1, y1, z1, x2, y2, z2) MIN3(MIN3(x0, y0, z0), MIN3(x1, y1, z1), MIN3(x2, y2, z2))

// SWP
#define K_SHIFT 1

#define LEFT_8(a,b)   ((a >> (8 -K_SHIFT))  | (b << K_SHIFT     ))
#define RIGHT_8(b,c)  ((b >> K_SHIFT     )  | (c << (8 -K_SHIFT)))

#define LEFT_32(a,b)  ((a >> (32-K_SHIFT))  | (b << K_SHIFT     ))
#define RIGHT_32(b,c) ((b >> K_SHIFT     )  | (c << (32-K_SHIFT)))

#define LEFT_64(a,b)  ((a >> (64-K_SHIFT))  | (b << K_SHIFT     ))
#define RIGHT_64(b,c) ((b >> K_SHIFT     )  | (c << (64-K_SHIFT)))

#ifdef __cplusplus
}
#endif

#endif // __MORPHO_H__
