/* ------------- */
/* --- swp.h --- */
/* ------------- */

/*
* Copyright (c) 2020 - 2021, Lionel Lacassagne, All rights reserved
* Sorbonne University, LIP6, CNRS
*/

#ifndef __SWP_TEST_H__
#define __SWP_TEST_H__

#ifdef __cplusplus
#pragma message ("C++")
extern "C" {
#endif

// bit to Register (8-bit version)
#define B8R(b7,b6,b5,b4,b3,b2,b1,b0) ((b7<<7) | (b6<<6) | (b5<<5) | (b4<<4) | (b3<<3) | (b2<<2) | (b1<<1)  | (b0<<0))

// register to bits (8-bit version)
#define R8B(x, b7,b6,b5,b4,b3,b2,b1,b0) b0=x&1; b1=(x>>1)&1; b2=(x>>2)&1; b3=(x>>3)&1; b4=(x>>4)&1; b5=(x>>5)&1; b6=(x>>6)&1; b7=(x>>7)&1
    
// bit to memory (8-bit version)
#define B8M(b0,b1,b2,b3,b4,b5,b6,b7) B8R(b7,b6,b5,b4,b3,b2,b1,b0)
//#define B8M(b0,b1,b2,b3,b4,b5,b6,b7) B8R(b0,b1,b2,b3,b4,b5,b6,b7)

// memory to bits (8-bit version)
#define M8B(x, b0,b1,b2,b3,b4,b5,b6,b7) R8B(x, b7,b6,b5,b4,b3,b2,b1,b0)
//#define M8B(x, b0,b1,b2,b3,b4,b5,b6,b7) R8B(x, b0,b1,b2,b3,b4,b5,b6,b7)

// bit to memory (16-bit version)
#define B16R(   bf,be,bd,bc,bb,ba,b9,b8, b7,b6,b5,b4,b3,b2,b1,b0) (B8R(bf,be,bd,bc,bb,ba,b9,b8) << 8) | B8R(b7,b6,b5,b4,b3,b2,b1,b0)

// register to bits (16-bit version)
#define R16B(x, bf,be,bd,bc,bb,ba,b9,b8, b7,b6,b5,b4,b3,b2,b1,b0) R8B( ((x >> 8)& 0xff), bf,be,bd,bc,bb,ba,b9,b8); R8B( (x & 0xff), b7,b6,b5,b4,b3,b2,b1,b0)

// bit to memory (16-bit version)
#define B16M(b0,b1,b2,b3,b4,b5,b6,b7, b8,b9,ba,bb,bc,bd,be,bf)  (B8M(b0,b1,b2,b3,b4,b5,b6,b7) | B8M(b8,b9,ba,bb,bc,bd,be,bf) << 8)
 
// memory to bits (16-bit version)
#define M16B(x, b0,b1,b2,b3,b4,b5,b6,b7, b8,b9,ba,bb,bc,bd,be,bf) M8B((x & 0xff), b0,b1,b2,b3,b4,b5,b6,b7); M8B(((x >> 8) & 0xff), b8,b9,ba,bb,bc,bd,be,bf)


// calcul d'indices (ne pas utiliser car ne couvre pas tous les cas)
void B2b(int i0, int i1, int card, int* b0, int* b1);
void b2B(int b0, int b1, int card, int* m0, int* m1);

// versions specialisees des fonctions precedentes
void Btob8 (int i0, int j1, int* b0, int* b1);
void Btob16(int i0, int j1, int* b0, int* b1);
void Btob32(int i0, int j1, int* b0, int* b1);
void Btob64(int i0, int j1, int* b0, int* b1);

void b8toB (int b0, int b1, int* m0, int* m1);
void b16toB(int b0, int b1, int* m0, int* m1);
void b32toB(int b0, int b1, int* m0, int* m1);
void b64toB(int b0, int b1, int* m0, int* m1);

// affichage 8 bits en registe / en memoire (version hard-coded sans boucle)
void printfR8(uint8_t x8r, char* name);
void printfM8(uint8_t x8m, char* name);

// affichage 16 bits en registe / en memoire (version hard-coded sans boucle)
void printfR16(uint16_t x16r, char* name);
void printfM16(uint16_t x16m, char* name);

// affichage 32 bit en registe / en memoire
// version avec boucle
void printfR32(uint32_t x, char* name);
void printfM32(uint32_t x, char* name);

// affichage 64 bit en registe / en memoire
// version avec boucle (fonctionne avec 8 /16 /32 / 64 bits
void printfR(uint64_t x, int n, char* name);
void printfM(uint64_t x, int n, char* name);

uint32_t str2Reg(char* str);
uint32_t str2Mem(char* str);

void set_ui8vector_str(uint8 *v, int i0, int i1, char *str);

int ishexa(char x);
// extension of isdigit to hexa number

int atob(char x); // atoi for bit
int atox(char x); // atoi for hexa (return 0 otherwise)

void set_ui8vector_str1 (uint8* v, int i0, int i1, char *str);
void set_ui8vector_str16(uint8* v, int i0, int i1, char *str);

void pack_ui4vector (uint8 *X1, int src_width, uint8  *Y4 );
void pack_ui8vector (uint8 *X1, int src_width, uint8  *Y8 );
void pack_ui16vector(uint8 *X1, int src_width, uint16 *Y16);
void pack_ui32vector(uint8 *X1, int src_width, uint32 *Y32);
void pack_ui64vector(uint8 *X1, int src_width, uint64 *Y64);

void unpack_ui4vector (uint8  *X4 , int src_width, uint8* Y1);
void unpack_ui8vector (uint8  *X8 , int src_width, uint8* Y1);
void unpack_ui16vector(uint16 *X16, int src_width, uint8* Y1);
void unpack_ui32vector(uint32 *X32, int src_width, uint8* Y1);
void unpack_ui64vector(uint64 *X64, int src_width, uint8* Y1);

void pack_ui4matrix (uint8  **X1, int src_height, int src_width, uint8  **Y4 );
void pack_ui8matrix (uint8  **X1, int src_height, int src_width, uint8  **Y8 );
void pack_ui16matrix(uint8  **X1, int src_height, int src_width, uint16 **Y16);
void pack_ui32matrix(uint8  **X1, int src_height, int src_width, uint32 **Y32);
void pack_ui64matrix(uint8  **X1, int src_height, int src_width, uint64 **Y64);

void unpack_ui4matrix (uint8  **X4 , int src_height, int src_width, uint8** Y1);
void unpack_ui8matrix (uint8  **X8 , int src_height, int src_width, uint8** Y1);
void unpack_ui16matrix(uint16 **X16, int src_height, int src_width, uint8** Y1);
void unpack_ui32matrix(uint32 **X32, int src_height, int src_width, uint8** Y1);
void unpack_ui64matrix(uint64 **X64, int src_height, int src_width, uint8** Y1);

void displayR_ui4vector (uint8  *v, int i0, int i1, char *name);
void displayR_ui8vector (uint8  *v, int i0, int i1, char *name);

void displayM_ui4vector (uint8  *v, int i0, int i1, char *name);
void displayM_ui8vector (uint8  *v, int i0, int i1, char *name);
void displayM_ui16vector(uint16 *v, int i0, int i1, char *name);
void displayM_ui32vector(uint32 *v, int i0, int i1, char *name);
void displayM_ui64vector(uint64 *v, int i0, int i1, char *name);

void displayR_ui4matrix (uint8  **m, int i0, int i1, int j0, int j1, char *name);
void displayR_ui8matrix (uint8  **m, int i0, int i1, int j0, int j1, char *name);

void displayM_ui4matrix (uint8  **m, int i0, int i1, int j0, int j1, char *name);
void displayM_ui8matrix (uint8  **m, int i0, int i1, int j0, int j1, char *name);
void displayM_ui16matrix(uint16 **m, int i0, int i1, int j0, int j1, char *name);
void displayM_ui32matrix(uint32 **m, int i0, int i1, int j0, int j1, char *name);
void displayM_ui64matrix(uint64 **m, int i0, int i1, int j0, int j1, char *name);

#ifdef __cplusplus
}
#endif

#endif // __SWP_H__
