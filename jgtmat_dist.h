#ifndef __JGTMAT_DIST_H__
#define __JGTMAT_DIST_H__
/// jgtmat_dist.h: calculate distances_from_jgtmat from jgtmat

#include "dataStructs.h"
//todo instead of datastructs, use:
// jgtmat.h
// dmat.h
// etc main data types in sep hdrs


// #include <math.h>
// #include <stdio.h>
// #include <limits>

// typedef struct jgtmat_t jgtmat_t;


// JGTMAT9
// calculations based on 3x3 matrix
// 		of pairwise genotype categories
// -------------------------------------
//
// 		MM	Mm	mm
// MM	A	D	G
// Mm	B	E	H
// mm	C	F	I
//
// linearized version:
// 		A D G B E H C F I
//      0 1 2 3 4 5 6 7 8
//
//
#define JGTMAT_LINEAR_IDXOF_A 0
#define JGTMAT_LINEAR_IDXOF_B 3
#define JGTMAT_LINEAR_IDXOF_C 6
#define JGTMAT_LINEAR_IDXOF_D 1
#define JGTMAT_LINEAR_IDXOF_E 4
#define JGTMAT_LINEAR_IDXOF_F 7
#define JGTMAT_LINEAR_IDXOF_G 2
#define JGTMAT_LINEAR_IDXOF_H 5
#define JGTMAT_LINEAR_IDXOF_I 8

inline void dmat_get_distances_from_jgtmat(jgtmat_t* jgtmat, dmat_t* dmat) {

    double** pm = jgtmat->pm;

    double* dm = NULL;
    if (dmat->n == 1) {
        dm = dmat->matrix[0];
    } else {
        NEVER; //TODO add multi matrix handling? either add matrix idx to use or loop through all n matrices in here
    }


    uint32_t transform = dmat->transform;

    double x, A, B, D, E, F, H, I;
    size_t p; // pair idx
    double* ppm = NULL; // pair pm
    for (p = 0; p < dmat->size; ++p) {

        x = 0.0;
        ppm = pm[p];
        A = ppm[JGTMAT_LINEAR_IDXOF_A];
        I = ppm[JGTMAT_LINEAR_IDXOF_I];
        B = ppm[JGTMAT_LINEAR_IDXOF_B];
        D = ppm[JGTMAT_LINEAR_IDXOF_D];
        E = ppm[JGTMAT_LINEAR_IDXOF_E];
        F = ppm[JGTMAT_LINEAR_IDXOF_F];
        H = ppm[JGTMAT_LINEAR_IDXOF_H];

        x = 1.0 - (A + I + ((B + D + E + F + H) / 2.0));

        if (DMAT_TRANSFORM_NONE != transform) {
            if (DMAT_TRANSFORM_SQUARE == transform) {
                x *= x;
            } else {
                NEVER;//TODO investigate other transformation methods
            }
        }

        dm[p] = x;
    }
    return;
}






#endif // __JGTMAT_DIST_H__
