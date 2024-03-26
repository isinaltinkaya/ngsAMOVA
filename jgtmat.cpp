#include "jgtmat.h"

/// @brief jgtmat_init - initialize a joint genotype matrix data structure
/// @param n - number of matrices
/// @return jgtmat_t* - pointer to the newly allocated and initialized jgtmat
jgtmat_t* jgtmat_init(const size_t n) {
    jgtmat_t* jgtmat = NULL;
    jgtmat = (jgtmat_t*)malloc(sizeof(jgtmat_t));
    ASSERT(jgtmat != NULL);

    jgtmat->n = n;
    jgtmat->pm = NULL;
    jgtmat->m = NULL;
    jgtmat->drop = NULL;

    if (ARG_DOJGTM_3GT == args->doJGTM) {
        jgtmat->size = 9;
    } else if (ARG_DOJGTM_10GT == args->doJGTM) {
        jgtmat->size = 100;
    } else {
        NEVER;
    }

    if (args->doEM) {
        jgtmat->pm = (double**)malloc(jgtmat->n * sizeof(double*));
        ASSERT(jgtmat->pm != NULL);
        for (size_t i = 0;i < jgtmat->n;++i) {
            jgtmat->pm[i] = NULL;
            jgtmat->pm[i] = (double*)malloc(jgtmat->size * sizeof(double));
            ASSERT(jgtmat->pm[i] != NULL);
            if (9 == jgtmat->size) {
                // use flat prior
                jgtmat->pm[i][0] = FRAC_1_9;
                jgtmat->pm[i][1] = FRAC_1_9;
                jgtmat->pm[i][2] = FRAC_1_9;
                jgtmat->pm[i][3] = FRAC_1_9;
                jgtmat->pm[i][4] = FRAC_1_9;
                jgtmat->pm[i][5] = FRAC_1_9;
                jgtmat->pm[i][6] = FRAC_1_9;
                jgtmat->pm[i][7] = FRAC_1_9;
                jgtmat->pm[i][8] = FRAC_1_9;
            } else if (100 == jgtmat->size) {
                // use flat prior
                for (size_t j = 0;j < jgtmat->size;++j) {
                    jgtmat->pm[i][j] = 0.01;
                }
            }
        }
    } else {
        jgtmat->m = (uint64_t**)malloc(jgtmat->n * sizeof(uint64_t*));
        ASSERT(jgtmat->m != NULL);
        for (size_t i = 0;i < jgtmat->n;++i) {
            jgtmat->m[i] = NULL;
            jgtmat->m[i] = (uint64_t*)malloc(jgtmat->size * sizeof(uint64_t));
            ASSERT(jgtmat->m[i] != NULL);
            for (size_t j = 0;j < jgtmat->size;++j) {
                jgtmat->m[i][j] = 0;
            }

        }
    }

    if (args->drop_pairs) {
        jgtmat->drop = (bool*)malloc(jgtmat->n * sizeof(bool));
        ASSERT(jgtmat->drop != NULL);
        for (size_t i = 0;i < jgtmat->n;++i) {
            jgtmat->drop[i] = false;
        }
    }

    return(jgtmat);
}

void jgtmat_destroy(jgtmat_t* jgtmat) {
    if (jgtmat->m != NULL) {
        for (size_t i = 0;i < jgtmat->n;++i) {
            FREE(jgtmat->m[i]);
        }
        FREE(jgtmat->m);
    }
    if (jgtmat->pm != NULL) {
        for (size_t i = 0;i < jgtmat->n;++i) {
            FREE(jgtmat->pm[i]);
        }
        FREE(jgtmat->pm);
    }
    if (jgtmat->drop != NULL) {
        FREE(jgtmat->drop);
    }
    FREE(jgtmat);
    return;
}


/// @brief jgtmat_print - print a joint genotype matrix data structure
/// @param jgtmat - pointer to the jgtmat to be printed
/// @param kbuf - pointer to a kstring_t to which the output will be written
void jgtmat_print(jgtmat_t* jgtmat,kstring_t* kbuf){
    //TODO also print info about number of sites shared, for both gt and gl analyses
}
