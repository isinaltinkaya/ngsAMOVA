#include "jgtmat.h"

size_t jgtmat_estimate_max_mem_use(const size_t n, const size_t size, const uint8_t n_gc) {

    LOG("Estimating maximum memory requirement for jgtmat with number of matrices: %ld, size: %ld, number of genotype categories: %d", n, size, n_gc);

    size_t mem = 0;

    // -> jgtmat_t jgtmat
    mem += sizeof(jgtmat_t);

    // -> double jgtmat->pm[n][size][n_gc]
    // -> uint64_t jgtmat->m[n][size][n_gc]
    bool has_pm = false;
    bool has_m = false;
    if (PROGRAM_WILL_USE_BCF_FMT_GL) {
        has_pm = true;
    }
    if (PROGRAM_WILL_USE_BCF_FMT_GT) {
        has_m = true;
    }
    if (has_pm) {
        mem += GET_ARR_SIZE_BYTES_3D(double, n, size, n_gc);
    }
    if (has_m) {
        mem += GET_ARR_SIZE_BYTES_3D(uint64_t, n, size, n_gc);
    }

    // -> uint64_t jgtmat->snsites[n][size]
    mem += GET_ARR_SIZE_BYTES_2D(uint64_t, n, size);

    // -> bool jgtmat->drop[n][size]
    mem += GET_ARR_SIZE_BYTES_2D(bool, n, size);

    return (mem);
}


jgtmat_t* jgtmat_init(const size_t n, const size_t size, const uint8_t n_gc) {

    jgtmat_t* jgtmat = NULL;
    jgtmat = (jgtmat_t*)malloc(sizeof(jgtmat_t));
    ASSERT(jgtmat != NULL);

    // -> init
    jgtmat->n = 0;
    jgtmat->size = 0;
    jgtmat->n_gc = 0;
    jgtmat->pm = NULL;
    jgtmat->m = NULL;
    jgtmat->snsites = NULL;
    jgtmat->drop = NULL;

    // -> set 

    jgtmat->n = n;
    jgtmat->size = size;
    jgtmat->n_gc = n_gc;

    if (args->doEM) {

        jgtmat->pm = (double***)malloc(jgtmat->n * sizeof(double**));
        ASSERT(jgtmat->pm != NULL);

        for (size_t i = 0;i < jgtmat->n;++i) {
            jgtmat->pm[i] = NULL;
            jgtmat->pm[i] = (double**)malloc(jgtmat->size * sizeof(double*));
            ASSERT(jgtmat->pm[i] != NULL);

            // if (9 == jgtmat->size) {
            //     // use flat prior
            //     jgtmat->pm[i][0] = FRAC_1_9;
            //     jgtmat->pm[i][1] = FRAC_1_9;
            //     jgtmat->pm[i][2] = FRAC_1_9;
            //     jgtmat->pm[i][3] = FRAC_1_9;
            //     jgtmat->pm[i][4] = FRAC_1_9;
            //     jgtmat->pm[i][5] = FRAC_1_9;
            //     jgtmat->pm[i][6] = FRAC_1_9;
            //     jgtmat->pm[i][7] = FRAC_1_9;
            //     jgtmat->pm[i][8] = FRAC_1_9;
            // } else if (100 == jgtmat->size) {
            //     // use flat prior
            //     for (size_t j = 0;j < jgtmat->size;++j) {
            //         jgtmat->pm[i][j] = 0.01;
            //     }
            // }

            double flat = 0.0;
            if (9 == jgtmat->n_gc) {
                flat = FRAC_1_9;
            } else if (100 == jgtmat->n_gc) {
                flat = 0.01;
            }

            for (size_t j = 0;j < jgtmat->size;++j) {
                jgtmat->pm[i][j] = NULL;
                jgtmat->pm[i][j] = (double*)malloc(jgtmat->n_gc * sizeof(double));
                ASSERT(jgtmat->pm[i][j] != NULL);
                for (size_t k = 0;k < jgtmat->n_gc;++k) {
                    jgtmat->pm[i][j][k] = flat;
                }
            }
        }
    } else {
        jgtmat->m = (uint64_t***)malloc(jgtmat->n * sizeof(uint64_t**));
        ASSERT(jgtmat->m != NULL);
        for (size_t i = 0;i < jgtmat->n;++i) {
            jgtmat->m[i] = NULL;
            jgtmat->m[i] = (uint64_t**)malloc(jgtmat->size * sizeof(uint64_t*));
            ASSERT(jgtmat->m[i] != NULL);
            for (size_t j = 0;j < jgtmat->size;++j) {
                jgtmat->m[i][j] = NULL;
                jgtmat->m[i][j] = (uint64_t*)malloc(jgtmat->n_gc * sizeof(uint64_t));
                ASSERT(jgtmat->m[i][j] != NULL);
                for (size_t k = 0;k < jgtmat->n_gc;++k) {
                    jgtmat->m[i][j][k] = 0;
                }
            }

        }
    }

    jgtmat->snsites = (uint64_t**)malloc(jgtmat->n * sizeof(uint64_t*));
    ASSERT(jgtmat->snsites != NULL);
    for (size_t i = 0;i < jgtmat->n;++i) {
        jgtmat->snsites[i] = NULL;
        jgtmat->snsites[i] = (uint64_t*)malloc(jgtmat->size * sizeof(uint64_t));
        ASSERT(jgtmat->snsites[i] != NULL);
        for (size_t j = 0;j < jgtmat->size;++j) {
            jgtmat->snsites[i][j] = 0;
        }
    }


    if (args->allow_mispairs) {
        jgtmat->drop = (bool**)malloc(jgtmat->n * sizeof(bool*));
        ASSERT(jgtmat->drop != NULL);
        for (size_t i = 0;i < jgtmat->n;++i) {
            jgtmat->drop[i] = NULL;
            jgtmat->drop[i] = (bool*)malloc(jgtmat->size * sizeof(bool));
            ASSERT(jgtmat->drop[i] != NULL);
            for (size_t j = 0;j < jgtmat->size;++j) {
                jgtmat->drop[i][j] = false;
            }
        }
    }

    return(jgtmat);
}

void jgtmat_destroy(jgtmat_t* jgtmat) {
    if (jgtmat->pm != NULL) {
        for (size_t i = 0;i < jgtmat->n;++i) {
            for (size_t j = 0;j < jgtmat->size;++j) {
                FREE(jgtmat->pm[i][j]);
            }
            FREE(jgtmat->pm[i]);
        }
        FREE(jgtmat->pm);
    }
    if (jgtmat->m != NULL) {
        for (size_t i = 0;i < jgtmat->n;++i) {
            for (size_t j = 0;j < jgtmat->size;++j) {
                FREE(jgtmat->m[i][j]);
            }
            FREE(jgtmat->m[i]);
        }
        FREE(jgtmat->m);
    }
    for (size_t i = 0;i < jgtmat->n;++i) {
        FREE(jgtmat->snsites[i]);
    }
    FREE(jgtmat->snsites);
    if (jgtmat->drop != NULL) {
        for (size_t i = 0;i < jgtmat->n;++i) {
            FREE(jgtmat->drop[i]);
        }
        FREE(jgtmat->drop);
    }
    FREE(jgtmat);
    return;
}


void jgtmat_print(jgtmat_t* jgtmat, outfile_t* outfile) {
    LOG("(--print-jgtmat) Printing joint genotype matrix to file: %s\n", outfile->fn);

    kstring_t* kbuf = &outfile->kbuf;

    // line 0: number of matrices
    ksprintf(kbuf, "%zu\n", jgtmat->n);

    // line 1: size
    ksprintf(kbuf, "%zu\n", jgtmat->size);

    // line 2: number of genotype categories
    ksprintf(kbuf, "%u\n", jgtmat->n_gc);

    // line 3: the type of the values in the matrix (0: counts, 1: probabilities)
    uint8_t type = 2; // init
    if (jgtmat->pm != NULL && jgtmat->m != NULL) {
        NEVER;
    }
    if (jgtmat->pm != NULL) {
        type = 1;
    }
    if (jgtmat->m != NULL) {
        type = 0;
    }
    ASSERT(type < 2);
    ksprintf(kbuf, "%u\n", type);

    // line 4: does the matrix contain drop data? (0: no, 1: yes)
    uint8_t hasdrop = (jgtmat->drop != NULL) ? 1 : 0;
    ksprintf(kbuf, "%u\n", hasdrop);

    // lines [5,X): drop info (if any)
    // if line 4 == 0, X=5 (no such lines)
    if (args->allow_mispairs) {
        for (size_t i = 0;i < jgtmat->n;++i) {
            for (size_t j = 0;j < jgtmat->size;++j) {
                ksprintf(kbuf, "%d\n", jgtmat->drop[i][j]);
            }
        }
    }

    // lines [X,Y): number of sites shared by the individuals for each individual pair
    // Y = X + (n * size)
    for (size_t i = 0;i < jgtmat->n;++i) {
        for (size_t j = 0;j < jgtmat->size;++j) {
            ksprintf(kbuf, "%zu\n", jgtmat->snsites[i][j]);
        }
    }

    // lines [Y,Z): matrix data
    // Z = Y + (n * size * n_gc)
    if (jgtmat->pm != NULL) {
        for (size_t i = 0;i < jgtmat->n;++i) {
            for (size_t j = 0;j < jgtmat->size;++j) {
                for (size_t k = 0;k < jgtmat->n_gc;++k) {
                    ksprintf(kbuf, "%.17g\n", jgtmat->pm[i][j][k]);
                }
            }
        }
    } else if (jgtmat->m != NULL) {
        for (size_t i = 0;i < jgtmat->n;++i) {
            for (size_t j = 0;j < jgtmat->size;++j) {
                for (size_t k = 0;k < jgtmat->n_gc;++k) {
                    ksprintf(kbuf, "%ld\n", jgtmat->m[i][j][k]);
                }
            }
        }
    } else {
        NEVER;
    }

    return;
}


jgtmat_t* jgtmat_read(void) {
    //TODO
    NEVER;
    jgtmat_t* jgtmat = NULL;
    return(jgtmat);
}

