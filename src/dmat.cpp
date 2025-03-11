#include <math.h>

#include "dmat.h"
#include "dataStructs.h"
#include "jgtmat.h"
#include "metadata.h"
#include "io.h"
#include "euclid.h"

bool dmat_has_identity_of_indiscernibles(dmat_t* dmat) {
    ASSERT(dmat->type == DMAT_TYPE_FULL);
    const size_t nelems = sqrt(dmat->size);
    size_t idx = 0;
    for (size_t matrix_i = 0;matrix_i < dmat->n;++matrix_i) {
        idx = 0;
        for (size_t i = 0; i < nelems;++i) {
            for (size_t j = 0; j < nelems;++j) {
                if (i == j) {
                    if (dmat->matrix[matrix_i][idx] != 0.0) {
                        return(false);
                    }
                }
                idx++;
            }
        }
    }
    return(true);
}

bool dmat_is_euclidean(dmat_t* dmat, double tol) {
    if (dmat->type != DMAT_TYPE_FULL) {
        ERROR("Only full distance matrices are supported for Euclidean check.");
    }
    if (dmat->n != 1) {
        WARN("Multiple distance matrices are provided to dmat_is_euclidean. Program will only check the first distance matrix for Euclidean property, and assume all matrices have the same property.");
    }
    size_t n = sqrt(dmat->size);
    bool isEuclidean = matrix_is_euclidean(dmat->matrix[0], n, tol);
    return(isEuclidean);
}

dmat_t* new_dmat_type_lted_to_full(dmat_t* lted_dmat) {
    // convert a lted matrix to full matrix
    // -> first figure out the number of elements in the full matrix
    const size_t lted_dmat_n = lted_dmat->n;
    const size_t lted_dmat_size = lted_dmat->size;

    const size_t lted_dmat_nelems = 1 + sqrt(1 + (8 * lted_dmat_size)) / 2;
    const size_t full_dmat_size = lted_dmat_nelems * lted_dmat_nelems;
    dmat_t* full_dmat = dmat_init(lted_dmat_n, lted_dmat_nelems, DMAT_TYPE_FULL, lted_dmat->method, lted_dmat->transform, lted_dmat->names, lted_dmat->names_src);

    // -> full_dmat->matrix all elements are already initted to 0 in dmat_init
    // so no need to set diagonal elements to 0

    // -> now copy the elements from lted to full
    size_t full_idx = 0;
    for (size_t matrix_i = 0;matrix_i < lted_dmat_n;++matrix_i) {
        const bool fill_drop = (lted_dmat->drop == NULL) ? false : lted_dmat->drop[matrix_i] != NULL;
        if (fill_drop) {
            full_dmat->drop[matrix_i] = (bool*)malloc(full_dmat_size * sizeof(bool));
            ASSERT(full_dmat->drop[matrix_i] != NULL);
        }
        full_idx = 0;
        for (size_t i = 0;i < lted_dmat_nelems;++i) {
            for (size_t j = 0;j < lted_dmat_nelems;++j) {
                if (i >= 1 && j < i) {
                    full_dmat->matrix[matrix_i][full_idx] = lted_dmat->matrix[matrix_i][MATRIX_GET_INDEX_LTED_IJ(i, j)];
                    if (fill_drop) {
                        full_dmat->drop[matrix_i][full_idx] = lted_dmat->drop[matrix_i][MATRIX_GET_INDEX_LTED_IJ(i, j)];
                    }
                }
                if (j >= 1 && i < j) {
                    full_dmat->matrix[matrix_i][full_idx] = lted_dmat->matrix[matrix_i][MATRIX_GET_INDEX_LTED_IJ(j, i)];
                    if (fill_drop) {
                        ASSERT(full_idx < full_dmat_size);
                        ASSERT(full_dmat->drop != NULL);
                        ASSERT(full_dmat->drop[matrix_i] != NULL);
                        full_dmat->drop[matrix_i][full_idx] = lted_dmat->drop[matrix_i][MATRIX_GET_INDEX_LTED_IJ(j, i)];
                    }
                }
                full_idx++;
            }
        }

    }
    return(full_dmat);
}

dmat_t* new_dmat_type_full_to_lted(dmat_t* full_dmat) {
    // convert a full matrix to lted matrix
    // -> first figure out the number of elements in the lted matrix
    const size_t full_dmat_n = full_dmat->n;
    const size_t full_dmat_size = full_dmat->size;
    const size_t full_dmat_nelems = sqrt(full_dmat_size);
    dmat_t* lted_dmat = dmat_init(full_dmat_n, full_dmat_nelems, DMAT_TYPE_LTED, full_dmat->method, full_dmat->transform, full_dmat->names, full_dmat->names_src);

    const bool fill_drop = lted_dmat->drop != NULL;

    // -> now copy the elements from full to lted
    size_t lted_idx = 0, full_idx = 0;
    for (size_t matrix_i = 0;matrix_i < full_dmat_n;++matrix_i) {
        lted_idx = 0;
        full_idx = 0;
        for (size_t i = 0;i < full_dmat_nelems;++i) {
            for (size_t j = 0;j < full_dmat_nelems;++j) {
                if (i >= 1 && j < i) {
                    lted_dmat->matrix[matrix_i][lted_idx] = full_dmat->matrix[matrix_i][full_idx];
                    if (fill_drop) {
                        lted_dmat->drop[matrix_i][lted_idx] = full_dmat->drop[matrix_i][full_idx];
                    }
                    lted_idx++;
                }
                full_idx++;
            }
        }
    }
    return(lted_dmat);
}

dmat_t* new_dmat_type_to_type(dmat_t* dmat, uint8_t required_type) {
    if (dmat->type == required_type) {
        return(dmat);
    }
    if (required_type == DMAT_TYPE_FULL) {
        if (dmat->type == DMAT_TYPE_LTED) {
            return(new_dmat_type_lted_to_full(dmat));
        } else {
            NEVER;
        }
    }
    if (dmat->type == DMAT_TYPE_FULL) {
        if (required_type == DMAT_TYPE_LTED) {
            return(new_dmat_type_full_to_lted(dmat));
        } else {
            NEVER;
        }
    }
    NEVER;
}

void dmat_type_to_type(dmat_t** dmat, uint8_t required_type) {
    const uint8_t dmat_type = (*dmat)->type;
    if (dmat_type == required_type) {
        return;
    }
    dmat_t* new_dmat = NULL;
    if (required_type == DMAT_TYPE_FULL) {
        if (dmat_type == DMAT_TYPE_LTED) {
            new_dmat = new_dmat_type_lted_to_full(*dmat);
        } else {
            NEVER;
        }
    }
    if (dmat_type == DMAT_TYPE_FULL) {
        if (required_type == DMAT_TYPE_LTED) {
            new_dmat = new_dmat_type_full_to_lted(*dmat);
        } else {
            NEVER;
        }
    }
    dmat_destroy(*dmat);
    *dmat = new_dmat;
    return;
}

size_t dmat_estimate_max_mem_use(const size_t n, const size_t nInd, const int allow_mispairs, const uint8_t type) {
    LOG("Estimating maximum memory requirement for dmat with number of distance matrices: %ld, number of individuals: %ld, allow_mispairs: %d, dmat type: %d", n, nInd, allow_mispairs, type);
    size_t mem = 0;

    size_t size;
    switch (type) {
    case DMAT_TYPE_LTED:
    case DMAT_TYPE_UTED:
        size = (nInd * (nInd - 1)) / 2;
        break;
    case DMAT_TYPE_LTID:
    case DMAT_TYPE_UTID:
        size = (nInd * (nInd + 1)) / 2;
        break;
    case DMAT_TYPE_FULL:
        size = nInd * nInd;
        break;
    default:
        NEVER;
    }

    // -> dmat_t dmat
    mem += sizeof(dmat_t);

    // -> strArray* dmat->names
    // mem += sizeof(strArray);
    //TODO also add the actual size of strArray after alloc
    // strArray_estimate_max_mem_use(nInds);

    // -> double dmat->matrix[n][size]
    mem += GET_ARR_SIZE_BYTES_2D(double, n, size);

    // -> bool dmat->drop[n][size]
    if (allow_mispairs != 0) {
        mem += GET_ARR_SIZE_BYTES_2D(bool, n, size);
    }

    return mem;
}

const char* get_dmat_method_str(const int method) {
    if (method < 0 || method > 9) {
        ERROR("Unrecognized distance matrix method: %d", method);
    }
    const char* dmat_method_to_str[10] = {
        "Dij",  // DMAT_METHOD_DIJ
        "Sij",  // DMAT_METHOD_SIJ
        "Fij",  // DMAT_METHOD_FIJ
        "IBS0", // DMAT_METHOD_IBS0
        "IBS1", // DMAT_METHOD_IBS1
        "IBS2", // DMAT_METHOD_IBS2
        "Kin",  // DMAT_METHOD_KIN
        "R0",   // DMAT_METHOD_R0       
        "R1",   // DMAT_METHOD_R1
        "Dxy"  // DMAT_METHOD_DXY
    };
    return (dmat_method_to_str[method]);
}


//TODO add testcase
dmat_t* new_dmat_pruned(dmat_t* dmat) {

    BEGIN_LOGSECTION_MSG("(--prune-dmat) Pruning distance matrix");

    if (dmat->n > 1) {
        ERROR("Pruning of distance matrices with more than one matrix is not possible, yet.");
    }

    // -> init
    dmat_t* pruned_dmat = NULL;
    pruned_dmat = (dmat_t*)malloc(sizeof(dmat_t));
    ASSERT(pruned_dmat != NULL);
    pruned_dmat->n = 0;
    pruned_dmat->size = 0;
    pruned_dmat->type = 0;
    pruned_dmat->transform = DMAT_INTPLUS_TRANSFORM_NONE;
    pruned_dmat->method = 0;
    pruned_dmat->names = NULL;
    pruned_dmat->names_src = 0;
    pruned_dmat->drop = NULL;

    //TODO implement dropping for multiple matrices
    // N.B. multiple matrices do not share the same drop
    const size_t ni = 0; // dmat->n - 1

    size_t idx;
    int maxDrop, maxDropIdx, totnPairsToPrune, nItems;
    nItems = dmat->names->len;


    int* nDrops = NULL;
    nDrops = (int*)malloc(nItems * sizeof(int));
    ASSERT(nDrops != NULL);

    int* rmIdx = NULL;
    rmIdx = (int*)malloc(nItems * sizeof(int));
    ASSERT(rmIdx != NULL);

    for (size_t i = 0;i < (size_t)nItems;++i) {
        nDrops[i] = 0;
        rmIdx[i] = -1;
    }

    pruned_dmat->n = dmat->n;
    // pruned_dmat->size =  // TBD
    pruned_dmat->type = dmat->type;
    pruned_dmat->transform = dmat->transform;
    pruned_dmat->method = dmat->method;
    // pruned_dmat->names = // TBD
    pruned_dmat->names_src = DMAT_NAMES_SRC_PRIVATE;
    pruned_dmat->matrix = (double**)malloc(pruned_dmat->n * sizeof(double*));
    ASSERT(pruned_dmat->matrix != NULL);
    // rest of matrix = // TBD
    // pruned_dmat->drop =  // NO NEED, NULL above

#if DEV==1
    bool has_drop = dmat->drop != NULL;
    ASSERT(has_drop);
    has_drop = false;
    for (size_t matrix_i = 0;matrix_i < dmat->n;++matrix_i) {
        if (dmat->drop[matrix_i] != NULL) {
            has_drop = true;
            break;
        }
    }
    ASSERT(has_drop);
#endif


    pruned_dmat->names = strArray_init();

    totnPairsToPrune = 0;
    idx = 0;
    for (size_t i = 1;i < (size_t)nItems;++i) {
        for (size_t j = 0;j < i;++j) {
            if (dmat->drop[ni][idx]) {
                nDrops[i]++;
                nDrops[j]++;
                totnPairsToPrune++;
            }
            ++idx;
        }
    }

    int new_nItems;
    int nIndsRemoved = 0;

    if (0 == totnPairsToPrune) {
        new_nItems = nItems;

    } else if (totnPairsToPrune > 0) {

        LOG("Program will remove %d out of %ld pairwise distances.", totnPairsToPrune, dmat->size);

        while (1) {

            maxDrop = 0;
            maxDropIdx = -1;
            for (size_t i = 0;i < (size_t)nItems;++i) {
                if (nDrops[i] > maxDrop) {
                    maxDrop = nDrops[i];
                    maxDropIdx = i;
                }
            }

            if (maxDrop == 0) {
                ASSERT(totnPairsToPrune == 0);
                break;
            }
            ASSERT(totnPairsToPrune > 0);

            // remove individual with max drop count
            rmIdx[nIndsRemoved] = maxDropIdx;
            ASSERT(maxDropIdx > -1);
            nDrops[maxDropIdx] = 0;

            LOG("Removing individual %s, who is in %d out of %d remaining dropped pairs ", dmat->names->d[rmIdx[nIndsRemoved]], maxDrop, totnPairsToPrune);

            idx = 0;
            for (size_t i = 1;i < (size_t)nItems;++i) {
                for (size_t j = 0;j < i;++j) {
                    if (i == (size_t)maxDropIdx) {
                        if (nDrops[j] > 0) {
                            nDrops[j]--;
                            --totnPairsToPrune;
                        }
                    } else if (j == (size_t)maxDropIdx) {
                        if (nDrops[i] > 0) {
                            nDrops[i]--;
                            --totnPairsToPrune;
                        }
                    }
                    ++idx;
                }
            }
            ++nIndsRemoved;
        }


#if DEV==1
        ASSERT(totnPairsToPrune == 0);
        for (size_t i = 0;i < (size_t)nItems;++i) {
            if (nDrops[i] > 0) {
                NEVER;
            }
            bool isRemoved;
            for (size_t k = 0;k < (size_t)nIndsRemoved;++k) {
                if ((int)i == rmIdx[k]) {
                    isRemoved = true;
                    break;
                }
            }
            if (isRemoved) {
                continue;
            }
            if (dmat->drop[0][i]) {
                NEVER;
            }
        }
#endif

        // copy the remaining data to the new pruned matrix
        new_nItems = nItems - nIndsRemoved;
        ASSERT(new_nItems > 0);
        LOG("Removed %d individuals from the distance matrix. Pruned matrix will have %d individuals", nIndsRemoved, new_nItems);

        kstring_t tmp = KS_INIT;
        ksprintf(&tmp, "List of removed individuals:\n");
        for (size_t i = 0;i < (size_t)nIndsRemoved;++i) {
            ksprintf(&tmp, "%s\n", dmat->names->d[rmIdx[i]]);
        }
        LOG("%s", tmp.s);
        ks_free(&tmp);

    } else {
        NEVER;
    }


    pruned_dmat->size = (new_nItems * (new_nItems - 1)) / 2;
    pruned_dmat->matrix[ni] = (double*)malloc(pruned_dmat->size * sizeof(double));
    ASSERT(pruned_dmat->matrix[0] != NULL);
    for (size_t i = 0;i < pruned_dmat->size;++i) {
        pruned_dmat->matrix[0][i] = 0.0;
    }

    bool isRemoved;

    size_t newidx = 0;
    for (size_t i = 0;i < (size_t)nItems;++i) {

        isRemoved = false;
        for (size_t k = 0;k < (size_t)nIndsRemoved;++k) {
            if ((int)i == rmIdx[k]) {
                isRemoved = true;
                break;
            }
        }
        if (isRemoved) {
            continue;
        }

        pruned_dmat->names->add(dmat->names->d[i]);

        for (size_t j = 0;j < i;++j) {

            isRemoved = false;
            for (size_t k = 0;k < (size_t)nIndsRemoved;++k) {
                if ((int)j == rmIdx[k]) {
                    isRemoved = true;
                    break;
                }
            }
            if (isRemoved) {
                continue;
            }

            ASSERT(dmat->type == DMAT_TYPE_LTED); // TODO add support for other types
            pruned_dmat->matrix[ni][newidx] = dmat->matrix[ni][MATRIX_GET_INDEX_LTED_IJ(i, j)];
            ++newidx;
        }
    }

    DEVASSERT(newidx == pruned_dmat->size);

    if (0 == nIndsRemoved) {
        WARN("(--prune-dmat) No individuals were removed from the distance matrix. The pruned matrix is identical to the original matrix.");

    }

    END_LOGSECTION;

    FREE(nDrops);
    FREE(rmIdx);

    return(pruned_dmat);
}

dmat_t* dmat_init(const size_t n, const size_t nInd, const uint8_t type, const uint32_t method, uint8_t transform, strArray* names, const uint8_t names_src) {

    dmat_t* dmat = NULL;
    dmat = (dmat_t*)malloc(sizeof(dmat_t));
    ASSERT(dmat != NULL);

    // -> init
    dmat->n = 0;
    dmat->size = 0;
    dmat->type = 0;
    dmat->transform = DMAT_INTPLUS_TRANSFORM_NONE;
    dmat->method = 0;
    dmat->names = NULL;
    dmat->names_src = 0;
    dmat->drop = NULL;


    // -> set 
    dmat->n = n;

    switch (type) {
    case DMAT_TYPE_LTED:
    case DMAT_TYPE_UTED:
        dmat->size = (nInd * (nInd - 1)) / 2;
        break;
    case DMAT_TYPE_LTID:
    case DMAT_TYPE_UTID:
        dmat->size = (nInd * (nInd + 1)) / 2;
        break;
    case DMAT_TYPE_FULL:
        dmat->size = nInd * nInd;
        break;
    default:
        NEVER;
    }

    dmat->type = type;
    dmat->transform = transform;
    dmat->method = method;


    dmat->names_src = names_src;
    if (dmat->names_src == DMAT_NAMES_SRC_IN_VCF_PARS_PTR) {
        dmat->names = names;
    } else if (dmat->names_src == DMAT_NAMES_SRC_IN_METADATA_NAMES_PTR) {
        dmat->names = names;
    } else if (dmat->names_src == DMAT_NAMES_SRC_PRIVATE) {
        // names is allocated and used internally in the program
    } else {
        NEVER;
    }

    dmat->matrix = (double**)malloc(dmat->n * sizeof(double*));
    ASSERT(dmat->matrix != NULL);

    //TODO allocate only if needed? and check has_drop?
    dmat->drop = NULL;
    //dmat->drop = (bool**)malloc(dmat->n * sizeof(bool*));
    //ASSERT(dmat->drop != NULL);
    //for (size_t i = 0; i < dmat->n; ++i) {
    //    dmat->drop[i] = NULL;
    //    dmat->drop[i] = (bool*)malloc(dmat->size * sizeof(bool));
    //    ASSERT(dmat->drop[i] != NULL);
    //    for (size_t j = 0; j < dmat->size; ++j) {
    //        dmat->drop[i][j] = false;
    //    }
    //}

    for (size_t i = 0; i < dmat->n;++i) {
        dmat->matrix[i] = NULL;
        dmat->matrix[i] = (double*)malloc(dmat->size * sizeof(double));
        ASSERT(dmat->matrix[i] != NULL);
        for (size_t j = 0;j < dmat->size;++j) {
            dmat->matrix[i][j] = 0.0;
        }
    }
    return(dmat);
}


void dmat_destroy(dmat_t* dmat) {
    for (size_t i = 0; i < dmat->n;++i) {
        FREE(dmat->matrix[i]);
    }
    FREE(dmat->matrix);

    if (dmat->drop != NULL) {
        // isNULL if dmat is pruned dmat
        for (size_t i = 0; i < dmat->n;++i) {
            FREE(dmat->drop[i]);
        }
        FREE(dmat->drop);
    }

    if (dmat->names_src == DMAT_NAMES_SRC_IN_DM_FILE) {
        strArray_destroy(dmat->names);
    } else if (dmat->names_src == DMAT_NAMES_SRC_IN_VCF_PARS_PTR) {
        dmat->names = NULL;
    } else if (dmat->names_src == DMAT_NAMES_SRC_IN_METADATA_NAMES_PTR) {
        dmat->names = NULL;
    } else if (dmat->names_src == DMAT_NAMES_SRC_NONE) {
        NEVER;
    } else if (dmat->names_src == DMAT_NAMES_SRC_PRIVATE) {
        strArray_destroy(dmat->names);
    } else {
        NEVER;
    }

    FREE(dmat);
    return;
}


inline void dmat_matrix_apply_transform(double** matrix, size_t n, size_t size, uint8_t type, bool** drop, uint8_t matrix_transform, uint8_t required_transform) {
    DEVASSERT(matrix != NULL);

    if (required_transform == matrix_transform) {
        return;
    }

    if (required_transform != DMAT_INTPLUS_TRANSFORM_NONE && required_transform != DMAT_INTPLUS_TRANSFORM_SQUARE) {
        // requested a complex transformation (not simply squaring or no transform)

        if (drop != NULL) {
            NEVER;
        }
        if (required_transform & DMAT_INTPLUS_TRANSFORM_SQUARE) {
            LOG("Both Cailliez and squaring transformations requested. Program will first apply Cailliez transformation and then square the distances.");
        }

        static const double tole = 2.220446e-16;
        static const bool cor_zero = true;
        size_t nelems;
        if (type == DMAT_TYPE_FULL) {
            nelems = sqrt(size);
        } else {
            NEVER;
        }
        for (size_t matrix_i = 0;matrix_i < n;++matrix_i) {
            double* cmatrix = (double*)malloc(size * sizeof(double));
            ASSERT(cmatrix != NULL);
            cailliez(matrix[matrix_i], nelems, cmatrix, tole, cor_zero);
            for (size_t i = 0;i < size;++i) {
                matrix[matrix_i][i] = cmatrix[i];
                // if requested, do the squaring last
                if (required_transform & DMAT_INTPLUS_TRANSFORM_SQUARE) {
                    matrix[matrix_i][i] *= matrix[matrix_i][i];
                }
            }
            FREE(cmatrix);
        }

    } else {
        if (required_transform == DMAT_INTPLUS_TRANSFORM_NONE) {
            // -> we need: no transform
            if (matrix_transform == DMAT_INTPLUS_TRANSFORM_SQUARE) {
                // -> we have: squared
                for (size_t m = 0;m < n;++m) {
                    for (size_t i = 0;i < size;++i) {
#if DEV==1
                        if (matrix[m][i] < 0.0) {
                            ERROR("Negative distance value detected after square root transformation.");
                        }
                        if (drop != NULL && drop[m][i] && !isnan(matrix[m][i])) {
                            ERROR("Distance value is not NaN for a pair with missing values.");
                        }
#endif
                        if (drop == NULL || drop[m] == NULL || !drop[m][i]) {
                            matrix[m][i] = sqrt(matrix[m][i]);
                        }
                    }
                }
            } else {
                // cannot undo any other transformation
                NEVER;
            }
        } else if (required_transform == DMAT_INTPLUS_TRANSFORM_SQUARE) {
            // -> we need: squared
            if (matrix_transform == DMAT_INTPLUS_TRANSFORM_NONE) {
                // -> we have: no transform
                for (size_t m = 0;m < n;++m) {
                    for (size_t i = 0;i < size;++i) {
#if DEV==1
                        if (matrix[m][i] < 0.0) {
                            ERROR("Negative distance value detected after square root transformation.");
                        }
                        if (drop != NULL && drop[m][i] && !isnan(matrix[m][i])) {
                            ERROR("Distance value is not NaN for a pair with missing values.");
                        }
#endif
                        if (drop == NULL || drop[m] == NULL || !drop[m][i]) {
                            matrix[m][i] *= matrix[m][i];
                        }
                    }
                }
            } else {
                NEVER;
            }
        } else {
            NEVER;
        }
    }

    return;
}

double** new_dmat_matrix_apply_transform(dmat_t* dmat, uint8_t required_transform) {

    ASSERT(dmat != NULL);
    if (dmat->transform == required_transform) {
        NEVER;
    }

    // transformation is needed, create a new matrix with applied transformation
    double** matrix = NULL;

    if (required_transform != DMAT_INTPLUS_TRANSFORM_NONE && required_transform != DMAT_INTPLUS_TRANSFORM_SQUARE) {
        // if requested a complex transformation (not simply squaring or no transform)
        if (dmat->drop != NULL) {
            ERROR("Complex transformation requested on a distance matrix that contains individual pairs with missing values. Please use --prune-dm 1 to prune the distance matrix before applying the transformation.");
        }
        if (dmat->type != DMAT_TYPE_FULL) {
            ERROR("Only full distance matrices are supported for complex transformations.");
        }
    }
    // create a new matrix of same size
    matrix = (double**)malloc(dmat->n * sizeof(double*));
    ASSERT(matrix != NULL);
    for (size_t i = 0;i < dmat->n;++i) {
        matrix[i] = NULL;
        matrix[i] = (double*)malloc(dmat->size * sizeof(double));
        ASSERT(matrix[i] != NULL);
        for (size_t j = 0;j < dmat->size;++j) {
            matrix[i][j] = dmat->matrix[i][j];
        }
    }

    dmat_matrix_apply_transform(matrix, dmat->n, dmat->size, dmat->type, dmat->drop, dmat->transform, required_transform);
    return(matrix);
}

void dmat_apply_transform(dmat_t* dmat, uint8_t required_transform) {
    if (dmat->transform == required_transform) {
        return;
    }
    // if a complex transformation is needed
    if (required_transform != DMAT_INTPLUS_TRANSFORM_NONE && required_transform != DMAT_INTPLUS_TRANSFORM_SQUARE) {
        if (dmat->type != DMAT_TYPE_FULL) {
            ERROR("Only full distance matrices are supported for complex transformations.");
        }
        double** matrix = new_dmat_matrix_apply_transform(dmat, required_transform);
        // a new matrix was created, free the old one and replace it with the new one
        for (size_t i = 0;i < dmat->n;++i) {
            FREE(dmat->matrix[i]);
        }
        FREE(dmat->matrix);
        dmat->matrix = matrix;
    } else {
        // if a simple transformation is needed
        // no need to create a new matrix, just apply the transformation to the existing matrix
        dmat_matrix_apply_transform(dmat->matrix, dmat->n, dmat->size, dmat->type, dmat->drop, dmat->transform, required_transform);
    }
    dmat->transform = required_transform;
    return;
}

dmat_t* new_dmat_apply_transform(dmat_t* dmat, uint8_t required_transform) {
    DEVASSERT(dmat->transform != required_transform);
    // if a complex transformation is needed
    if (required_transform != DMAT_INTPLUS_TRANSFORM_NONE && required_transform != DMAT_INTPLUS_TRANSFORM_SQUARE) {
        if (dmat->type != DMAT_TYPE_FULL) {
            ERROR("Only full distance matrices are supported for complex transformations.");
        }
    }
    double** matrix = new_dmat_matrix_apply_transform(dmat, required_transform);

    size_t nelems;
    if (dmat->type == DMAT_TYPE_FULL) {
        nelems = sqrt(dmat->size);
    } else if (dmat->type == DMAT_TYPE_LTED) {
        nelems = (1 + sqrt(1 + 8 * dmat->size)) / 2;
    } else {
        NEVER;
    }
    dmat_t* new_dmat = dmat_init(dmat->n, nelems, dmat->type, dmat->method, required_transform, dmat->names, dmat->names_src);
    const bool fill_drop = new_dmat->drop != NULL;
    for (size_t matrix_i = 0;matrix_i < new_dmat->n;++matrix_i) {
        for (size_t i = 0;i < new_dmat->size;++i) {
            new_dmat->matrix[matrix_i][i] = matrix[matrix_i][i];
            if (fill_drop) {
                new_dmat->drop[matrix_i][i] = dmat->drop[matrix_i][i];
            }
        }
        FREE(matrix[matrix_i]);
    }
    FREE(matrix);
    return(new_dmat);
}

dmat_t* dmat_read(const char* in_dm_fn, uint8_t required_transform, metadata_t* metadata) {

    LOG("(--in-dm) Reading distance matrix from file: %s", in_dm_fn);

    FILE* fp = NULL;

    bool isGzFile = false;
    if (IO::isGzFile(args->in_dm_fn) == 1) {
        isGzFile = true;
        kstring_t cmd = KS_INIT;
        ksprintf(&cmd, "zcat %s", args->in_dm_fn);
        fp = popen(cmd.s, "r");
    } else {
        fp = IO::getFile(args->in_dm_fn, "r");
    }

    dmat_t* dmat = NULL;
    dmat = (dmat_t*)malloc(sizeof(dmat_t));
    ASSERT(dmat != NULL);

    // -> init
    dmat->n = 0;
    dmat->size = 0;
    dmat->type = 0;
    dmat->transform = DMAT_INTPLUS_TRANSFORM_NONE;
    dmat->method = 0;
    dmat->names = NULL;
    dmat->names_src = 0;
    dmat->drop = NULL;

    // -> set
    if (metadata != NULL) {
        dmat->names_src = DMAT_NAMES_SRC_IN_METADATA_NAMES_PTR;
    } else {
        dmat->names_src = DMAT_NAMES_SRC_IN_DM_FILE;
    }

    int lineno = 0;

    BEGIN_LOGSECTION;

    // line 0: number of matrices
    if (fscanf(fp, "%lu", &dmat->n) != 1) {
        // ERROR()
        ERROR("Could not read the number of matrices line (line 1) from the distance matrix input file %s", args->in_dm_fn);
    }
    if (0 == dmat->n) {
        ERROR("Distance matrix input file must contain at least 1 distance matrices (line 1 > 0)");
    }

    LOG("Reading %ld distance matrices from the distance matrix file", dmat->n);

    ++lineno;

    // line 1: type
    char in_type[5] = { '\0' };

    if (fscanf(fp, "%s", in_type) != 1) {
        ERROR("Could not read the type line (line 2) from the distance matrix input file %s", args->in_dm_fn);
    }
    if (in_type[4] != '\0') {
        ERROR("Unrecognized distance matrix type: %s", in_type);
    }

    if (in_type[0] == 'F' && in_type[1] == 'U' && in_type[2] == 'L' && in_type[2] == 'L') {
        dmat->type = DMAT_TYPE_FULL;
        LOG("Input distance matrix type is detected as FULL: Full Matrix");
    } else if (in_type[0] == 'L' && in_type[1] == 'T' && in_type[2] == 'E' && in_type[3] == 'D') {
        dmat->type = DMAT_TYPE_LTED;
        LOG("Input distance matrix type is detected as LTED: Lower Triangular Matrix (Excluding Diagonal)");
    } else if (in_type[0] == 'L' && in_type[1] == 'T' && in_type[2] == 'I' && in_type[3] == 'D') {
        dmat->type = DMAT_TYPE_LTID;
        LOG("Input distance matrix type is detected as LTID: Lower Triangular Matrix (Including Diagonal)");
    } else if (in_type[0] == 'U' && in_type[1] == 'T' && in_type[2] == 'E' && in_type[3] == 'D') {
        dmat->type = DMAT_TYPE_UTED;
        LOG("Input distance matrix type is detected as UTED: Upper Triangular Matrix (Excluding Diagonal)");
    } else if (in_type[0] == 'U' && in_type[1] == 'T' && in_type[2] == 'I' && in_type[3] == 'D') {
        dmat->type = DMAT_TYPE_UTID;
        LOG("Input distance matrix type is detected as UTID: Upper Triangular Matrix (Including Diagonal)");
    } else {
        ERROR("Unrecognized distance matrix type: %s", in_type);
    }

    ++lineno;

    // line 2: transformation
    int tmp_transform;
    if (fscanf(fp, "%d", &tmp_transform) != 1) {
        ERROR("Could not read the transform line (line 3) from the distance matrix input file %s", args->in_dm_fn);
    }
    if (tmp_transform < 0 || tmp_transform >= 255) {
        ERROR("Unrecognized distance matrix transformation: %d", tmp_transform);
    }
    dmat->transform = (uint8_t)tmp_transform;

    if (DMAT_INTPLUS_TRANSFORM_NONE == dmat->transform) {
        LOG("Input distance matrix transform is detected as: None");
    } else {
        kstring_t tmp = KS_INIT;
        ksprintf(&tmp, "Input distance matrix transform is detected as: ");
        if (dmat->transform & DMAT_INTPLUS_TRANSFORM_SQUARE) {
            ksprintf(&tmp, "Squared");
        }
        if (dmat->transform & DMAT_INTPLUS_TRANSFORM_QUASIEUCLID) {
            ksprintf(&tmp, ", Quasi-Euclid Transformed");
        }
        if (dmat->transform & DMAT_INTPLUS_TRANSFORM_LINGOES) {
            ksprintf(&tmp, ", Lingoes Transformed");
        }
        if (dmat->transform & DMAT_INTPLUS_TRANSFORM_CAILLIEZ) {
            ksprintf(&tmp, ", Cailliez Transformed");
        }

        LOG("Input distance matrix transform is detected as: %s", tmp.s);
        ks_free(&tmp);
    }


    ++lineno;

    // line 3: method
    if (fscanf(fp, "%u", &dmat->method) != 1) {
        ERROR("Could not read the method line (line 4) from the distance matrix input file %s", args->in_dm_fn);
    }


    LOG("Input distance matrix method is detected as: %s (%d)", get_dmat_method_str(dmat->method), dmat->method);

    ++lineno;

    // line 4: number of names
    int in_nInd;
    if (fscanf(fp, "%d", &in_nInd) != 1) {
        ERROR("Could not read the second line from the distance matrix input file %s", args->in_dm_fn);
    }

    if (in_nInd > 0) {
        LOG("Found %d names in the distance matrix file", in_nInd);
    } else {
        ERROR("Number of names in the distance matrix file detected as %d. It must be greater than 0.", in_nInd);
    }

    ++lineno;

    // line 5: number of distances 
    int in_nIndCmb;
    if (fscanf(fp, "%d", &in_nIndCmb) != 1) {
        ERROR("Could not read the third line from the distance matrix input file %s", args->in_dm_fn);
    }


    if (in_nIndCmb > 0) {
        LOG("Found %d distances in the distance matrix file", in_nIndCmb);
    } else {
        ERROR("Number of distances in the distance matrix file detected as %d. It must be greater than 0.", in_nIndCmb);
    }

    dmat->size = in_nIndCmb;

    if (dmat->type == DMAT_TYPE_LTED) {
        if ((int)dmat->size != ((in_nInd * (in_nInd - 1)) / 2)) {
            ERROR("Number of distances in the input distance matrix (%ld) does not match the expected number of pairwise distances given the number of individuals (%d) and the distance matrix type (LTED).", dmat->size, in_nInd);
        }
    } else if (dmat->type == DMAT_TYPE_UTED) {
        if ((int)dmat->size != ((in_nInd * (in_nInd - 1)) / 2)) {
            ERROR("Number of distances in the input distance matrix (%ld) does not match the expected number of pairwise distances given the number of individuals (%d) and the distance matrix type (UTED).", dmat->size, in_nInd);
        }
    } else if (dmat->type == DMAT_TYPE_LTID) {
        if ((int)dmat->size != ((in_nInd * (in_nInd + 1)) / 2)) {
            ERROR("Number of distances in the input distance matrix (%ld) does not match the expected number of pairwise distances given the number of individuals (%d) and the distance matrix type (LTID).", dmat->size, in_nInd);
        }
    } else if (dmat->type == DMAT_TYPE_UTID) {
        if ((int)dmat->size != ((in_nInd * (in_nInd + 1)) / 2)) {
            ERROR("Number of distances in the input distance matrix (%ld) does not match the expected number of pairwise distances given the number of individuals (%d) and the distance matrix type (UTID).", dmat->size, in_nInd);
        }
    } else if (dmat->type == DMAT_TYPE_FULL) {
        if ((int)dmat->size != (in_nInd * in_nInd)) {
            ERROR("Number of distances in the input distance matrix (%ld) does not match the expected number of pairwise distances given the number of individuals (%d) and the distance matrix type (FULL).", dmat->size, in_nInd);
        }
    } else {
        NEVER;
    }

    ++lineno;


    /// ------------------------------------------------------------
    // lines 6:(6+nInd-1) names

    if (dmat->names_src == DMAT_NAMES_SRC_IN_DM_FILE) {
        dmat->names = strArray_alloc(in_nInd);

        char tmp[1024] = { '\0' };
        while (lineno < 6 + in_nInd) {
            tmp[0] = '\0';
            fscanf(fp, "%s\n", tmp);
            if ('\0' == tmp[0]) {
                ERROR("Found bad name at line %d of distance matrix file %s", lineno + 1, args->in_dm_fn);
            }
            dmat->names->add(tmp);
            ++lineno;
        }

        ASSERT((int)dmat->names->len == in_nInd);
        LOG("Read %ld names from distance matrix file", dmat->names->len);

        ++lineno;

    } else if (dmat->names_src == DMAT_NAMES_SRC_IN_METADATA_NAMES_PTR) {

        if (in_nInd != (int)metadata->indNames->len) {
            ERROR("Number of names in the distance matrix file (%d) does not match the number of individuals in the metadata file (%ld). Please edit your metadata file to match the names in the distance matrix file.", in_nInd, metadata->indNames->len);
        }

        char tmp[1024] = { '\0' };
        bool badorder = false;
        size_t mtdidx;
        while (lineno < 6 + in_nInd) {
            tmp[0] = '\0';
            fscanf(fp, "%s\n", tmp);
            if ('\0' == tmp[0]) {
                ERROR("Found bad name at line %d of distance matrix file %s", lineno + 1, args->in_dm_fn);
            }
            if (metadata->indNames->find(tmp, &mtdidx)) {
                if ((int)mtdidx != lineno - 6) {
                    DEVPRINT("name:%s idx:%d mtdidx:%ld mtdname:%s", tmp, lineno - 6, mtdidx, metadata->indNames->d[mtdidx]);
                    badorder = true;
                    break;
                }
            } else {
                ERROR("Name %s in the distance matrix file is not found in the metadata file.", tmp);
            }
            ++lineno;
        }
        if (badorder) {
            ERROR("Names in the distance matrix file are not in the same order as in the metadata file. Please edit your metadata file to match the names in the distance matrix file.");
        }

        dmat->names = metadata->indNames;

    } else {
        NEVER;
    }

    /// ------------------------------------------------------------
    // lines [6-X): dmat->names
    // X = 6 + (line 4)
    // lines (6+nInd):(6+nInd+nIndCmb-1) distances
    dmat->matrix = (double**)malloc(dmat->n * sizeof(double*));
    ASSERT(dmat->matrix != NULL);

    for (size_t i = 0; i < dmat->n;++i) {
        dmat->matrix[i] = NULL;
        dmat->matrix[i] = (double*)malloc(dmat->size * sizeof(double));
        ASSERT(dmat->matrix[i] != NULL);
        for (size_t j = 0;j < dmat->size;++j) {
            dmat->matrix[i][j] = 0.0;
        }
    }

    // lines [X, Y): distances
    if (dmat->n > 1) {
        //TODO multidm reading still reads into matrix[0] ???
        NEVER;
    }
    size_t k = 0;
    while ((int)lineno < (int)(7 + ((in_nInd + in_nIndCmb) * dmat->n))) {
        fscanf(fp, "%lf\n", &dmat->matrix[0][k]);
        ++lineno;
        ++k;
    }


    ASSERT((int)dmat->size == in_nIndCmb);
    LOG("Read %ld distances from distance matrix file", dmat->size);

    // -> set drop based on nan values

    // bool indicating whether we need dmat->drop
    // we need dmat->drop if any of the pairs are dropped  for matrix i
    // if we don't need to drop any pairs, we don't need dmat->drop[matrix_i]
    bool need_drop;
    bool any_need_drop = false;
    for (size_t i = 0; i < dmat->n; ++i) {
        need_drop = false;
        for (size_t j = 0; j < dmat->size; ++j) {
            if (isnan(dmat->matrix[i][j])) {
                if (need_drop == false) {
                    if (any_need_drop == false) {
                        // not even allocated mem for dmat->drop yet
                        dmat->drop = NULL;
                        dmat->drop = (bool**)malloc(dmat->n * sizeof(bool*));
                        ASSERT(dmat->drop != NULL);
                        for (size_t k = 0;k < dmat->n;++k) {
                            dmat->drop[k] = NULL;
                        }
                    }
                    // not allocated mem for dmat->drop[matrix_i] yet
                    dmat->drop[i] = NULL;
                    dmat->drop[i] = (bool*)malloc(dmat->size * sizeof(bool));
                    ASSERT(dmat->drop[i] != NULL);

                    // fill everything until this point with false
                    for (size_t k = 0;k < j;++k) {
                        dmat->drop[i][k] = false;
                    }
                }
                dmat->drop[i][j] = true;
                need_drop = true;
                any_need_drop = true;
            } else {
                if (need_drop) {
                    // allocated mem, fill with false
                    dmat->drop[i][j] = false;

                    // else not allocated mem yet, no need to drop;
                    // bc maybe this matrix has no missing values
                }
            }
        }
    }

    dmat_apply_transform(dmat, required_transform);

    if (isGzFile) {
        PCLOSE(fp);
    } else {
        FCLOSE(fp);
    }

    // -> finished reading the distance matrix file

    END_LOGSECTION;

    return(dmat);
}


void dmat_print(dmat_t* dmat, outfile_t* outfile) {

    LOG("(--print-dm) Printing distance matrix to file: %s", outfile->fn);

    kstring_t* kbuf = &outfile->kbuf;

    //TODO instead of printing numbers directly without explanation, make it a bit more user friendly, change number to text when possible
    // line 0: number of matrices
    ksprintf(kbuf, "%ld\n", dmat->n);

    // line 1: type
    if (DMAT_TYPE_LTED == dmat->type) {
        ksprintf(kbuf, "%s\n", "LTED");
    } else if (DMAT_TYPE_LTID == dmat->type) {
        ksprintf(kbuf, "%s\n", "LTID");
    } else if (DMAT_TYPE_UTED == dmat->type) {
        ksprintf(kbuf, "%s\n", "UTED");
    } else if (DMAT_TYPE_UTID == dmat->type) {
        ksprintf(kbuf, "%s\n", "UTID");
    } else if (DMAT_TYPE_FULL == dmat->type) {
        ksprintf(kbuf, "%s\n", "FULL");
    } else {
        NEVER;
    }

    // line 2: transformation
    ksprintf(kbuf, "%d\n", dmat->transform);

    // line 3: method
    ksprintf(kbuf, "%d\n", dmat->method);

    // line 4: number of dmat->names
    ksprintf(kbuf, "%ld\n", dmat->names->len);

    // line 5: number of distances 
    ksprintf(kbuf, "%ld\n", dmat->size);

    // lines [6-X): dmat->names
    // X = 6 + (line 4)
    for (size_t i = 0;i < dmat->names->len;++i) {
        ksprintf(kbuf, "%s\n", dmat->names->d[i]);
    }

    double* matrix = NULL;

    for (size_t mi = 0;mi < dmat->n;++mi) {
        matrix = dmat->matrix[mi];

        // lines [X, Y): distances
        for (size_t p = 0; p < dmat->size; ++p) {
            ksprintf(kbuf, "%.17g\n", matrix[p]);
        }
    }

    return;
}

void dmat_print_verbose(dmat_t* dmat, outfile_t* outfile) {

    kstring_t* kbuf = &outfile->kbuf;

    ksprintf(kbuf, "## This file was produced by: %s\n", PROGRAM_VERSION_INFO);
    ksprintf(kbuf, "## Command: %s\n", PROGRAM_COMMAND);
    // number of matrices
    ksprintf(kbuf, "## Number of distance matrices: %ld\n", dmat->n);

    ksprintf(kbuf, "## -> Following information is shared for all distance matrices in this file:\n");

    // type
    ksprintf(kbuf, "# Type of the distance matrix: ");
    if (DMAT_TYPE_LTED == dmat->type) {
        ksprintf(kbuf, "%s\n", "Lower Triangular Matrix (Excluding Diagonal)");
    } else if (DMAT_TYPE_LTID == dmat->type) {
        ksprintf(kbuf, "%s\n", "Lower Triangular Matrix (Including Diagonal)");
    } else if (DMAT_TYPE_UTED == dmat->type) {
        ksprintf(kbuf, "%s\n", "Upper Triangular Matrix (Excluding Diagonal)");
    } else if (DMAT_TYPE_UTID == dmat->type) {
        ksprintf(kbuf, "%s\n", "Upper Triangular Matrix (Including Diagonal)");
    } else if (DMAT_TYPE_FULL == dmat->type) {
        ksprintf(kbuf, "%s\n", "Full Matrix");
    } else {
        NEVER;
    }

    // transformation
    ksprintf(kbuf, "# Transformation applied to the distances: ");
    if (DMAT_INTPLUS_TRANSFORM_NONE == dmat->transform) {
        ksprintf(kbuf, "%s", "None");
    } else {
        kstring_t tmp = KS_INIT;
        if (dmat->transform & DMAT_INTPLUS_TRANSFORM_SQUARE) {
            ksprintf(&tmp, "Squared");
        }
        if (dmat->transform & DMAT_INTPLUS_TRANSFORM_QUASIEUCLID) {
            ksprintf(&tmp, ", Quasi-Euclid Transformed");
        }
        if (dmat->transform & DMAT_INTPLUS_TRANSFORM_LINGOES) {
            ksprintf(&tmp, ", Lingoes Transformed");
        }
        if (dmat->transform & DMAT_INTPLUS_TRANSFORM_CAILLIEZ) {
            ksprintf(&tmp, ", Cailliez Transformed");
        }
        ksprintf(kbuf, "%s", tmp.s);
        ks_free(&tmp);
    }

    ksprintf(kbuf, " (--dm-transform %d)", dmat->transform);
    kputc('\n', kbuf);


    // method
    ksprintf(kbuf, "# Method used to calculate the distances: ");
    ksprintf(kbuf, "%s", get_dmat_method_str(dmat->method));
    ksprintf(kbuf, " (--dm-method %d)", DMAT_METHOD_DIJ);
    kputc('\n', kbuf);

    const size_t dmat_n_names = dmat->names->len;
    // number of dmat->names
    ksprintf(kbuf, "# Number of items in the distance matrix: %ld\n", dmat_n_names);

    // number of distances 
    ksprintf(kbuf, "# Number of distances in the distance matrix: %ld\n", dmat->size);

    // dmat->names
    ksprintf(kbuf, "# NAME, Names of the items in the distance matrix:\n");
    ksprintf(kbuf, "%-6s\t%-8s\t%-25s\n", "# NAME", "[2]INDEX", "[3]NAME");
    for (size_t i = 0;i < dmat_n_names;++i) {
        ksprintf(kbuf, "%-6s\t%-8ld\t%-25s\n", "NAME", i, dmat->names->d[i]);
    }

    double* matrix = NULL;

    for (size_t mi = 0;mi < dmat->n;++mi) {

        // we don't expect a distance matrix file to have more than 99999999 matrices; and 8+2 for parantheses+1 for \0 = 11
        char mi_str[11] = { 0 };
        snprintf(mi_str, sizeof(mi_str), "(%ld)", mi);


        ksprintf(kbuf, "## -> Following information is specific to distance matrix with MATRIX_INDEX=%ld in this file:\n", mi);

        ksprintf(kbuf, "# DIST(MATRIX_INDEX), Pairwise distances in the distance matrix:\n");
        ksprintf(kbuf, "# DIST%-11.10s\t%-13s\t%-19s\t%-19s\t%-25s\t%-25s\t%-25s\n", mi_str, "[4]PAIR_INDEX", "[5]FIRST_ITEM_INDEX", "[6]SECOND_ITEM_INDEX", "[7]FIRST_ITEM_NAME", "[8]SECOND_ITEM_NAME", "[9]DISTANCE");

        matrix = dmat->matrix[mi];

        if (dmat->type == DMAT_TYPE_LTED) {

            size_t idx = 0;
            for (size_t i = 0;i < dmat_n_names;++i) {
                for (size_t j = 0;j < i;++j) {
                    ksprintf(kbuf, "DIST%-13.10s\t%-13ld\t%-19ld\t%-19ld\t%-25s\t%-25s\t%-25.17g\n", mi_str, idx, i, j, dmat->names->d[i], dmat->names->d[j], matrix[idx]);
                    ++idx;
                }
            }

            ksprintf(kbuf, "# DMAT(MATRIX_INDEX), Distance matrix in matrix format:\n");
            ksprintf(kbuf, "# DMAT%-11.10s\t%-25s\t", mi_str, "[2]ITEMS");
            for (size_t i = 0;i < dmat_n_names;++i) {
                ksprintf(kbuf, "%-25s", dmat->names->d[i]);
                if (i != dmat_n_names - 1) {
                    kputc('\t', kbuf);
                }
            }
            kputc('\n', kbuf);
            idx = 0;
            for (size_t i = 0;i < dmat_n_names;++i) {
                ksprintf(kbuf, "DMAT%-13.10s\t", mi_str);
                ksprintf(kbuf, "%-25s\t", dmat->names->d[i]);
                for (size_t j = 0;j < i;++j) {
                    ksprintf(kbuf, "%-25.17g", matrix[idx]);
                    if (j != i - 1) {
                        kputc('\t', kbuf);
                    }
                    ++idx;
                }
                kputc('\n', kbuf);
            }
        } else if (dmat->type == DMAT_TYPE_FULL) {

            size_t idx = 0;
            for (size_t i = 0;i < dmat_n_names;++i) {
                for (size_t j = 0;j < dmat_n_names;++j) {
                    ksprintf(kbuf, "DIST%-13.10s\t%-13ld\t%-19ld\t%-19ld\t%-25s\t%-25s\t%-25.17g\n", mi_str, idx, i, j, dmat->names->d[i], dmat->names->d[j], matrix[idx]);
                    ++idx;
                }
            }

            ksprintf(kbuf, "# DMAT(MATRIX_INDEX), Distance matrix in matrix format:\n");
            ksprintf(kbuf, "# DMAT%-11.10s\t%-25s\t", mi_str, "[2]ITEMS");
            for (size_t i = 0;i < dmat_n_names;++i) {
                ksprintf(kbuf, "%-25s", dmat->names->d[i]);
                if (i != dmat_n_names - 1) {
                    kputc('\t', kbuf);
                }
            }
            kputc('\n', kbuf);
            idx = 0;
            for (size_t i = 0;i < dmat_n_names;++i) {
                ksprintf(kbuf, "DMAT%-13.10s\t", mi_str);
                ksprintf(kbuf, "%-25s\t", dmat->names->d[i]);
                for (size_t j = 0;j < dmat_n_names;++j) {
                    ksprintf(kbuf, "%-25.17g", matrix[idx]);
                    if (j != dmat_n_names - 1) {
                        kputc('\t', kbuf);
                    }
                    ++idx;
                }
                kputc('\n', kbuf);
            }
        } else {
            ERROR("Verbose printing is not supported for distance matrix type %d", dmat->type);
        }

    }

    return;
}


void dmat_calculate_distances(jgtmat_t* jgtmat, dmat_t* dmat, uint8_t required_transform) {

    double* dm = NULL;

    uint32_t method = dmat->method;
    const size_t nPairs = dmat->size;

    //TODO match pair indexes in jgtmat with dmat
    ASSERT(dmat->type == DMAT_TYPE_LTED); // TODO add matching with others, currently jgtmat is always LTED so dmat must be LTED too

    int nPairs_after_drop;

    // bool indicating whether we need dmat->drop
    // we need dmat->drop if any of the pairs are dropped  for matrix i
    // if we don't need to drop any pairs, we don't need dmat->drop[matrix_i]
    bool need_drop;
    bool any_need_drop = false;
    for (size_t i = 0;i < dmat->n;++i) {
        need_drop = false;

        nPairs_after_drop = nPairs;

        dm = dmat->matrix[i];


        if (PROGRAM_WILL_USE_BCF_FMT_GT) {

            double sum, x, A, B, C, D, E, F, G, H, I;
            uint64_t** jgtmat_mi = jgtmat->m[i];

            for (size_t p = 0; p < nPairs; ++p) {

                if (jgtmat->drop != NULL && jgtmat->drop[i][p]) {

                    if (need_drop == false) {
                        if (any_need_drop == false) {
                            // not even allocated mem for dmat->drop yet
                            dmat->drop = NULL;
                            dmat->drop = (bool**)malloc(dmat->n * sizeof(bool*));
                            ASSERT(dmat->drop != NULL);
                            for (size_t k = 0;k < dmat->n;++k) {
                                dmat->drop[k] = NULL;
                            }
                        }
                        // not allocated mem for dmat->drop[matrix_i] yet
                        dmat->drop[i] = NULL;
                        dmat->drop[i] = (bool*)malloc(dmat->size * sizeof(bool));
                        ASSERT(dmat->drop[i] != NULL);

                        // fill everything until this point with false
                        for (size_t k = 0;k < p;++k) {
                            dmat->drop[i][k] = false;
                        }
                    }
                    dmat->drop[i][p] = true;
                    need_drop = true;
                    any_need_drop = true;

                    dm[p] = NAN;
                    nPairs_after_drop--;

                    continue;
                } else {
                    if (need_drop) {
                        // allocated mem, fill with false
                        dmat->drop[i][p] = false;

                        // else not allocated mem yet, no need to drop;
                        // bc maybe this matrix has no missing values
                    }
                }

                uint64_t* jgtmat_mip = jgtmat_mi[p];
                A = jgtmat_mip[JGTMAT_LINEAR_IDXOF_A];
                B = jgtmat_mip[JGTMAT_LINEAR_IDXOF_B];
                C = jgtmat_mip[JGTMAT_LINEAR_IDXOF_C];
                D = jgtmat_mip[JGTMAT_LINEAR_IDXOF_D];
                E = jgtmat_mip[JGTMAT_LINEAR_IDXOF_E];
                F = jgtmat_mip[JGTMAT_LINEAR_IDXOF_F];
                G = jgtmat_mip[JGTMAT_LINEAR_IDXOF_G];
                H = jgtmat_mip[JGTMAT_LINEAR_IDXOF_H];
                I = jgtmat_mip[JGTMAT_LINEAR_IDXOF_I];
                sum = A + B + C + D + E + F + G + H + I;

                switch (method) {
                case DMAT_METHOD_DIJ:
                    x = JGTMAT_GET_DIJ;
                    break;
                case DMAT_METHOD_SIJ:
                    x = JGTMAT_GET_SIJ;
                    break;
                case DMAT_METHOD_FIJ:
                    x = JGTMAT_GET_FIJ;
                    break;
                case DMAT_METHOD_IBS0:
                    x = JGTMAT_GET_IBS0;
                    break;
                case DMAT_METHOD_IBS1:
                    x = JGTMAT_GET_IBS1;
                    break;
                case DMAT_METHOD_IBS2:
                    x = JGTMAT_GET_IBS2;
                    break;
                case DMAT_METHOD_R0:
                    x = JGTMAT_GET_R0;
                    break;
                case DMAT_METHOD_R1:
                    x = JGTMAT_GET_R1;
                    break;
                case DMAT_METHOD_KIN:
                    x = JGTMAT_GET_KIN;
                    break;
                default:
                    NEVER;
                }
                x /= sum;

                dm[p] = x;
            }

        } else if (PROGRAM_WILL_USE_BCF_FMT_GL) {

            double x, A, B, C, D, E, F, G, H, I;
            double** jgtmat_pmi = jgtmat->pm[i];

            for (size_t p = 0; p < nPairs; ++p) {
                double* jgtmat_pmip = jgtmat_pmi[p];

                if (jgtmat->drop != NULL && jgtmat->drop[i][p]) {


                    if (need_drop == false) {
                        if (any_need_drop == false) {
                            // not even allocated mem for dmat->drop yet
                            dmat->drop = NULL;
                            dmat->drop = (bool**)malloc(dmat->n * sizeof(bool*));
                            ASSERT(dmat->drop != NULL);
                            for (size_t k = 0;k < dmat->n;++k) {
                                dmat->drop[k] = NULL;
                            }
                        }
                        // not allocated mem for dmat->drop[matrix_i] yet
                        dmat->drop[i] = NULL;
                        dmat->drop[i] = (bool*)malloc(dmat->size * sizeof(bool));
                        ASSERT(dmat->drop[i] != NULL);

                        // fill everything until this point with false
                        for (size_t k = 0;k < p;++k) {
                            dmat->drop[i][k] = false;
                        }
                    }
                    dmat->drop[i][p] = true;
                    need_drop = true;
                    any_need_drop = true;

                    dm[p] = NAN;
                    nPairs_after_drop--;

                    continue;
                } else {
                    if (need_drop) {
                        // allocated mem, fill with false
                        dmat->drop[i][p] = false;

                        // else not allocated mem yet, no need to drop;
                        // bc maybe this matrix has no missing values
                    }
                }

                A = jgtmat_pmip[JGTMAT_LINEAR_IDXOF_A];
                B = jgtmat_pmip[JGTMAT_LINEAR_IDXOF_B];
                C = jgtmat_pmip[JGTMAT_LINEAR_IDXOF_C];
                D = jgtmat_pmip[JGTMAT_LINEAR_IDXOF_D];
                E = jgtmat_pmip[JGTMAT_LINEAR_IDXOF_E];
                F = jgtmat_pmip[JGTMAT_LINEAR_IDXOF_F];
                G = jgtmat_pmip[JGTMAT_LINEAR_IDXOF_G];
                H = jgtmat_pmip[JGTMAT_LINEAR_IDXOF_H];
                I = jgtmat_pmip[JGTMAT_LINEAR_IDXOF_I];


                switch (method) {
                case DMAT_METHOD_DIJ:
                    x = JGTMAT_GET_DIJ;
                    break;
                case DMAT_METHOD_SIJ:
                    x = JGTMAT_GET_SIJ;
                    break;
                case DMAT_METHOD_FIJ:
                    x = JGTMAT_GET_FIJ;
                    break;
                case DMAT_METHOD_IBS0:
                    x = JGTMAT_GET_IBS0;
                    break;
                case DMAT_METHOD_IBS1:
                    x = JGTMAT_GET_IBS1;
                    break;
                case DMAT_METHOD_IBS2:
                    x = JGTMAT_GET_IBS2;
                    break;
                case DMAT_METHOD_R0:
                    x = JGTMAT_GET_R0;
                    break;
                case DMAT_METHOD_R1:
                    x = JGTMAT_GET_R1;
                    break;
                case DMAT_METHOD_KIN:
                    x = JGTMAT_GET_KIN;
                    break;
                default:
                    NEVER;
                }

                dm[p] = x;
            }

        } else {
            NEVER;
        }

        if (nPairs_after_drop < args->min_n_pairs) {
            if (0 == i) {
                // original run
                ERROR("Number of pairs after dropping pairs according to --min-pairsites is %d, which is less than the minimum number of pairs required by --min-npairs (%d). Please adjust the --min-pairsites and --min-npairs parameters.", nPairs_after_drop, args->min_n_pairs);
            } else {
                // bootstrap run
                ERROR("(Bootstrap replicate %ld) Number of pairs after dropping pairs according to --min-pairsites is %d, which is less than the minimum number of pairs required by --min-npairs (%d). Please adjust the --min-pairsites and --min-npairs parameters.", i, nPairs_after_drop, args->min_n_pairs);
            }
        }

    }

    dmat_apply_transform(dmat, required_transform);

    return;
}
