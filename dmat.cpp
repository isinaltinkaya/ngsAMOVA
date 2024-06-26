#include "dmat.h"

#include "dataStructs.h"
#include "jgtmat.h"
#include "metadata.h"
#include "io.h"

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
/// @details if(has_drop); remove the individual with most dropped pairs from the matrix
/// aim: remove as little individuals as possible while ensuring that there is no missing data in the matrix
dmat_t* dmat_prune_remove_dropped_distances(dmat_t* dmat) {

    BEGIN_LOGSECTION_MSG("(--prune-dmat) Pruning distance matrix");

    if (dmat->n > 1) {
        ERROR("Pruning of distance matrices with more than one matrix is not possible, yet.");
    }
    const size_t ni = 0; // dmat->n - 1

    // -> init
    dmat_t* prunedmat = NULL;
    prunedmat = (dmat_t*)malloc(sizeof(dmat_t));
    ASSERT(prunedmat != NULL);
    prunedmat->n = 0;
    prunedmat->size = 0;
    prunedmat->type = 0;
    prunedmat->transform = 0;
    prunedmat->method = 0;
    prunedmat->names = NULL;
    prunedmat->names_src = 0;
    prunedmat->drop = NULL;
    prunedmat->has_drop = false;

    size_t idx;
    int maxDrop, maxDropIdx, totnPairsToPrune, nItems;
    nItems = dmat->names->len;


    int nDrops[nItems];
    int rmIdx[nItems];
    for (size_t i = 0;i < (size_t) nItems;++i) {
        nDrops[i] = 0;
        rmIdx[i] = -1;
    }

    prunedmat->n = dmat->n;
    // prunedmat->size =  // TBD
    prunedmat->type = dmat->type;
    prunedmat->transform = dmat->transform;
    prunedmat->method = dmat->method;
    // prunedmat->names = // TBD
    prunedmat->names_src = DMAT_NAMES_SRC_PRIVATE;
    prunedmat->matrix = (double**)malloc(prunedmat->n * sizeof(double*));
    ASSERT(prunedmat->matrix != NULL);
    // rest of matrix = // TBD
    // prunedmat->drop =  // NO NEED, NULL above
    // prunedmat->has_drop = false; // NO NEED, FALSE above


    prunedmat->names = strArray_init();

    totnPairsToPrune = 0;
    idx = 0;
    for (size_t i = 1;i < (size_t) nItems;++i) {
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
            for (size_t i = 0;i < (size_t) nItems;++i) {
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
            for (size_t i = 1;i < (size_t) nItems;++i) {
                for (size_t j = 0;j < i;++j) {
                    if (i == (size_t) maxDropIdx) {
                        if (nDrops[j] > 0) {
                            nDrops[j]--;
                            --totnPairsToPrune;
                        }
                    } else if (j == (size_t) maxDropIdx) {
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
        for (size_t i = 0;i < (size_t) nItems;++i) {
            if (nDrops[i] > 0) {
                NEVER;
            }
            bool isRemoved;
            for (size_t k = 0;k < (size_t) nIndsRemoved;++k) {
                if ((int) i == rmIdx[k]) {
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

        kstring_t tmp = KS_INITIALIZE;
        ksprintf(&tmp, "List of removed individuals:\n");
        for (size_t i = 0;i < (size_t) nIndsRemoved;++i) {
            ksprintf(&tmp, "%s\n", dmat->names->d[rmIdx[i]]);
        }
        LOG("%s", tmp.s);
        ks_free(&tmp);

    } else {
        NEVER;
    }


    prunedmat->size = (new_nItems * (new_nItems - 1)) / 2;
    prunedmat->matrix[ni] = (double*)malloc(prunedmat->size * sizeof(double));
    ASSERT(prunedmat->matrix[0] != NULL);
    for (size_t i = 0;i < prunedmat->size;++i) {
        prunedmat->matrix[0][i] = 0.0;
    }

    bool isRemoved;

    size_t newidx = 0;
    for (size_t i = 0;i < (size_t) nItems;++i) {

        isRemoved = false;
        for (size_t k = 0;k < (size_t) nIndsRemoved;++k) {
            if ((int)i == rmIdx[k]) {
                isRemoved = true;
                break;
            }
        }
        if (isRemoved) {
            continue;
        }

        prunedmat->names->add(dmat->names->d[i]);

        for (size_t j = 0;j < i;++j) {

            isRemoved = false;
            for (size_t k = 0;k < (size_t) nIndsRemoved;++k) {
                if ((int)j == rmIdx[k]) {
                    isRemoved = true;
                    break;
                }
            }
            if (isRemoved) {
                continue;
            }

            prunedmat->matrix[ni][newidx] = dmat->matrix[ni][MATRIX_GET_INDEX_LTED_IJ(i, j)];
            ++newidx;
        }
    }

    DEVASSERT(newidx == prunedmat->size);

    if (0 == nIndsRemoved) {
        WARN("(--prune-dmat) No individuals were removed from the distance matrix. The pruned matrix is identical to the original matrix.");

    }

    END_LOGSECTION;

    return(prunedmat);
}


/// @brief dmat_init - initialize a distance matrix data structure
/// @param nInd      - number of individuals
/// @param type      - type of the distance matrix
/// @param method    - method used to calculate the distances in the matrix
/// @param transform - transformation applied to the distances in the matrix
/// @param names     - array of names of the items in the distance matrix
/// @param names_src - source of the names array
/// @return dmat_t* - pointer to the newly allocated and initialized distance matrix
dmat_t* dmat_init(const size_t nInd, const uint8_t type, const uint32_t method, const uint32_t transform, strArray* names, const uint8_t names_src) {

    dmat_t* dmat = NULL;
    dmat = (dmat_t*)malloc(sizeof(dmat_t));
    ASSERT(dmat != NULL);

    // -> init
    dmat->n = 0;
    dmat->size = 0;
    dmat->type = 0;
    dmat->transform = 0;
    dmat->method = 0;
    dmat->names = NULL;
    dmat->names_src = 0;
    dmat->drop = NULL;
    dmat->has_drop = false;


    // -> set 

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

    dmat->n = (args->nBootstraps > 0) ? (1 + args->nBootstraps) : 1;

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

    dmat->drop = (bool**)malloc(dmat->n * sizeof(bool*));
    ASSERT(dmat->drop != NULL);
    for (size_t i = 0; i < dmat->n; ++i) {
        dmat->drop[i] = NULL;
        dmat->drop[i] = (bool*)malloc(dmat->size * sizeof(bool));
        ASSERT(dmat->drop[i] != NULL);
        for (size_t j = 0; j < dmat->size; ++j) {
            dmat->drop[i][j] = false;
        }
    }

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

/// @brief dmat_read - read a distance matrix from a file
/// @param in_dm_fn 
/// @param required_transform 
/// @param metadata 
/// @return 
dmat_t* dmat_read(const char* in_dm_fn, const uint32_t required_transform, metadata_t* metadata) {

    LOG("(--in-dm) Reading distance matrix from file: %s", in_dm_fn);

    FILE* fp = NULL;

    bool isGzFile = false;
    if (IO::isGzFile(args->in_dm_fn) == 1) {
        isGzFile = true;
        kstring_t cmd = KS_INITIALIZE;
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
    dmat->transform = 0;
    dmat->method = 0;
    dmat->names = NULL;
    dmat->names_src = 0;
    dmat->drop = NULL;
    dmat->has_drop = false;

    // -> set

    if (metadata != NULL) {
        dmat->names_src = DMAT_NAMES_SRC_IN_METADATA_NAMES_PTR;
    } else {
        dmat->names_src = DMAT_NAMES_SRC_IN_DM_FILE;
    }

    int lineno = 0;

    BEGIN_LOGSECTION;

    // line 0: number of matrices
    if (fscanf(fp, "%ld", &dmat->n) != 1) {
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
    if (fscanf(fp, "%d", &dmat->transform) != 1) {
        ERROR("Could not read the transform line (line 3) from the distance matrix input file %s", args->in_dm_fn);
    }

    if (DMAT_TRANSFORM_NONE == dmat->transform) {
        LOG("Input distance matrix transform is detected as: None");
    } else if (DMAT_TRANSFORM_SQUARE == dmat->transform) {
        LOG("Input distance matrix transform is detected as: Squared");
    } else {
        ERROR("Unrecognized distance matrix transformation: %d", dmat->transform);
    }


    ++lineno;

    // line 3: method
    if (fscanf(fp, "%d", &dmat->method) != 1) {
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

        if (in_nInd != (int) metadata->indNames->len) {
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
                if ((int) mtdidx != lineno - 6) {
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

    size_t j = 0;
    while ((int) lineno < (int) (7 + ((in_nInd + in_nIndCmb) * dmat->n))) {
        fscanf(fp, "%lf\n", &dmat->matrix[0][j]);
        ++lineno;
        ++j;
    }


    ASSERT((int)dmat->size == in_nIndCmb);
    LOG("Read %ld distances from distance matrix file", dmat->size);

    // -> set drop based on nan values
    dmat->drop = (bool**)malloc(dmat->n * sizeof(bool*));
    ASSERT(dmat->drop != NULL);
    for (size_t i = 0; i < dmat->n; ++i) {
        dmat->drop[i] = NULL;
        dmat->drop[i] = (bool*)malloc(dmat->size * sizeof(bool));
        ASSERT(dmat->drop[i] != NULL);
        for (size_t j = 0; j < dmat->size; ++j) {
            if (isnan(dmat->matrix[i][j])) {
                dmat->drop[i][j] = true;
                dmat->has_drop = true;
            } else {
                dmat->drop[i][j] = false;
            }
        }
    }


    if (required_transform != dmat->transform) {

        if (required_transform == DMAT_TRANSFORM_NONE) {
            // -> we need: no transform

            if (dmat->transform == DMAT_TRANSFORM_SQUARE) {
                // -> we have: squared
                LOG("Input distance matrix transform (Squared, %d) is not the same as the required transform (None, %d). Will take square roots of the input values to obtain the required transform.", dmat->transform, required_transform);
                for (size_t m = 0;m < dmat->n;++m) {
                    for (size_t i = 0;i < dmat->size;++i) {
                        dmat->matrix[m][i] = sqrt(dmat->matrix[m][i]);
                    }
                }
                dmat->transform = DMAT_TRANSFORM_NONE;
            } else {
                NEVER;
            }

        } else if (required_transform == DMAT_TRANSFORM_SQUARE) {
            // -> we need: squared

            if (dmat->transform == DMAT_TRANSFORM_NONE) {
                // -> we have: no transform
                LOG("Input distance matrix transform (None, %d) is not the same as the required transform (Squared, %d). Will take squares of the input values to obtain the required transform.", dmat->transform, required_transform);
                for (size_t m = 0;m < dmat->n;++m) {
                    for (size_t i = 0;i < dmat->size;++i) {
                        if (!dmat->drop) {
                            dmat->matrix[m][i] *= dmat->matrix[m][i];
                        }
                    }
                }
                dmat->transform = DMAT_TRANSFORM_SQUARE;
            } else {
                NEVER;
            }
        } else {
            NEVER;
        }
    }

    if (isGzFile) {
        PCLOSE(fp);
    } else {
        FCLOSE(fp);

    }


    // -> finished reading the distance matrix file

    END_LOGSECTION;
    return(dmat);
}


/// @brief dmat_print - print a distance matrix to a kstring
/// @param dmat - pointer to the distance matrix to be printed
/// @param kbuf - pointer to the kstring to which the distance matrix will be printed
void dmat_print(dmat_t* dmat, outfile_t* outfile) {

    LOG("(--print-dm) Printing distance matrix to: %s", outfile->fn);

    kstring_t* kbuf = &outfile->kbuf;

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

    // lines 6-X: dmat->names
    for (size_t i = 0;i < dmat->names->len;++i) {
        ksprintf(kbuf, "%s\n", dmat->names->d[i]);
    }

    double* matrix = NULL;

    for (size_t mi = 0;mi < dmat->n;++mi) {
        matrix = dmat->matrix[mi];

        // lines (X+1)-Y: distances 
        for (size_t p = 0; p < dmat->size; ++p) {
            ksprintf(kbuf, "%.17g\n", matrix[p]);
        }
    }

    return;
}

void dmat_print_verbose(dmat_t* dmat, outfile_t* outfile){

    if (dmat->type != DMAT_TYPE_LTED) {
        ERROR("Verbose printing is only supported for Lower Triangular Matrix (Excluding Diagonal) distance matrices.");
    }

    kstring_t* kbuf = &outfile->kbuf;

    // line 0: number of matrices
    ksprintf(kbuf, "# This file was produced by: %s\n", PROGRAM_VERSION_INFO);
    // ksprintf(kbuf, "# The command was: %s\n", PROGRAM_COMMAND); // TODO
    ksprintf(kbuf, "# Number of distance matrices: %ld\n", dmat->n);

    // line 1: type
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

    // line 2: transformation
    ksprintf(kbuf, "# Transformation applied to the distances: ");
    if (DMAT_TRANSFORM_NONE == dmat->transform) {
        ksprintf(kbuf, "%s", "None");
    } else if (DMAT_TRANSFORM_SQUARE == dmat->transform) {
        ksprintf(kbuf, "%s", "Squared");
    } else {
        NEVER;
    }

    ksprintf(kbuf, " (--dm-transform %d)", DMAT_TRANSFORM_NONE);
    kputc('\n', kbuf);


    // line 3: method
    ksprintf(kbuf, "# Method used to calculate the distances: ");
    ksprintf(kbuf, "%s", get_dmat_method_str(dmat->method));
    ksprintf(kbuf, " (--dm-method %d)", DMAT_METHOD_DIJ);
    kputc('\n', kbuf);

    // line 4: number of dmat->names
    ksprintf(kbuf, "# Number of items in the distance matrix: %ld\n", dmat->names->len);

    // line 5: number of distances 
    ksprintf(kbuf, "# Number of distances in the distance matrix: %ld\n", dmat->size);

    // lines 6-X: dmat->names
    ksprintf(kbuf, "# NAME, Names of the items in the distance matrix:\n");
    ksprintf(kbuf, "# NAME\t[2]item_index\t[3]item_name\n");
    for (size_t i = 0;i < dmat->names->len;++i) {
        // ksprintf(kbuf, "%s\n", dmat->names->d[i]);
        ksprintf(kbuf, "NAME\t%ld\t%s\n", i, dmat->names->d[i]);
    }

    double* matrix = NULL;

    ksprintf(kbuf, "# DIST, Pairwise distances in the distance matrix:\n");
    ksprintf(kbuf, "# DIST\t[2]matrix_index\t[3]pair_index\t[4]first_item_index\t[5]second_item_index\t[6]first_item_name\t[7]second_item_name\t[8]distance\n");

    for (size_t mi = 0;mi < dmat->n;++mi) {
        matrix = dmat->matrix[mi];


        size_t idx;
        idx = 0;
        for (size_t i = 0;i < dmat->names->len;++i) {
            for (size_t j = 0;j < i;++j) {
                ksprintf(kbuf, "DIST\t%ld\t%ld\t%ld\t%ld\t%s\t%s\t%.17g\n", mi, idx, i, j, dmat->names->d[i], dmat->names->d[j], matrix[idx]);
                ++idx;
            }
        }

        // lines (X+1)-Y: distances 
        // for (size_t p = 0; p < dmat->size; ++p) {
            // ksprintf(kbuf, "%.17g\n", matrix[p]);

        // }
    // print as matrix (pretty)
        ksprintf(kbuf, "# DMAT, Distance matrix in matrix format:\n");
        ksprintf(kbuf, "%-25s\t%-25s\t", "DMAT", "ITEMS");
        for (size_t i = 0;i < dmat->names->len;++i) {
            ksprintf(kbuf, "%-25s", dmat->names->d[i]);
            if (i != dmat->names->len - 1) {
                kputc('\t', kbuf);
            }
        }
        kputc('\n', kbuf);
        idx = 0;
        for (size_t i = 0;i < dmat->names->len;++i) {
            ksprintf(kbuf, "%-25s\t", "DMAT");
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
    }

    return;
}


/// @brief dmat_calculate_distances - calculate the distances using jgtmat and store them in dmat
/// @param jgtmat - pointer to the joint genotype matrix data structure to be used for calculating the distances
/// @param dmat   - pointer to the distance matrix data structure to store the calculated distances
void dmat_calculate_distances(jgtmat_t* jgtmat, dmat_t* dmat) {

    double* dm = NULL;

    uint32_t transform = dmat->transform;
    uint32_t method = dmat->method;
    const size_t nPairs = dmat->size;

    bool* drop = jgtmat->drop;

    int nPairs_after_drop;
    size_t idx;

    for (size_t i = 0;i < dmat->n;++i) {
        // N.B. i==0 original run in dmat

        nPairs_after_drop = nPairs;

        idx = i * nPairs;
        dm = dmat->matrix[i];

        if (PROGRAM_WILL_USE_BCF_FMT_GT) {

            double sum, x, A, B, C, D, E, F, G, H, I;
            uint64_t** m = jgtmat->m;
            uint64_t* pm = NULL;

            for (size_t p = 0; p < nPairs; ++p) {

                if (drop != NULL && drop[idx + p]) {
                    dmat->drop[i][p] = true;
                    dmat->has_drop = true;
                    dm[p] = NAN;
                    nPairs_after_drop--;
                    continue;
                } else {
                    dmat->drop[i][p] = false;
                }

                pm = m[idx + p];
                A = pm[JGTMAT_LINEAR_IDXOF_A];
                B = pm[JGTMAT_LINEAR_IDXOF_B];
                C = pm[JGTMAT_LINEAR_IDXOF_C];
                D = pm[JGTMAT_LINEAR_IDXOF_D];
                E = pm[JGTMAT_LINEAR_IDXOF_E];
                F = pm[JGTMAT_LINEAR_IDXOF_F];
                G = pm[JGTMAT_LINEAR_IDXOF_G];
                H = pm[JGTMAT_LINEAR_IDXOF_H];
                I = pm[JGTMAT_LINEAR_IDXOF_I];
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


                if (DMAT_TRANSFORM_NONE != transform) {
                    if (DMAT_TRANSFORM_SQUARE == transform) {
                        x *= x;
                    }
                }
                dm[p] = x;
            }

        } else if (PROGRAM_WILL_USE_BCF_FMT_GL) {

            double x, A, B, C, D, E, F, G, H, I;
            double** m = jgtmat->pm;
            double* pm = NULL; // pair pm

            for (size_t p = 0; p < nPairs; ++p) {

                if (drop != NULL && drop[idx + p]) {
                    dmat->drop[i][p] = true;
                    dmat->has_drop = true;
                    dm[p] = NAN;
                    nPairs_after_drop--;
                    continue;
                } else {
                    dmat->drop[i][p] = false;
                }

                pm = m[idx + p];
                A = pm[JGTMAT_LINEAR_IDXOF_A];
                B = pm[JGTMAT_LINEAR_IDXOF_B];
                C = pm[JGTMAT_LINEAR_IDXOF_C];
                D = pm[JGTMAT_LINEAR_IDXOF_D];
                E = pm[JGTMAT_LINEAR_IDXOF_E];
                F = pm[JGTMAT_LINEAR_IDXOF_F];
                G = pm[JGTMAT_LINEAR_IDXOF_G];
                H = pm[JGTMAT_LINEAR_IDXOF_H];
                I = pm[JGTMAT_LINEAR_IDXOF_I];


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

                if (DMAT_TRANSFORM_NONE != transform) {
                    if (DMAT_TRANSFORM_SQUARE == transform) {
                        x *= x;
                    }
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

    return;
}

