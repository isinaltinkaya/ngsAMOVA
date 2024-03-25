#include "dmat.h"
#include "dataStructs.h"
#include "jgtmat.h"

dmat_t* dmat_init(const size_t nInd, const uint8_t type, const uint32_t method, const uint32_t transform, strArray* names, const uint8_t names_src) {

    dmat_t* ret = NULL;
    ret = (dmat_t*)malloc(sizeof(dmat_t));
    ASSERT(ret != NULL);

    switch (type) {
    case DMAT_TYPE_LTED:
    case DMAT_TYPE_UTED:
        ret->size = (nInd * (nInd - 1)) / 2;
        break;
    case DMAT_TYPE_LTID:
    case DMAT_TYPE_UTID:
        ret->size = (nInd * (nInd + 1)) / 2;
        break;
    case DMAT_TYPE_FULL:
        ret->size = nInd * nInd;
        break;
    default:
        NEVER;
    }

    ret->type = type;
    ret->transform = transform;
    ret->method = method;

    ret->n = (args->nBootstraps > 0) ? (1 + args->nBootstraps) : 1;

    ret->names_src = names_src;
    if (ret->names_src == DMAT_NAMES_SRC_IN_VCF_PARS_PTR) {
        ret->names = names;
    } else if (ret->names_src == DMAT_NAMES_SRC_IN_METADATA_NAMES_PTR) {
        ret->names = names;
    } else if (ret->names_src == DMAT_NAMES_SRC_PRIVATE) {
        // names is allocated and used internally in the program
    } else {
        NEVER;
    }

    ret->matrix = NULL;
    ret->matrix = (double**)malloc(ret->n * sizeof(dmat_t*));
    ASSERT(ret->matrix != NULL);

    for (size_t i = 0; i < ret->n;++i) {
        ret->matrix[i] = NULL;
        ret->matrix[i] = (double*)malloc(ret->size * sizeof(double));
        ASSERT(ret->matrix[i] != NULL);
        for (size_t j = 0;j < ret->size;++j) {
            ret->matrix[i][j] = 0.0;
        }
    }
    return(ret);
}

void dmat_destroy(dmat_t* d) {
    for (size_t i = 0; i < d->n;++i) {
        FREE(d->matrix[i]);
    }
    FREE(d->matrix);

    if (d->names_src == DMAT_NAMES_SRC_IN_DM_FILE) {
        strArray_destroy(d->names);
    } else if (d->names_src == DMAT_NAMES_SRC_IN_VCF_PARS_PTR) {
        d->names = NULL;
    } else if (d->names_src == DMAT_NAMES_SRC_IN_METADATA_NAMES_PTR) {
        d->names = NULL;
    } else if (d->names_src == DMAT_NAMES_SRC_NONE) {
        NEVER;
    } else if (d->names_src == DMAT_NAMES_SRC_PRIVATE) {
        strArray_destroy(d->names);
    } else {
        NEVER;
    }

    FREE(d);
    return;
}

dmat_t* dmat_read(const char* in_dm_fn, const uint32_t required_transform, metadataStruct* metadata) {

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

    dmat_t* ret = NULL;
    ret = (dmat_t*)malloc(sizeof(dmat_t));
    ASSERT(ret != NULL);
    ret->type = 0;
    ret->transform = 0;
    ret->method = 0;
    ret->size = 0;
    ret->names = NULL;

    if (metadata != NULL) {
        ret->names_src = DMAT_NAMES_SRC_IN_METADATA_NAMES_PTR;
    } else {
        ret->names_src = DMAT_NAMES_SRC_IN_DM_FILE;
    }
    ret->matrix = NULL;
    ret->n = 0;

    int lineno = 0;

    BEGIN_LOGSECTION;

    // line 0: number of matrices
    if (fscanf(fp, "%ld", &ret->n) != 1) {
        // ERROR()
        ERROR("Could not read the number of matrices line (line 1) from the distance matrix input file %s", args->in_dm_fn);
    }
    if (0 == ret->n) {
        ERROR("Distance matrix input file must contain at least 1 distance matrices (line 1 > 0)");
    }

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
        ret->type = DMAT_TYPE_FULL;
        LOG("Input distance matrix type is detected as FULL: Full Matrix");
    } else if (in_type[0] == 'L' && in_type[1] == 'T' && in_type[2] == 'E' && in_type[3] == 'D') {
        ret->type = DMAT_TYPE_LTED;
        LOG("Input distance matrix type is detected as LTED: Lower Triangular Matrix (Excluding Diagonal)");
    } else if (in_type[0] == 'L' && in_type[1] == 'T' && in_type[2] == 'I' && in_type[3] == 'D') {
        ret->type = DMAT_TYPE_LTID;
        LOG("Input distance matrix type is detected as LTID: Lower Triangular Matrix (Including Diagonal)");
    } else if (in_type[0] == 'U' && in_type[1] == 'T' && in_type[2] == 'E' && in_type[3] == 'D') {
        ret->type = DMAT_TYPE_UTED;
        LOG("Input distance matrix type is detected as UTED: Upper Triangular Matrix (Excluding Diagonal)");
    } else if (in_type[0] == 'U' && in_type[1] == 'T' && in_type[2] == 'I' && in_type[3] == 'D') {
        ret->type = DMAT_TYPE_UTID;
        LOG("Input distance matrix type is detected as UTID: Upper Triangular Matrix (Including Diagonal)");
    } else {
        ERROR("Unrecognized distance matrix type: %s", in_type);
    }

    ++lineno;

    // line 2: transformation
    if (fscanf(fp, "%d", &ret->transform) != 1) {
        ERROR("Could not read the transform line (line 3) from the distance matrix input file %s", args->in_dm_fn);
    }

    if (DMAT_TRANSFORM_NONE == ret->transform) {
        LOG("Input distance matrix transform is detected as: None");
    } else if (DMAT_TRANSFORM_SQUARE == ret->transform) {
        LOG("Input distance matrix transform is detected as: Squared");
    } else {
        ERROR("Unrecognized distance matrix transformation: %d", ret->transform);
    }


    ++lineno;

    // line 3: method
    if (fscanf(fp, "%d", &ret->method) != 1) {
        ERROR("Could not read the method line (line 4) from the distance matrix input file %s", args->in_dm_fn);
    }

    if (DMAT_METHOD_DIJ == ret->method) {
        LOG("Input distance matrix method is detected as: Dij");
    } else if (DMAT_METHOD_FIJ == ret->method) {
        LOG("Input distance matrix method is detected as: Fij");
    } else {
        ERROR("Unrecognized distance matrix method: %d", ret->method);
    }


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

    ret->size = in_nIndCmb;

    if (ret->type == DMAT_TYPE_LTED) {
        if ((int)ret->size != ((in_nInd * (in_nInd - 1)) / 2)) {
            ERROR("Number of distances in the input distance matrix (%ld) does not match the expected number of pairwise distances given the number of individuals (%d) and the distance matrix type (LTED).", ret->size, in_nInd);
        }
    } else if (ret->type == DMAT_TYPE_UTED) {
        if ((int)ret->size != ((in_nInd * (in_nInd - 1)) / 2)) {
            ERROR("Number of distances in the input distance matrix (%ld) does not match the expected number of pairwise distances given the number of individuals (%d) and the distance matrix type (UTED).", ret->size, in_nInd);
        }
    } else if (ret->type == DMAT_TYPE_LTID) {
        if ((int)ret->size != ((in_nInd * (in_nInd + 1)) / 2)) {
            ERROR("Number of distances in the input distance matrix (%ld) does not match the expected number of pairwise distances given the number of individuals (%d) and the distance matrix type (LTID).", ret->size, in_nInd);
        }
    } else if (ret->type == DMAT_TYPE_UTID) {
        if ((int)ret->size != ((in_nInd * (in_nInd + 1)) / 2)) {
            ERROR("Number of distances in the input distance matrix (%ld) does not match the expected number of pairwise distances given the number of individuals (%d) and the distance matrix type (UTID).", ret->size, in_nInd);
        }
    } else if (ret->type == DMAT_TYPE_FULL) {
        if ((int)ret->size != (in_nInd * in_nInd)) {
            ERROR("Number of distances in the input distance matrix (%ld) does not match the expected number of pairwise distances given the number of individuals (%d) and the distance matrix type (FULL).", ret->size, in_nInd);
        }
    } else {
        NEVER;
    }

    ++lineno;


    /// ------------------------------------------------------------
    // lines 6:(6+nInd-1) names

    if (ret->names_src == DMAT_NAMES_SRC_IN_DM_FILE) {
        ret->names = strArray_alloc(in_nInd);

        char tmp[1024] = { '\0' };
        while (lineno < 6 + in_nInd) {
            tmp[0] = '\0';
            fscanf(fp, "%s\n", tmp);
            if ('\0' == tmp[0]) {
                ERROR("Found bad name at line %d of distance matrix file %s", lineno + 1, args->in_dm_fn);
            }
            ret->names->add(tmp);
            ++lineno;
        }

        ASSERT((int)ret->names->len == in_nInd);
        LOG("Read %ld names from distance matrix file", ret->names->len);

        ++lineno;

    } else if (ret->names_src == DMAT_NAMES_SRC_IN_METADATA_NAMES_PTR) {

        if (in_nInd != metadata->indNames->len) {
            ERROR("Number of names in the distance matrix file (%ld) does not match the number of individuals in the metadata file (%ld). Please edit your metadata file to match the names in the distance matrix file.", in_nInd, metadata->indNames->len);
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
                if (mtdidx != lineno - 6) {
                    DEVPRINT("name:%s idx:%ld mtdidx:%ld mtdname:%s", tmp, lineno - 6, mtdidx, metadata->indNames->d[mtdidx]);
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

        ret->names = metadata->indNames;

    } else {
        NEVER;
    }

    /// ------------------------------------------------------------
    // lines (6+nInd):(6+nInd+nIndCmb-1) distances

    ret->matrix = (double**)malloc(ret->n * sizeof(dmat_t*));
    ASSERT(ret->matrix != NULL);

    for (size_t i = 0; i < ret->n;++i) {
        ret->matrix[i] = NULL;
        ret->matrix[i] = (double*)malloc(ret->size * sizeof(double));
        ASSERT(ret->matrix[i] != NULL);
        for (size_t j = 0;j < ret->size;++j) {
            ret->matrix[i][j] = 0.0;
        }
    }

    size_t j = 0;
    if (ret->n == 1) {
        while (lineno < 7 + in_nInd + in_nIndCmb) {
            fscanf(fp, "%lf\n", &ret->matrix[0][j]);
            ++lineno;
            ++j;
        }
    } else {
        NEVER;//TODO
    }


    ASSERT((int)ret->size == in_nIndCmb);
    LOG("Read %ld distances from distance matrix file", ret->size);



    if (required_transform != ret->transform) {

        if (required_transform == DMAT_TRANSFORM_NONE) {
            // -> we need: no transform

            if (ret->transform == DMAT_TRANSFORM_SQUARE) {
                // -> we have: squared
                LOG("Input distance matrix transform (Squared, %d) is not the same as the required transform (None, %d). Will take square roots of the input values to obtain the required transform.", ret->transform, required_transform);
                for (size_t m = 0;m < ret->n;++m) {
                    for (size_t i = 0;i < ret->size;++i) {
                        ret->matrix[m][i] = sqrt(ret->matrix[m][i]);
                    }
                }
                ret->transform = DMAT_TRANSFORM_NONE;
            } else {
                NEVER;
            }

        } else if (required_transform == DMAT_TRANSFORM_SQUARE) {
            // -> we need: squared

            if (ret->transform == DMAT_TRANSFORM_NONE) {
                // -> we have: no transform
                LOG("Input distance matrix transform (None, %d) is not the same as the required transform (Squared, %d). Will take squares of the input values to obtain the required transform.", ret->transform, required_transform);
                for (size_t m = 0;m < ret->n;++m) {
                    for (size_t i = 0;i < ret->size;++i) {
                        ret->matrix[m][i] *= ret->matrix[m][i];
                    }
                }
                ret->transform = DMAT_TRANSFORM_SQUARE;
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
    return(ret);
}

void dmat_write(dmat_t* dmat) {
    LOG("Writing distance matrix to %s\n", outFiles->out_dm_fs->fn);
    outFiles->out_dm_fs->kbuf = kbuf_init();
    dmat_print(dmat, outFiles->out_dm_fs->kbuf);
    outFiles->out_dm_fs->kbuf_write();
}


void dmat_print(dmat_t* dmat, kstring_t* kstr) {

    // line 0: number of matrices
    ksprintf(kstr, "%ld\n", dmat->n);

    // line 1: type
    if (DMAT_TYPE_LTED == dmat->type) {
        ksprintf(kstr, "%s\n", "LTED");
    } else if (DMAT_TYPE_LTID == dmat->type) {
        ksprintf(kstr, "%s\n", "LTID");
    } else if (DMAT_TYPE_UTED == dmat->type) {
        ksprintf(kstr, "%s\n", "UTED");
    } else if (DMAT_TYPE_UTID == dmat->type) {
        ksprintf(kstr, "%s\n", "UTID");
    } else if (DMAT_TYPE_FULL == dmat->type) {
        ksprintf(kstr, "%s\n", "FULL");
    } else {
        NEVER;
    }

    // line 2: transformation
    ksprintf(kstr, "%d\n", dmat->transform);

    // line 3: method
    ksprintf(kstr, "%d\n", dmat->method);

    // line 4: number of dmat->names
    ksprintf(kstr, "%ld\n", dmat->names->len);

    // line 5: number of distances 
    ksprintf(kstr, "%ld\n", dmat->size);

    // lines 6-X: dmat->names
    for (size_t i = 0;i < dmat->names->len;++i) {
        ksprintf(kstr, "%s\n", dmat->names->d[i]);
    }

    double* matrix = NULL;

    for(size_t mi=0;mi<dmat->n;++mi){
        matrix=dmat->matrix[mi];

        // lines (X+1)-Y: distances 
        for (size_t p = 0; p < dmat->size; ++p) {
            // ksprintf(kstr, "%f\n", matrix[p]);
            ksprintf(kstr,"%.17g\n",matrix[p]);
        }
    }

    return;
}

void dmat_calculate_distances(jgtmat_t* jgtmat, dmat_t* dmat) {


    double* dm = NULL;

    uint32_t transform = dmat->transform;
    uint32_t method= dmat->method;
    const size_t nPairs = dmat->size;


    size_t idx;
    for (size_t i = 0;i < dmat->n;++i) {
        // N.B. i==0 original run in dmat

        idx = i * nPairs;
        dm = dmat->matrix[i];

        if (PROGRAM_WILL_USE_BCF_FMT_GT) {

            double sum, x, A, B, C, D, E, F, G, H, I;
            uint64_t** m = jgtmat->m;
            uint64_t* pm = NULL; 

            for (size_t p = 0; p < nPairs; ++p) {

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

    }


    return;
}