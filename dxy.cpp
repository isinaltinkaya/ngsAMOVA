#include "dxy.h"
#include "dataStructs.h"
#include "dmat.h"

void dxy_destroy(dxy_t* dxy) {
    FREE(dxy->d);
    FREE(dxy->g1names_p);
    FREE(dxy->g2names_p);
    FREE(dxy->levelnames_p);
    for (size_t i = 0;i < dxy->nLevels;++i) {
        if (dxy->dm[i] != NULL) {
            dmat_destroy(dxy->dm[i]);
        }
    }
    FREE(dxy->dm);
    FREE(dxy);
    return;
}

 dxy_t* dxy_init(metadataStruct* mtd) {
    dxy_t* dxy = NULL;
    dxy = (dxy_t*)malloc(sizeof(dxy_t));
    ASSERT(dxy != NULL);

    // -> init
    dxy->n = 0;
    dxy->d = NULL;
    dxy->dm = NULL;

    // -> set
    dxy->nLevels = mtd->nLevels - 1; // excluding the "individual" level

    dxy->dm = (dmat_t**)malloc(dxy->nLevels * sizeof(dmat_t*));
    ASSERT(dxy->dm != NULL);

    size_t ndxy = 0; // number of pairwise group comparisons to be made
    for (size_t lvl = 0; lvl < dxy->nLevels; ++lvl) {
        size_t n_groups_at_level = mtd->level2groupIndices[lvl]->len;
        if (n_groups_at_level == 1) {
            dxy->dm[lvl] = NULL;
            continue;
        }
        ndxy += (n_groups_at_level * (n_groups_at_level - 1)) / 2;
        dxy->dm[lvl] = dmat_init(n_groups_at_level, DMAT_TYPE_LTED, DMAT_METHOD_DXY, DMAT_TRANSFORM_NONE, NULL, DMAT_NAMES_SRC_PRIVATE);
        dxy->dm[lvl]->names = strArray_alloc(n_groups_at_level);
        for (size_t g = 0; g < n_groups_at_level; ++g) {
            dxy->dm[lvl]->names->add(mtd->groupNames->d[mtd->level2groupIndices[lvl]->d[g]]);
        }
    }

    if (ndxy == 0) {
        ERROR("Could not find any pairwise group comparisons to be made. Please check the metadata file.");
    }

    dxy->n = ndxy;

    dxy->d = (double*)malloc(dxy->n * sizeof(double));
    ASSERT(dxy->d != NULL);

    dxy->g1names_p = (char**)malloc(dxy->n * sizeof(char*));
    ASSERT(dxy->g1names_p != NULL);

    dxy->g2names_p = (char**)malloc(dxy->n * sizeof(char*));
    ASSERT(dxy->g2names_p != NULL);

    dxy->levelnames_p = (char**)malloc(dxy->n * sizeof(char*));
    ASSERT(dxy->levelnames_p != NULL);

    return(dxy);
}


// dxy file format: comma-separated list of pairwise dxy values
// group1Name,group2Name,levelID,dxyValue
dxy_t* dxy_read(paramStruct* pars, dmat_t* dmat, metadataStruct* mtd) {
    dxy_t* dxy = new dxy_t();

    // number of lines in dxy file == number of pairwise dxy values
    int n_vals = 0;

    char* line = (char*)malloc(FGETS_BUF_SIZE);
    ASSERT(line != NULL);

    char dxy_buf[FGETS_BUF_SIZE];

    FILE* in_dxy_fp = fopen(args->in_dxy_fn, "r");

    // skip the first line (header)
    ASSERT(fgets(dxy_buf, FGETS_BUF_SIZE, in_dxy_fp) != NULL); // error or unexpected eof; so handle both with assert

    int col = 0;
    while (fgets(dxy_buf, FGETS_BUF_SIZE, in_dxy_fp)) {
        col = 0;
        char* tok = strtok(dxy_buf, ",");
        while (tok != NULL) {
            // while (n_vals >= (int)dxy->_d) {
                // dxy->expand();
                //TODO
            // }

            switch (col) {
            case 0:
                // group1's name
                dxy->g1names_p[n_vals] = strdup(tok);
                break;
            case 1:
                // group2's name
                dxy->g2names_p[n_vals] = strdup(tok);
                break;
            case 2:
                // levelName
                dxy->levelnames_p[n_vals] = strdup(tok);
                break;
            case 3:
                // dxyValue
                dxy->d[n_vals] = atof(tok);
                break;
            default:
                fprintf(stderr, "\n[ERROR][dxy_read]\tToo many columns in dxy file %s.\n", args->in_dxy_fn);
                exit(1);
            }

            tok = strtok(NULL, ",");
            ++col;
        }
        ++n_vals;
    }

    dxy->n = n_vals;

    FREE(line);
    FCLOSE(in_dxy_fp);

    return dxy;
}


void dxy_print(dxy_t* dxy) {
    fprintf(stderr, "\n[INFO]\t-> Writing the dxy results to %s.\n", outFiles->out_dxy_fs->fn);
    outFiles->out_dxy_fs->kbuf = kbuf_init();
    ksprintf(outFiles->out_dxy_fs->kbuf, "hierarchy_level,group1_id,group2_id,dxy\n");
    for (size_t i = 0; i < dxy->n; i++) {
        ksprintf(outFiles->out_dxy_fs->kbuf, "%s,%s,%s,%f\n", dxy->levelnames_p[i], dxy->g1names_p[i], dxy->g2names_p[i], dxy->d[i]);
    }
    outFiles->out_dxy_fs->kbuf_write();


    double* matrix = NULL;
    dmat_t* dmi = NULL;
    IO::outputStruct* out_dxy_i = NULL;
    for (size_t i = 0;i < dxy->nLevels;++i) {

        dmi = dxy->dm[i];
        if (dmi == NULL) {
            continue;
        }

        // suffix = "levelName_dxy.csv"
        char* suffix = (char*)malloc(strlen(dxy->levelnames_p[i]) + 10);
        ASSERT(suffix != NULL);
        sprintf(suffix, "_%s_dxy.csv", dxy->levelnames_p[i]);
        out_dxy_i = new IO::outputStruct(args->out_fnp, suffix, OUTFC::NONE);
        LOG("Writing DXY distance matrix to %s\n", out_dxy_i->fn);
        out_dxy_i->kbuf = kbuf_init();
        dmat_print(dmi, out_dxy_i->kbuf);
        out_dxy_i->kbuf_write();
        delete(out_dxy_i);
        FREE(suffix);
    }


    return;
}

dxy_t* dxy_get(paramStruct* pars, dmat_t* dmat, metadataStruct* mtd) {

    dxy_t* dxy = dxy_init(mtd);

    double* dm = NULL;
    if (dmat->n == 1) {
        dm = dmat->matrix[dmat->n - 1];
    } else {
        NEVER;//TODO
    }


    size_t dxy_i = 0;
    double val = 0.0;

    for (size_t lvl = 0; lvl < dxy->nLevels; ++lvl) {

        size_t n_groups_at_level = mtd->level2groupIndices[lvl]->len;
        if (n_groups_at_level == 1) {
            LOG("[DXY] Skipping level %ld. Reason: only one group found at this level.\n", lvl);
            continue;
        }

        // -> perform dxy for each pair of groups at this level

        // -> only use the individual pairs in the distance matrix where one individual is from group 1 and the other is from group 2
        for (size_t g1 = 1; g1 < n_groups_at_level; ++g1) {
            for (size_t g2 = 0; g2 < g1; ++g2) {

                size_t g1id = mtd->level2groupIndices[lvl]->d[g1];
                size_t g2id = mtd->level2groupIndices[lvl]->d[g2];

                size_t Ng1 = mtd->group2indIndices[g1id]->len;
                size_t Ng2 = mtd->group2indIndices[g2id]->len;

                val = 0.0;
                for (size_t i1 = 1; i1 < Ng1; ++i1) {
                    for (size_t i2 = 0; i2 < Ng2; ++i2) {
                        size_t i1id = mtd->group2indIndices[g1id]->d[i1];
                        size_t i2id = mtd->group2indIndices[g2id]->d[i2];

                        double dist = dm[MATRIX_GET_INDEX_LTED_IJ_UNORDERED(i1id, i2id)];
                        val += dist;
                    }
                }
                val = val / (double)(Ng1 * Ng2);
                dxy->d[dxy_i] = val;
                dxy->g1names_p[dxy_i] = mtd->groupNames->d[g1id];
                dxy->g2names_p[dxy_i] = mtd->groupNames->d[g2id];
                dxy->levelnames_p[dxy_i] = mtd->levelNames->d[lvl];
                ++dxy_i;
            }
        }

    }

    DEVASSERT(dxy_i == dxy->n); // sanity check

    dxy_print(dxy);

    return(dxy);
}
