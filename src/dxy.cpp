#include "dxy.h"

#include "dataStructs.h"
#include "metadata.h"


// TODO add block bootstrapping for dxy

void dxy_destroy(dxy_t* dxy) {
    FREE(dxy->d);
    FREE(dxy->g1names_p);
    FREE(dxy->g2names_p);
    FREE(dxy->levelnames_p);
    for (size_t i = 0;i < dxy->nLevels;++i) {
        if (dxy->dmat[i] != NULL) {
            dmat_destroy(dxy->dmat[i]);
        }
    }
    FREE(dxy->dmat);
    FREE(dxy);
    return;
}

static dxy_t* dxy_init(metadata_t* mtd) {
    dxy_t* dxy = NULL;
    dxy = (dxy_t*)malloc(sizeof(dxy_t));
    ASSERT(dxy != NULL);

    // -> init
    dxy->n = 0;
    dxy->d = NULL;
    dxy->dmat = NULL;

    // -> set
    dxy->nLevels = mtd->nLevels - 1; // excluding the "individual" level

    dxy->dmat = (dmat_t**)malloc(dxy->nLevels * sizeof(dmat_t*));
    ASSERT(dxy->dmat != NULL);

    size_t ndxy = 0; // number of pairwise group comparisons to be made
    for (size_t lvl = 0; lvl < dxy->nLevels; ++lvl) {
        size_t n_groups_at_level = mtd->level2groupIndices[lvl]->len;
        if (n_groups_at_level == 1) {
            dxy->dmat[lvl] = NULL;
            continue;
        }
        ndxy += (n_groups_at_level * (n_groups_at_level - 1)) / 2;
        dxy->dmat[lvl] = dmat_init(1, n_groups_at_level, DMAT_TYPE_LTED, DMAT_METHOD_DXY, DMAT_INTPLUS_TRANSFORM_NONE, NULL, DMAT_NAMES_SRC_PRIVATE);
        dxy->dmat[lvl]->names = strArray_alloc(n_groups_at_level);
        for (size_t g = 0; g < n_groups_at_level; ++g) {
            dxy->dmat[lvl]->names->add(mtd->groupNames->d[mtd->level2groupIndices[lvl]->d[g]]);
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



static void dxy_print(dxy_t* dxy, outfile_t* outfile) {
    fprintf(stderr, "\n[INFO]\t-> Writing the dxy results to %s.\n", outfile->fn);
    kstring_t* kbuf = &outfile->kbuf;
    ksprintf(kbuf, "hierarchy_level,group1_id,group2_id,dxy\n");
    for (size_t i = 0; i < dxy->n; i++) {
        ksprintf(kbuf, "%s,%s,%s,%.17g\n", dxy->levelnames_p[i], dxy->g1names_p[i], dxy->g2names_p[i], dxy->d[i]);
    }

    return;
}

dxy_t* dxy_get(dmat_t* dmat, metadata_t* mtd) {

    if (dmat->drop != NULL) {
        ERROR("Distance matrix with dropped pairs is not supported in dXY. Please prune the distance matrix via --prune-dmat option.");
    }


    dxy_t* dxy = dxy_init(mtd);

    double* dm = NULL;

    if (dmat->n == 1) {
        dm = dmat->matrix[0];
    } else {
        NEVER;//TODO
    }


    double dist;
    size_t i1id, i2id;
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
                        i1id = mtd->group2indIndices[g1id]->d[i1];
                        i2id = mtd->group2indIndices[g2id]->d[i2];
                        dist = dm[MATRIX_GET_INDEX_LTED_IJ_UNORDERED(i1id, i2id)];
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

    if (args->print_dxy) {
        outfile_t* outfile = outfile_init("dxy", "csv", args->print_dxy_ctype);
        dxy_print(dxy, outfile);
        outfile_write(outfile);
        outfile_destroy(outfile);
    }

    return(dxy);
}
