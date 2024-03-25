#include "amova.h"
#include "dmat.h"



/// @brief calculate part of AMOVA statistics where the calculations are independent of the distance matrix but depend on the metadata
/// @param amova amovaStruct
/// @param mtd metadataStruct
/// @return void
/// @details
/// perform shared AMOVA calculations depends on the metadata but does not depend on the distance matrix
/// therefore this function can be called once for each metadataStruct
/// i.e. one call is enough for all amova replicates
/// 
///     
void amova_run_shared(amovaStruct* amova) {

    metadataStruct* mtd = amova->metadata;
    const int nInd = mtd->nInd;

    const size_t nLevels = (size_t)mtd->nLevels; // L
    size_t lvlidx; // 0-based lvl

    /// ------------------------------------------------------------------- ///
    /// -> DEGREES OF FREEDOM (amova->df)
    /// @note 
    ///
    /// df_i = k_i - k_{i-1}  , 0 < i < L
    /// df_total = N - 1
    /// k_0 = 1
    /// k_i = # groups at level i
    /// k_L = N
    ///
    /// 0-based array indices:
    /// df[i-1] = df_i
    /// df[L] = df_total
    ///

    ///
    /// among   |   within  |  df
    /// -----   |   ------  |  ---
    /// 1       |   T       |  n_groups_at_level(1) - 1
    /// 2       |   1       |  n_groups_at_level(2) - n_groups_at_level(1)
    /// ...     |   ...     |  ...
    /// L-1     |   L-2     |  n_groups_at_level(L-1) - n_groups_at_level(L-2)
    /// L       |   L-1     |  n_groups_at_level(0)(==nInd) - n_groups_at_level(L-1)


    lvlidx = 0;
    // $ df_1 = k_1 - k_0 $
    amova->df[lvlidx] = mtd->level2groupIndices[lvlidx]->len - 1;
    ++lvlidx;

    while (lvlidx < nLevels - 1) {
        // $ df_{lvlidx+1} = k_{lvlidx+1} - k_{lvlidx} $
        amova->df[lvlidx] = mtd->level2groupIndices[lvlidx]->len - mtd->level2groupIndices[lvlidx - 1]->len;
        ++lvlidx;
    }

    // $ df_L = N - k_L $
    amova->df[lvlidx] = nInd - mtd->level2groupIndices[nLevels - 2]->len;
    ++lvlidx;

    // $ df_{total} = N - 1 $
    amova->df_total = nInd - 1;

    double sum;

    size_t iti;
    size_t itj;
    double val;

    size_t idx;

    // -> get vmat
    // store UTID vmat matrix as LTID with i=itj and j=iti

    idx = 0;
    for (iti = 0;iti < nLevels;++iti) {
        for (itj = iti;itj < nLevels;++itj) {

            val = 0.0;
            do {

                //
                // v_ij special cases:
                // [case 1] i==j (iti==itj)
                //          v_ij = N
                // [case 2] j==L (itj==nLevels - 1)
                //          v_ij = G_i


                if (iti == itj) {
                    val = (double)nInd;
                    break;
                }

                size_t G_iti = mtd->level2groupIndices[iti]->len;
                if (itj == nLevels - 1) {
                    val = (double)G_iti;
                    break;
                }

                for (size_t g_iti = 0;g_iti < G_iti;++g_iti) {
                    double innersum = 0.0;

                    size_t group_g_iti = mtd->level2groupIndices[iti]->d[g_iti];
                    size_t N_g_iti = mtd->group2indIndices[group_g_iti]->len;
                    size_t nSubgroups_of_group_g_iti = mtd->group2subgroupIndices[group_g_iti]->len;
                    for (size_t g_itj = 0;g_itj < nSubgroups_of_group_g_iti;++g_itj) {
                        size_t group_g_itj = mtd->group2subgroupIndices[group_g_iti]->d[g_itj];
                        if (itj != mtd->group2levelIndices->d[group_g_itj]) {
                            continue;
                        }
                        // only use subgroups of g_iti that are from itj-th level
                        int N_g_itj = (int)mtd->group2indIndices[group_g_itj]->len;
                        innersum += SQUARE(N_g_itj);
                    }
                    innersum = innersum / (double)N_g_iti;
                    val += innersum;
                }
            } while (0);

            amova->vmat[idx] = val;
            ++idx;

        }
    }

    // -> get cmat

    idx = 0;
    for (iti = 0;iti < nLevels;++iti) {
        for (itj = iti;itj < nLevels;++itj) {

            val = 0.0;

            do {
                // c_ij cases:
                // [case 0] i > j (iti > itj) 
                //          c_ij = 0
                // [case 1] i==1 (iti==0)
                //          c_ij = (1/df_i) * (v_ij - \sum_{g_j=1}^{G_j} N_{g_j}^2 / N)
                // [case 2] 1 < i <= j
                //          c_ij = (1/df_i) * (v_ij - v_{i-1,j})

                if (iti > itj) {
                    break;
                }

                double left;
                double right;

                // left = amova->vmat[MATRIX_GET_INDEX_UTID_IJ(iti, itj, nLevels)];
                left = amova->vmat[idx];

                if (iti == 0) {

                    right = 0.0;

                    if (itj == nLevels - 1) {
                        right = 1.0;
                    } else {

                        size_t G_itj = mtd->level2groupIndices[itj]->len;
                        for (size_t g_itj = 0;g_itj < G_itj;++g_itj) {
                            size_t group_g_itj = mtd->level2groupIndices[itj]->d[g_itj];
                            int N_g_itj = mtd->group2indIndices[group_g_itj]->len;
                            right += SQUARE(N_g_itj);
                        }
                        right = right / (double)nInd;
                    }

                } else {
                    // 1 < i <= j (0 < iti <= itj)
                    right = amova->vmat[MATRIX_GET_INDEX_UTID_IJ(iti - 1, itj, nLevels)];

                }

                val = (1.0 / (double)amova->df[iti]) * (left - right);

            } while (0);

            // amova->cmat[MATRIX_GET_INDEX_UTID_IJ(iti, itj, nLevels)] = val;
            amova->cmat[idx] = val;

            ++idx;
        }
    }


    for (iti = 0;iti < nLevels;++iti) {
        // itj==iti
        const size_t idx = MATRIX_GET_INDEX_UTID_IJ(iti, iti, nLevels);
        amova->lmat[idx] = 1.0 / amova->cmat[idx];
    }


    // go reverse to avoid recursive function call
    // so the lmat values needed for the inner p loop are already calculated
    iti = nLevels;
    while (1) {
        if (iti == 0) {
            break;
        }
        --iti;

        for (size_t itj = iti + 1;itj < nLevels;++itj) {
            sum = 0.0;
            for (size_t p = iti + 1;p <= itj;++p) {
                sum += amova->cmat[MATRIX_GET_INDEX_UTID_IJ(iti, p, nLevels)] * amova->lmat[MATRIX_GET_INDEX_UTID_IJ(p, itj, nLevels)];
            }
            amova->lmat[MATRIX_GET_INDEX_UTID_IJ(iti, itj, nLevels)] = -1.0 * (sum / amova->cmat[MATRIX_GET_INDEX_UTID_IJ(iti, iti, nLevels)]);

        }
    }

}

struct simple_pthread_data_t {
    int job_id;
    void* shared_data;
    void* private_data;
};


void* amova_run_private(void* data) {
    DEVASSERT(data != NULL);

    const size_t runidx = (size_t)((simple_pthread_data_t*)data)->job_id;

    amovaStruct* amova = (amovaStruct*)((simple_pthread_data_t*)data)->shared_data;
    dmat_t* dm = (dmat_t*)((simple_pthread_data_t*)data)->private_data;
    double* matrix = dm->matrix[runidx];

    //TODO investigate the squared transform
    if (DMAT_TRANSFORM_SQUARE == dm->transform) {
        // ok
    } else if (DMAT_TRANSFORM_NONE == dm->transform) {
        for (size_t i = 0;i < dm->size;++i) {
            matrix[i] = SQUARE(matrix[i]);
        }
    }


    metadataStruct* mtd = amova->metadata;
    double* ss_total = amova->ss_total + runidx;
    double* ssd_total = amova->ssd_total + runidx;
    double* msd_total = amova->msd_total + runidx;

    size_t iti;
    size_t itj;

    double* ss = amova->ss[runidx];
    double* ssd = amova->ssd[runidx];
    double* msd = amova->msd[runidx];
    double* sigmasq = amova->sigmasq[runidx];
    double* phi_xt = amova->phi_xt[runidx];
    double* phi_xy = (amova->phi_xy == NULL ? NULL : amova->phi_xy[runidx]);

    const int nInd = mtd->nInd;
    const size_t nLevels = (size_t)mtd->nLevels;

    double sum;
    size_t groupIndex;
    size_t* pairsInGroup = NULL;
    size_t nPairsInGroup;
    int nIndsInGroup;

    size_t p;

    size_t lvlidx; // 0-based level (e.g. 1 for level 2 in Ind~Level1/Level2 and 2 for level Ind)


    /// -----------------------------------------------------------------------
    /// -> SUM OF SQUARES WITHIN
    ///

    lvlidx = 0;
    // except within ind level (lvl=L; lvlidx=L-1)
    while (lvlidx < nLevels - 1) {
        for (size_t g = 0; g < mtd->level2groupIndices[lvlidx]->len; ++g) {
            sum = 0.0;

            groupIndex = mtd->level2groupIndices[lvlidx]->d[g];
            nIndsInGroup = mtd->group2indIndices[groupIndex]->len;
            pairsInGroup = mtd->group2pairIndices[groupIndex]->d;
            nPairsInGroup = mtd->group2pairIndices[groupIndex]->len;

            for (p = 0; p < nPairsInGroup; ++p) {
                DEVASSERT(pairsInGroup[p] < dm->size);
                sum += matrix[pairsInGroup[p]];
            }

            ss[lvlidx] += sum / nIndsInGroup;

        }
        ++lvlidx;
    }

    // -> ss total
    for (size_t j = 0; j < dm->size; ++j) {
        *ss_total += matrix[j];
    }
    *ss_total = *ss_total / nInd;

    /// -----------------------------------------------------------------------
    /// -> SUM OF SQUARED DEVIATIONS (ssd)
    /// @note ssd[i] = SSD at the ith row of the AMOVA table
    ///
    /// among   |   within  |  ssd (ss_within - ss_among)
    /// -----   |   ------  |  ---
    /// 1       |   T       | ss_total[0] - ss_within_lvl_1
    /// 2       |   1       | ss_within_lvl_1 - ss_within_lvl_2
    /// ...     |   ...     | ...
    /// L-1     |   L-2     | ss_within_lvl_L-2 - ss_within_lvl_L-1
    /// L       |   L-1     | ss_within_lvl_L-1 - 0*
    /// 
    /// * (- 0 because currently within ind is disabled; when enabled, ss[L-1 (last level)] - ss[Individuals])

    lvlidx = 0;
    // ss total - ss within lvl 1
    ssd[lvlidx] = *ss_total - ss[lvlidx];
    *ssd_total = ssd[lvlidx];
    ++lvlidx;

    while (lvlidx < nLevels - 1) {
        // ss within level lvlidx-1 - ss within level lvlidx
        ssd[lvlidx] = ss[lvlidx - 1] - ss[lvlidx];
        *ssd_total += ssd[lvlidx];
        ++lvlidx;
    }
    // ss within lvl l-1 - 0*
    ssd[lvlidx] = ss[lvlidx - 1];
    *ssd_total += ssd[lvlidx];



    /// -----------------------------------------------------------------------
    /// -> MEAN SQUARED DEVIATIONS

    lvlidx = 0;
    while (lvlidx < nLevels) {
        msd[lvlidx] = ssd[lvlidx] / amova->df[lvlidx];
        ++lvlidx;
    }

    *msd_total = *ssd_total / amova->df_total;

    /// -----------------------------------------------------------------------
    /// -> SIGMA SQUARED (VARIANCE COMPONENTS)

    size_t idx;
    idx = 0;
    for (iti = 0;iti < nLevels;++iti) {
        sigmasq[iti] = 0.0;
        for (itj = iti; itj < nLevels; ++itj) {
            // sigmasq[iti] += msd[itj] * amova->lmat[MATRIX_GET_INDEX_UTID_IJ(iti, itj, nLevels)];
            sigmasq[iti] += msd[itj] * amova->lmat[idx];
            ++idx;
        }
    }
    amova->sigmasq_total[0] = 0.0;
    for (size_t i = 0; i < nLevels; ++i) {
        amova->sigmasq_total[0] += sigmasq[i];
    }

    /// -----------------------------------------------------------------------
    /// -> PHI STATISTICS

    // -> phi_xy
    // only run if nLevels > 2
    if (phi_xy != NULL) {
        for (iti = 1; iti < nLevels - 1;++iti) {
            sum = 0.0;
            for (itj = iti; itj < nLevels;++itj) {
                sum += sigmasq[itj];
            }
            phi_xy[iti - 1] = sigmasq[iti] / sum;
        }
    }
    // phi_xt
    for (iti = 0;iti < nLevels - 1;++iti) {
        sum = 0.0;
        for (itj = 0; itj <= iti;++itj) {
            sum += sigmasq[itj];
        }
        phi_xt[iti] = sum / amova->sigmasq_total[0];
    }

    return(NULL);

}


void amovaStruct_print_as_csv(amovaStruct* amova, metadataStruct* mtd, const char* bootstrap_results) {

    //  header
    //  type,label,value
    //  SSD,Among_region,0.1234
    //  fprintf(fp, "type,label,value\n");


    const size_t nLevels = (size_t)mtd->nLevels;

    outFiles->out_amova_fs->kbuf = kbuf_init();
    kstring_t* kbuf = outFiles->out_amova_fs->kbuf;

    ksprintf(kbuf, "df,Total,%d\n", amova->df_total);
    ksprintf(kbuf, "SSD,Total,%f\n", amova->ssd_total[0]);
    ksprintf(kbuf, "MSD,Total,%f\n", amova->msd_total[0]);


    // among       idx |   within    idx
    // -----   0-based |   ------    0-based 
    // i=1       0     |   T         -    
    // i=2       1     |   j=i-1=1   0
    // ...       ...   |   ...       ...
    // L-1       L-2   |   L-2       L-3
    // L         L-1   |   L-1       L-2

    size_t amonglvlidx;

    amonglvlidx = 0;
    ksprintf(kbuf, "df,Among_%s_within_Total,%d\n", mtd->levelNames->d[amonglvlidx], amova->df[amonglvlidx]);
    ksprintf(kbuf, "SSD,Among_%s_within_Total,%f\n", mtd->levelNames->d[amonglvlidx], amova->ssd[0][amonglvlidx]);
    ksprintf(kbuf, "MSD,Among_%s_within_Total,%f\n", mtd->levelNames->d[amonglvlidx], amova->msd[0][amonglvlidx]);

    ++amonglvlidx;

    while (amonglvlidx < nLevels) {

        ksprintf(kbuf, "df,Among_%s_within_%s,%d\n", mtd->levelNames->d[amonglvlidx], mtd->levelNames->d[amonglvlidx - 1], amova->df[amonglvlidx]);
        ksprintf(kbuf, "SSD,Among_%s_within_%s,%f\n", mtd->levelNames->d[amonglvlidx], mtd->levelNames->d[amonglvlidx - 1], amova->ssd[0][amonglvlidx]);
        ksprintf(kbuf, "MSD,Among_%s_within_%s,%f\n", mtd->levelNames->d[amonglvlidx], mtd->levelNames->d[amonglvlidx - 1], amova->msd[0][amonglvlidx]);

        ++amonglvlidx;
    }

    // -> phi_xy
    // only run if nLevels > 2
    if (amova->phi_xy != NULL) {
        for (size_t iti = 1; iti < nLevels - 1;++iti) {
            ksprintf(kbuf, "Phi,%s_in_%s,%f\n", mtd->levelNames->d[iti], mtd->levelNames->d[iti - 1], amova->phi_xy[0][iti - 1]);
        }
    }

    // phi_xt
    for (size_t iti = 0;iti < (size_t)(mtd->nLevels - 1);++iti) {
        ksprintf(kbuf, "Phi,%s_in_Total,%f\n", mtd->levelNames->d[iti], amova->phi_xt[0][iti]);
    }


    for (size_t iti = 0;iti < nLevels;++iti) {
        for (size_t itj = iti;itj < nLevels;++itj) {
            ksprintf(kbuf, "Variance_coefficient,c_%ld_%ld,%f\n", iti, itj, amova->cmat[MATRIX_GET_INDEX_UTID_IJ(iti, itj, nLevels)]);
        }
    }

    for (size_t i = 0;i < nLevels;++i) {
        ksprintf(kbuf, "Variance_component,%s,%f\n", mtd->levelNames->d[i], amova->sigmasq[0][i]);
        ksprintf(kbuf, "Percentage_variance,%s,%f\n", mtd->levelNames->d[i], (amova->sigmasq[0][i] / amova->sigmasq_total[0]) * 100.0);
    }

    if (NULL != bootstrap_results) {
        ksprintf(kbuf, "%s", bootstrap_results);
    }

    outFiles->out_amova_fs->kbuf_write();

}


void amovaStruct_print_as_table(amovaStruct* amova, metadataStruct* mtd) {

    kstring_t kbuf = KS_INITIALIZE;
    ksprintf(&kbuf, "=== AMOVA ======================================================================\n");
    ksprintf(&kbuf, "Formula: %s\n\n", args->formula);
    ksprintf(&kbuf, "Source of variation%-30sd.f.%-6sSSD%-10sMSD\n", " ", " ", " ");
    ksprintf(&kbuf, "--------------------------------------------------------------------------------\n");
    size_t amonglvlidx = 0;
    kstring_t tmp = KS_INITIALIZE;
    ksprintf(&tmp, "Among %s within %s", mtd->levelNames->d[amonglvlidx], "Total");
    ksprintf(&kbuf, "%-49s%-10d%-13f%-13f\n", tmp.s, amova->df[amonglvlidx], amova->ssd[0][amonglvlidx], amova->msd[0][amonglvlidx]);
    ++amonglvlidx;
    while (amonglvlidx < mtd->nLevels) {
        ks_clear(&tmp);
        ksprintf(&tmp, "Among %s within %s", mtd->levelNames->d[amonglvlidx], mtd->levelNames->d[amonglvlidx - 1]);
        ksprintf(&kbuf, "%-49s%-10d%-13f%-13f\n", tmp.s, amova->df[amonglvlidx], amova->ssd[0][amonglvlidx], amova->msd[0][amonglvlidx]);
        ++amonglvlidx;
    }
    ksprintf(&kbuf, "\n\n");


    ksprintf(&kbuf, "\nVariance coefficients:\n");
    for (size_t iti = 0;iti < mtd->nLevels;++iti) {
        for (size_t itj = iti;itj < mtd->nLevels;++itj) {
            ksprintf(&kbuf, "c_%ld_%ld\t%f\n", iti, itj, amova->cmat[MATRIX_GET_INDEX_UTID_IJ(iti, itj, mtd->nLevels)]);
        }
    }

    ksprintf(&kbuf, "\n\n");

    // print variance components
    ksprintf(&kbuf, "Variance components:\n");
    for (size_t i = 0;i < mtd->nLevels;++i) {
        ksprintf(&kbuf, "%s\t%f\t%f%%\n", mtd->levelNames->d[i], amova->sigmasq[0][i], (amova->sigmasq[0][i] / amova->sigmasq_total[0]) * 100.0);
    }

    ksprintf(&kbuf, "================================================================================\n");

    ksprintf(&kbuf, "\n\n");
    fprintf(stdout, "%s\n", kbuf.s);

    ks_free(&kbuf);
    ks_free(&tmp);

    return;
}

void amovaStruct_destroy(amovaStruct* amova) {

    amova->metadata = NULL;

    for (size_t i = 0;i < amova->nRuns;++i) {
        FREE(amova->ss[i]);
        FREE(amova->ssd[i]);
        FREE(amova->msd[i]);
        FREE(amova->sigmasq[i]);
        FREE(amova->phi_xt[i]);
        if (amova->phi_xy != NULL) {
            FREE(amova->phi_xy[i]);
        }
    }

    FREE(amova->df);
    FREE(amova->ss);
    FREE(amova->ss_total);
    FREE(amova->ssd);
    FREE(amova->ssd_total);
    FREE(amova->msd);
    FREE(amova->msd_total);
    FREE(amova->vmat);
    FREE(amova->cmat);
    FREE(amova->lmat);
    FREE(amova->sigmasq);
    FREE(amova->sigmasq_total);
    FREE(amova->phi_xt);
    if (amova->phi_xy != NULL) {
        FREE(amova->phi_xy);
    }

    FREE(amova);

}

amovaStruct* amovaStruct_init(metadataStruct* mtd, const int nAmovaRuns) {

    amovaStruct* ret = (amovaStruct*)malloc(sizeof(amovaStruct));
    ret->metadata = mtd;

    const size_t nLevels = (size_t)mtd->nLevels;

    const size_t nRuns = (size_t)nAmovaRuns;
    ret->nRuns = nRuns;


    const size_t nCmat = (nLevels * (nLevels + 1)) / 2;
    ret->vmat = NULL;
    ret->vmat = (double*)malloc((nCmat) * sizeof(double));
    ASSERT(ret->vmat != NULL);
    ret->cmat = NULL;
    ret->cmat = (double*)malloc((nCmat) * sizeof(double));
    ASSERT(ret->cmat != NULL);
    ret->lmat = NULL;
    ret->lmat = (double*)malloc((nCmat) * sizeof(double));
    ASSERT(ret->cmat != NULL);
    for (size_t i = 0; i < nCmat; ++i) {
        ret->cmat[i] = 0.0;
        ret->lmat[i] = 0.0;
        ret->vmat[i] = 0.0;
    }


    ret->df = NULL;
    ret->df = (int*)malloc((nLevels) * sizeof(int));
    ASSERT(ret->df != NULL);

    ret->df_total = 0;

    ret->ss = NULL;
    ret->ss = (double**)malloc((nRuns) * sizeof(double*));
    ASSERT(ret->ss != NULL);

    ret->ss_total = NULL;
    ret->ss_total = (double*)malloc((nRuns) * sizeof(double));
    ASSERT(ret->ss_total != NULL);

    ret->ssd = NULL;
    ret->ssd = (double**)malloc((nRuns) * sizeof(double*));
    ASSERT(ret->ssd != NULL);

    ret->ssd_total = NULL;
    ret->ssd_total = (double*)malloc((nRuns) * sizeof(double));
    ASSERT(ret->ssd_total != NULL);

    ret->msd = NULL;
    ret->msd = (double**)malloc((nRuns) * sizeof(double*));
    ASSERT(ret->msd != NULL);

    ret->msd_total = NULL;
    ret->msd_total = (double*)malloc((nRuns) * sizeof(double));
    ASSERT(ret->msd_total != NULL);

    ret->sigmasq = NULL;
    ret->sigmasq = (double**)malloc((nRuns) * sizeof(double*));
    ASSERT(ret->sigmasq != NULL);

    ret->sigmasq_total = NULL;
    ret->sigmasq_total = (double*)malloc((nRuns) * sizeof(double));
    ASSERT(ret->sigmasq_total != NULL);

    ret->phi_xt = NULL;
    ret->phi_xt = (double**)malloc((nRuns) * sizeof(double*));
    ASSERT(ret->phi_xt != NULL);


    ret->phi_xy = NULL;
    if (nLevels > 2) {
        ret->phi_xy = (double**)malloc((nRuns) * sizeof(double*));
        ASSERT(ret->phi_xy != NULL);
    }

    for (size_t i = 0;i < nRuns;++i) {
        ret->ss[i] = (double*)malloc((nLevels) * sizeof(double));
        ASSERT(ret->ss[i] != NULL);

        ret->ss_total[i] = 0.0;

        ret->ssd[i] = (double*)malloc((nLevels) * sizeof(double));
        ASSERT(ret->ssd[i] != NULL);

        ret->ssd_total[i] = 0.0;

        ret->msd[i] = (double*)malloc((nLevels) * sizeof(double));
        ASSERT(ret->msd[i] != NULL);

        ret->msd_total[i] = 0.0;

        ret->sigmasq[i] = (double*)malloc((nLevels) * sizeof(double));
        ASSERT(ret->sigmasq[i] != NULL);

        ret->sigmasq_total[i] = 0.0;

        ret->phi_xt[i] = (double*)malloc((nLevels - 1) * sizeof(double));
        ASSERT(ret->phi_xt[i] != NULL);

        if (ret->phi_xy != NULL) {
            ret->phi_xy[i] = NULL;
            ret->phi_xy[i] = (double*)malloc((nLevels - 2) * sizeof(double));
            ASSERT(ret->phi_xy[i] != NULL);
            for (size_t j = 0;j < nLevels - 2;++j) {
                ret->phi_xy[i][j] = 0.0;
            }
        }


        for (size_t j = 0;j < nLevels;++j) {
            ret->ss[i][j] = 0.0;
            ret->ssd[i][j] = 0.0;
            ret->msd[i][j] = 0.0;
            ret->sigmasq[i][j] = 0.0;
            if (j < nLevels - 1) {
                ret->phi_xt[i][j] = 0.0;
            }
        }

    }

    return(ret);
}


amovaStruct* amovaStruct_get(paramStruct* pars, metadataStruct* mtd) {


    // nRuns = n bootstrap runs + 1 (the original run)
    const int nRuns = (args->nBootstraps > 0) ? (args->nBootstraps + 1) : 1;


    amovaStruct* amova = amovaStruct_init(mtd, nRuns);

    const int maxnThreads = (args->nThreads == 0) ? 1 : args->nThreads;

    amova_run_shared(amova);

    pthread_t threads[nRuns];
    simple_pthread_data_t data[nRuns];

    int nJobsAlive = 0;

    size_t run_to_wait = 0;
    size_t runidx = 0;

    while (runidx < nRuns) {

        while (1) {
            if (nJobsAlive < maxnThreads) {
                break;
            }
            while (nJobsAlive >= maxnThreads) {
                // wait for the run that was sent first
                if (0 != pthread_join(threads[run_to_wait], NULL)) {
                    ERROR("Problem with joining the thread.");
                }
                ++run_to_wait;
                nJobsAlive--;
            }
        }

        data[runidx].job_id = runidx;
        data[runidx].shared_data = (void*)amova;
        data[runidx].private_data = (void*)pars->dm;

        if (0 != pthread_create(&threads[runidx], NULL, amova_run_private, (void*)&data[runidx])) {
            ERROR("Problem with the spawning thread.");
        }
        nJobsAlive++;

        ++runidx;
    }

    while (nJobsAlive > 0) {
        if (0 != pthread_join(threads[run_to_wait], NULL)) {
            ERROR("Problem with joining the thread.");
        }
        ++run_to_wait;
        nJobsAlive--;
    }

    kstring_t kbuf_csv = KS_INITIALIZE;
    kstring_t kbuf_table = KS_INITIALIZE;

    if (nRuns > 1) {

        double mean, sd, margin_of_error, ci_lower, ci_upper;

        double ci = args->bootstrap_ci;
        int nReps = nRuns - 1;

        ksprintf(&kbuf_csv, "Block_Bootstrapping,nReplicates,%d\n", nReps);
        ksprintf(&kbuf_csv, "Block_Bootstrapping,Confidence_Interval,%f\n", ci);
        ksprintf(&kbuf_table, "Block Bootstrapping:\n");
        ksprintf(&kbuf_table, "Number of replicates: %d\n", nReps);
        ksprintf(&kbuf_table, "Confidence interval: %f\n", ci);

        for (size_t i = 0;i < (mtd->nLevels - 1);++i) {

            mean = 0.0;
            for (size_t r=0;r < nReps;++r) {
                mean += amova->phi_xt[r][i];
            }
            mean = mean / (double)nReps;

            sd = 0.0;
            for (size_t r=0;r < nReps;++r) {
                sd += pow(amova->phi_xt[r][i] - mean, 2);
            }
            sd = sqrt(sd / (nReps-1)); // sample stdev 

            margin_of_error= 1.96 * (sd / sqrt((double)nReps));
            ci_lower = mean - margin_of_error;
            ci_upper = mean + margin_of_error;

            ksprintf(&kbuf_csv, "Block_Bootstrapping_Phi_Mean,%s_in_Total,%f\n", mtd->levelNames->d[i], mean);
            ksprintf(&kbuf_csv, "Block_Bootstrapping_Phi_SD,%s_in_Total,%f\n", mtd->levelNames->d[i], sd);
            ksprintf(&kbuf_csv, "Block_Bootstrapping_Phi_LowerCI,%s_in_Total,%f\n", mtd->levelNames->d[i], ci_lower);
            ksprintf(&kbuf_csv, "Block_Bootstrapping_Phi_UpperCI,%s_in_Total,%f\n", mtd->levelNames->d[i], ci_upper);

            ksprintf(&kbuf_table, "Phi(%s in Total): %f\n", mtd->levelNames->d[i], mean);
            ksprintf(&kbuf_table, "Phi(%s in Total) SD: %f\n", mtd->levelNames->d[i], sd);
            ksprintf(&kbuf_table, "Phi(%s in Total) Lower CI: %f\n", mtd->levelNames->d[i], ci_lower);
            ksprintf(&kbuf_table, "Phi(%s in Total) Upper CI: %f\n", mtd->levelNames->d[i], ci_upper);

        }

        if (amova->phi_xy != NULL) {
            for (size_t i = 0;i < (mtd->nLevels - 2);++i) {

                mean = 0.0;
                for (size_t r=0;r < nReps;++r) {
                    mean += amova->phi_xy[r][i];
                }
                mean = mean / (double)nReps;

                sd = 0.0;
                for (size_t r=0;r < nReps;++r) {
                    sd += pow(amova->phi_xy[r][i] - mean, 2);
                }
                sd = sqrt(sd / (nReps-1)); // sample stdev

                margin_of_error= 1.96 * (sd / sqrt((double)nReps));
                ci_lower = mean - margin_of_error;
                ci_upper = mean + margin_of_error;

                ksprintf(&kbuf_csv, "Block_Bootstrapping_Phi_Mean,%s_in_%s,%f\n", mtd->levelNames->d[i + 1], mtd->levelNames->d[i], mean);
                ksprintf(&kbuf_csv, "Block_Bootstrapping_Phi_SD,%s_in_%s,%f\n", mtd->levelNames->d[i + 1], mtd->levelNames->d[i], sd);
                ksprintf(&kbuf_csv, "Block_Bootstrapping_Phi_LowerCI,%s_in_%s,%f\n", mtd->levelNames->d[i + 1], mtd->levelNames->d[i], ci_lower);
                ksprintf(&kbuf_csv, "Block_Bootstrapping_Phi_UpperCI,%s_in_%s,%f\n", mtd->levelNames->d[i + 1], mtd->levelNames->d[i], ci_upper);

                ksprintf(&kbuf_table, "Phi(%s in %s): %f\n", mtd->levelNames->d[i + 1], mtd->levelNames->d[i], mean);
                ksprintf(&kbuf_table, "Phi(%s in %s) SD: %f\n", mtd->levelNames->d[i + 1], mtd->levelNames->d[i], sd);
                ksprintf(&kbuf_table, "Phi(%s in %s) Lower CI: %f\n", mtd->levelNames->d[i + 1], mtd->levelNames->d[i], ci_lower);
                ksprintf(&kbuf_table, "Phi(%s in %s) Upper CI: %f\n", mtd->levelNames->d[i + 1], mtd->levelNames->d[i], ci_upper);

            }
        }
    }

    amovaStruct_print_as_csv(amova, mtd, kbuf_csv.s);
    if (kbuf_csv.l > 0) {
        ks_free(&kbuf_csv);
    }

    if (args->printAmovaTable == 1) {
        amovaStruct_print_as_table(amova, mtd);
    }

    return (amova);
}

