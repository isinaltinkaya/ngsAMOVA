#include "amova.h"


/// @brief calculate part of AMOVA statistics where the calculations are independent of the distance matrix but depend on the metadata
/// @param amv amovaStruct
/// @param mtd metadataStruct
/// @return void
/// @details
/// perform shared AMOVA calculations depends on the metadata but does not depend on the distance matrix
/// therefore this function can be called once for each metadataStruct
/// i.e. one call is enough for all amova replicates
/// 
///     
inline void amova_run_shared(amovaStruct* amv) {

    metadataStruct* mtd = amv->metadata;
    const int nInd = mtd->nInd;

    const size_t nLevels = (size_t)mtd->nLevels; // L
    size_t lvlidx; // 0-based lvl

    /// ------------------------------------------------------------------- ///
    /// -> DEGREES OF FREEDOM (amv->df)
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
    amv->df[lvlidx] = mtd->level2groupIndices[lvlidx]->len - 1;
    ++lvlidx;

    while (lvlidx < nLevels - 1) {
        // $ df_{lvlidx+1} = k_{lvlidx+1} - k_{lvlidx} $
        amv->df[lvlidx] = mtd->level2groupIndices[lvlidx]->len - mtd->level2groupIndices[lvlidx - 1]->len;
        ++lvlidx;
    }

    // $ df_L = N - k_L $
    amv->df[lvlidx] = nInd - mtd->level2groupIndices[nLevels - 2]->len;
    ++lvlidx;

    // $ df_{total} = N - 1 $
    amv->df_total = nInd - 1;

    double sum;

    size_t iti;
    size_t itj;
    double val;


    // -> get vmat

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

            amv->vmat[GET_UPTRID_MATRIX_IJ(iti, itj)] = val;

        }
    }

    // -> get cmat

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

                left = amv->vmat[GET_UPTRID_MATRIX_IJ(iti, itj)];

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
                    right = amv->vmat[GET_UPTRID_MATRIX_IJ(iti - 1, itj)];

                }

                val = (1.0 / (double)amv->df[iti]) * (left - right);

            } while (0);

            amv->cmat[GET_UPTRID_MATRIX_IJ(iti, itj)] = val;

        }
    }


    for (iti = 0;iti < nLevels;++iti) {
        // itj==iti
        const size_t idx = GET_UPTRID_MATRIX_IJ(iti, iti);
        amv->lmat[idx] = 1.0 / amv->cmat[idx];
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
                sum += amv->cmat[GET_UPTRID_MATRIX_IJ(iti, p)] * amv->lmat[GET_UPTRID_MATRIX_IJ(p, itj)];
            }
            amv->lmat[GET_UPTRID_MATRIX_IJ(iti, itj)] = -1.0 * (sum / amv->cmat[GET_UPTRID_MATRIX_IJ(iti, iti)]);

        }
    }

}

struct simple_pthread_data_t {
    int job_id;
    void* shared_data;
    void* private_data;
};


void* amova_run_private(void* data) {

    amovaStruct* amv = (amovaStruct*)((simple_pthread_data_t*)data)->shared_data;
    distanceMatrixStruct* dm = (distanceMatrixStruct*)((simple_pthread_data_t*)data)->private_data;
    const size_t runidx = (size_t)((simple_pthread_data_t*)data)->job_id;
    metadataStruct* mtd = amv->metadata;
    double* ss_total = amv->ss_total + runidx;
    double* ssd_total = amv->ssd_total + runidx;
    double* msd_total = amv->msd_total + runidx;

    size_t iti;
    size_t itj;

    double* ss = amv->ss[runidx];
    double* ssd = amv->ssd[runidx];
    double* msd = amv->msd[runidx];
    double* sigmasq = amv->sigmasq[runidx];
    double* phi_xt = amv->phi_xt[runidx];
    double* phi_xy = (amv->phi_xy == NULL ? NULL : amv->phi_xy[runidx]);

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
                sum += dm->M[pairsInGroup[p]];
            }

            ss[lvlidx] += sum / nIndsInGroup;

        }
        ++lvlidx;
    }

    // -> ss total
    for (int j = 0; j < dm->nIndCmb; ++j) {
        *ss_total += dm->M[j];
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
        msd[lvlidx] = ssd[lvlidx] / amv->df[lvlidx];
        ++lvlidx;
    }

    // msd total = ssd total / df total
    *msd_total = *ssd_total / amv->df_total;

    /// -----------------------------------------------------------------------
    /// -> SIGMA SQUARED (VARIANCE COMPONENTS)

    for (iti = 0;iti < nLevels;++iti) {
        sigmasq[iti] = 0.0;
        for (itj = iti; itj < nLevels; ++itj) {
            sigmasq[iti] += msd[itj] * amv->lmat[GET_UPTRID_MATRIX_IJ(iti, itj)];
        }
    }
    amv->sigmasq_total[0] = 0.0;
    for (size_t i = 0; i < nLevels; ++i) {
        amv->sigmasq_total[0] += sigmasq[i];
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
        phi_xt[iti] = sum / amv->sigmasq_total[0];
    }


    return(NULL);

}


void amovaStruct_print_as_csv(amovaStruct* amv, metadataStruct* mtd) {

    //  header
    //  type,label,value
    //  SSD,Among_region,0.1234
    //  fprintf(fp, "type,label,value\n");


    const size_t nLevels = (size_t)mtd->nLevels;

    outFiles->out_amova_fs->kbuf = kbuf_init();
    kstring_t* kbuf = outFiles->out_amova_fs->kbuf;

    ksprintf(kbuf, "df,Total,%d\n", amv->df_total);
    ksprintf(kbuf, "SSD,Total,%f\n", amv->ssd_total[0]);
    ksprintf(kbuf, "MSD,Total,%f\n", amv->msd_total[0]);


    // among       idx |   within    idx
    // -----   0-based |   ------    0-based 
    // i=1       0     |   T         -    
    // i=2       1     |   j=i-1=1   0
    // ...       ...   |   ...       ...
    // L-1       L-2   |   L-2       L-3
    // L         L-1   |   L-1       L-2

    size_t amonglvlidx;

    amonglvlidx = 0;
    ksprintf(kbuf, "df,Among_%s_within_Total,%d\n", mtd->levelNames->d[amonglvlidx], amv->df[amonglvlidx]);
    ksprintf(kbuf, "SSD,Among_%s_within_Total,%f\n", mtd->levelNames->d[amonglvlidx], amv->ssd[0][amonglvlidx]);
    ksprintf(kbuf, "MSD,Among_%s_within_Total,%f\n", mtd->levelNames->d[amonglvlidx], amv->msd[0][amonglvlidx]);

    ++amonglvlidx;

    while (amonglvlidx < nLevels) {

        ksprintf(kbuf, "df,Among_%s_within_%s,%d\n", mtd->levelNames->d[amonglvlidx], mtd->levelNames->d[amonglvlidx - 1], amv->df[amonglvlidx]);
        ksprintf(kbuf, "SSD,Among_%s_within_%s,%f\n", mtd->levelNames->d[amonglvlidx], mtd->levelNames->d[amonglvlidx - 1], amv->ssd[0][amonglvlidx]);
        ksprintf(kbuf, "MSD,Among_%s_within_%s,%f\n", mtd->levelNames->d[amonglvlidx], mtd->levelNames->d[amonglvlidx - 1], amv->msd[0][amonglvlidx]);

        ++amonglvlidx;
    }

    // -> phi_xy
    // only run if nLevels > 2
    if (amv->phi_xy != NULL) {
        for (size_t iti = 1; iti < nLevels - 1;++iti) {
            ksprintf(kbuf, "Phi,%s_in_%s,%f\n", mtd->levelNames->d[iti], mtd->levelNames->d[iti - 1], amv->phi_xy[0][iti - 1]);
        }
    }

    // phi_xt
    for (size_t iti = 0;iti < mtd->nLevels - 1;++iti) {
        ksprintf(kbuf, "Phi,%s_in_Total,%f\n", mtd->levelNames->d[iti], amv->phi_xt[0][iti]);
    }


    for (size_t iti = 0;iti < nLevels;++iti) {
        for (size_t itj = iti;itj < nLevels;++itj) {
            ksprintf(kbuf, "Variance_coefficient,c_%ld_%ld,%f\n", iti, itj, amv->cmat[GET_UPTRID_MATRIX_IJ(iti, itj)]);
        }
    }

    for (size_t i = 0;i < nLevels;++i) {
        ksprintf(kbuf, "Variance_component,%s,%f\n", mtd->levelNames->d[i], amv->sigmasq[0][i]);
        ksprintf(kbuf, "Percentage_variance,%s,%f\n", mtd->levelNames->d[i], (amv->sigmasq[0][i] / amv->sigmasq_total[0]) * 100.0);
    }


    outFiles->out_amova_fs->kbuf_write();

}


void amovaStruct_print_as_table(FILE* fp, amovaStruct* amv, metadataStruct* mtd) {
    //TODO UPDATEME

//     const size_t nLevels = (size_t)mtd->nLevels;

//     fprintf(fp, "\n");
//     fprintf(fp, "\n\n");
//     fprintf(fp, "==========================================  AMOVA  ==========================================");

//     fprintf(fp, "\n\n");
//     fprintf(fp, "Formula: %s ", args->formula);
//     fprintf(fp, "\n");
//     fprintf(fp, "Source of variation\t\t\t\t\td.f.\tSSD\t\tMSD");
//     fprintf(fp, "\n");
//     fprintf(fp, "---------------------------------------------------------------------------------------------");
//     fprintf(fp, "\n\n");
//     fprintf(fp, "\n");


//     size_t x = 1;
//     fprintf(fp, "Among %-15s\t\t\t\t\t%d\t%f\t%f", mtd->levelNames->get(x), amv->df[x], amv->ssd[x], amv->msd[x]);

//     while (x < mtd->nLevels) {
//         fprintf(fp, "\n");
//         fprintf(fp, "Among %s within %-25s\t%d\t%f\t%f", mtd->levelNames->get(x + 1), mtd->levelNames->get(x), amv->df[x], amv->ssd[x], amv->msd[x]);
//         fprintf(fp, "\n");
//         ++x;
//     }

//     fprintf(fp, "\n");
//     fprintf(fp, "Among %s within %-25s\t%d\t%f\t%f", mtd->levelNames->get(0), mtd->levelNames->get(x), amv->df[x], amv->ssd[x], amv->msd[x]);
//     fprintf(fp, "\n");


//     ++x;
//     fprintf(fp, "\n");
//     fprintf(fp, "Total\t\t\t\t\t\t\t%d\t%f\t%f", amv->df[x], amv->ssd[x], amv->msd[x]);

//     fprintf(fp, "\n\n\n");
//     fprintf(fp, "Variance components:\n\n\tVariance Component\tSigma^2\t\t%% of total");

//     for (size_t i = 0;i < nLevels;++i) {
//         fprintf(fp, "\n\t%-20s", mtd->levelNames->d[i + 1]);
//         fprintf(fp, "\t%f", amv->sigmasq[0][i]);
//         // fprintf(fp, "\t%f", pct_sigmasq[i]);
//     }

//     // Lowest level (i.e. Individual)
//     fprintf(fp, "\n\t%-20s", mtd->levelNames->d[0]);
//     fprintf(fp, "\t%f", amv->sigmasq[0][nLevels - 1]);
//     // fprintf(fp, "\t%f", pct_sigmasq[_ncoef - 1]);

//     fprintf(fp, "\n\n\n");
//     fprintf(fp, "\nVariance coefficients:\n\n\t");
//     if (mtd->nLevels == 2) {
//         fprintf(fp, "%f", amv->ncoef[0]);
//     } else {
//         for (size_t i = 0;i < nLevels;++i) {
//             fprintf(fp, "%f\t", amv->ncoef[i]);
//         }
//     }

//     fprintf(fp, "\n\n\n");
//     fprintf(fp, "Phi-statistic:\n\n");
//     for (size_t i = 0; i < ((2 * nLevels) - 3); ++i) {
//         fprintf(fp, "\t%f", amv->phi[i]);
//     }
//     fprintf(fp, "\n\n");
//     fprintf(fp, "=============================================================================================");
//     fprintf(fp, "\n\n");
}


amovaStruct* amovaStruct_get(distanceMatrixStruct* dm, metadataStruct* mtd, blobStruct* blobs) {

    // nRuns = n bootstrap runs + 1 (the original run)
    const size_t nRuns = (size_t)args->nBootstraps + 1;

    amovaStruct* amova = amovaStruct_init(mtd, nRuns);

    const size_t maxnThreads = (args->nThreads == 0) ? 1 : args->nThreads;

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
        data[runidx].private_data = (void*)dm;

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

    //     //TODO HERE BURADAKALDIN!
        // const int nReplicates = blob->bootstraps->nReplicates;
            // THREADS[i] = new amovaBootstrapThreads(blob->bootstraps->replicates[i]->amova, blob->bootstraps->replicates[i]->distanceMatrix, mtd);

        // //TODO change these with new vals
        // blob->bootstraps->nPhiValues = (mtd->nLevels * 2) - 3;
        // blob->bootstraps->phiValues = (double**)malloc(blob->bootstraps->nPhiValues * sizeof(double*));

        // for (int i = 0; i < blob->bootstraps->nPhiValues; i++) {
        //     blob->bootstraps->phiValues[i] = (double*)malloc(nReplicates * sizeof(double));

        //     for (int r = 0; r < nReplicates; r++) {
        //         // blob->bootstraps->phiValues[i][r] = THREADS[r]->amova->phi[i];
        //     }
        // }

    amovaStruct_print_as_csv(amova, mtd);

    // if (0 < args->nBootstraps) {
    //     spawnThreads_amovaBootstrap(metadata, blobs);
    //     blobs->bootstraps->print_confidenceInterval(stderr);
    // }

// if (0 < args->nBootstraps) {
//     ASSERT(blobs != NULL);


//     spawnThreads_amovaBootstrap(mtd, blobs);
//     blobs->bootstraps->print_confidenceInterval(stderr);
// }


// if (0 == doAmova(amova, dm, metadata)) {
//     fprintf(stderr, "\n[INFO]\t-> Finished running AMOVA");
//     if (0 < args->nBootstraps) {
//         fprintf(stderr, " for the original dataset.\n");
//     } else {
//         fprintf(stderr, ".\n");
//     }

// if (args->printAmovaTable == 1) {
    // amova->print_as_table(stdout);
    //TODO make this to print file instead
// }

    return (amova);
}

