#include "amova.h"

/// @brief thread handler for doAmova
void *t_doAmova(void *p) {
    amovaBootstrapThreads *THREAD = (amovaBootstrapThreads *)p;
    if (doAmova(THREAD) != 0) {
        NEVER;
    }
    return (0);
}

int doAmova(amovaBootstrapThreads *THREAD) {
    return (doAmova(THREAD->amova, THREAD->distanceMatrix, THREAD->metadata));
}

void spawnThreads_amovaBootstrap(metadataStruct *metadata, blobStruct *blob) {
    const int nReplicates = blob->bootstraps->nReplicates;

    pthread_t bootstrapThreads[nReplicates];
    amovaBootstrapThreads **THREADS = new amovaBootstrapThreads *[nReplicates];

    int nJobs_sent = 0;

    for (int i = 0; i < nReplicates; ++i) {
        blob->bootstraps->replicates[i]->amova = new amovaStruct(metadata);
        THREADS[i] = new amovaBootstrapThreads(blob->bootstraps->replicates[i]->amova, blob->bootstraps->replicates[i]->distanceMatrix, metadata);
        if (nJobs_sent == args->nThreads) {
            int t = 0;
            while (nJobs_sent > 0) {
                t = i - nJobs_sent;
                if (pthread_join(bootstrapThreads[t], NULL) != 0) {
                    ERROR("Problem with joining thread.");
                } else {
                    nJobs_sent--;
                }
            }
        }
        if (pthread_create(&bootstrapThreads[i], NULL, t_doAmova, THREADS[i]) == 0) {
            nJobs_sent++;
        } else {
            ERROR("Problem with spawning thread.");
        }
    }

    int t = 0;
    while (nJobs_sent > 0) {
        t = nReplicates - nJobs_sent;
        if (pthread_join(bootstrapThreads[t], NULL) != 0) {
            ERROR("Problem with joining thread.");
        } else {
            nJobs_sent--;
        }
    }

    blob->bootstraps->nPhiValues = THREADS[0]->amova->_phi;
    blob->bootstraps->phiValues = (double **)malloc(blob->bootstraps->nPhiValues * sizeof(double *));

    for (int i = 0; i < blob->bootstraps->nPhiValues; i++) {
        blob->bootstraps->phiValues[i] = (double *)malloc(nReplicates * sizeof(double));

        for (int r = 0; r < nReplicates; r++) {
            blob->bootstraps->phiValues[i][r] = THREADS[r]->amova->phi[i];
        }
    }

    for (int r = 0; r < nReplicates; r++) {
        delete THREADS[r];
    }
    delete[] THREADS;
}

amovaStruct *amovaStruct_get(distanceMatrixStruct *dm, metadataStruct *metadata) {
    amovaStruct *amova = new amovaStruct(metadata);

    ASSERT(NULL != amova);

    if (0 == doAmova(amova, dm, metadata)) {
        fprintf(stderr, "\n[INFO]\t-> Finished running AMOVA");
        if (0 < args->nBootstraps) {
            fprintf(stderr, " for the original dataset.\n");
        } else {
            fprintf(stderr, ".\n");
        }

        // eval_amovaStruct(amova);
        // return 0 on successful run
        // return 1 on error (either during analyses or during evaluation)
        if (args->printAmovaTable == 1) {
            amova->print_as_table(stdout, metadata);
        }
        amova->print_as_csv(metadata);
    }
    return (amova);
}

// calculate the degrees of freedom for each level
void get_degreesOfFreedom(amovaStruct *amova, metadataStruct *metadata) {
    const int nInd = metadata->nInd;

    // highest level df is nStrata - 1
    // e.g. Individual ~ Region / Population / Subpopulation
    // if there are 3 regions, then nStrata is number of unique regions - 1 = 3-1 = 2
    amova->df[0] = metadata->nGroups[0] - 1;

    // total df
    amova->df[amova->nAmovaLevels - 1] = nInd - 1;

    if (metadata->nLevels == 1) {
        // df for lowest amova level
        // e.g. Individual ~ Region / Population / Subpopulation
        // nInd - sum(for each Population, nStrata (num_Subpopulation) in Population)
        amova->df[1] = nInd - metadata->nGroups[0];
    } else if (metadata->nLevels == 2) {
        int sum0 = 0;
        for (int i = 0; i < metadata->nGroups[0]; i++) {
            sum0 += metadata->countNSubgroupAtLevel(0, i, 1);
        }
        int sum1 = metadata->nGroups[0];
        amova->df[1] = sum0 - sum1;
        amova->df[2] = nInd - sum0;
    } else {
        fprintf(stderr, "[ERROR] NOT IMPLEMENTED YET");
        ASSERT(0 == 1);
    }
}

void get_sumOfSquares(amovaStruct *amova, metadataStruct *metadata, distanceMatrixStruct *dm) {
    // calculate the sum of squares within for each level
    //
    for (size_t i = 0; i < amova->_ncoef; i++) {
        // within level sum of squares
        amova->ss[i] = calculate_SumOfSquares_Within(i, amova, dm, metadata);
    }
    // ----------------------------------------------- //
    // SUM OF SQUARED DEVIATIONS
    // ----------------------------------------------- //

    if (metadata->nLevels == 1) {
        amova->ssd[0] = amova->ss[1] - amova->ss[0];
        amova->ssd[1] = amova->ss[0];
        amova->ssd[2] = amova->ss[1];
    } else if (metadata->nLevels == 2) {
        amova->ssd[3] = amova->ss[2];
        amova->ssd[1] = amova->ss[0] - amova->ss[1];
        amova->ssd[2] = amova->ss[1];
        amova->ssd[0] = amova->ssd[3] - (amova->ssd[0] + amova->ssd[1] + amova->ssd[2]);
    } else {
        fprintf(stderr, "[ERROR] NOT IMPLEMENTED YET");
        ASSERT(0 == 1);
    }
    // get msd from ssd and df
    for (size_t i = 0; i < (size_t)amova->nAmovaLevels; i++) {
        amova->msd[i] = amova->ssd[i] / amova->df[i];
    }
}

void get_varianceCoefficients(amovaStruct *amova, metadataStruct *metadata, distanceMatrixStruct *dm) {
    if (metadata->nLevels == 1) {
        // calculate variance coefficients (n) for 1 level AMOVA

        double n_gi = 0.0;

        // n = [ N - \sum_{g \in G} ( N^2_{g}/N) ) ]  /   G - 1
        for (int sti = 0; sti < metadata->nGroups[0]; sti++) {
            ASSERT(metadata->nIndPerStrata[0] != 0);
            n_gi += SQUARE(metadata->nIndPerStrata[0][sti]);
        }
        n_gi = n_gi / (double)dm->nInd;
        amova->ncoef[0] = (double)((double)dm->nInd - (double)n_gi) / (double)(metadata->nGroups[0] - 1);

    } else if (metadata->nLevels == 2) {
        // calculate variance coefficients for 2 level AMOVA
        //
        // eq.9a-c in Excoffier1992
        // estimate the method of moments: n, n', and n''
        // ncoef[2] = {n, n', n''} == {n, n1, n2}

        // n =
        // 		(NUL) \frac{ \sum^G_{g=1} \sum^{I_g}_{i=1} N_{ig} -
        // 		(NUR)		\sum^G_{g=1} ( \frac{ \sum^{I_g}_{i=1} N^2_{ig}}{\sum^{I_g}_{i=1} N_{ig}} ) }
        // 		(NL)	{ \sum^G_{g=1} (I_g - 1) }
        //
        // n = (NUL-NUR) / NL
        // 		NUL = n formula upper left
        // 		NUR = n formula upper right
        // 		NL = n formula lower
        //
        // NUL=0;
        // for each region
        // 		for each population in region
        // 			NUL += number of individuals in population (nIndPerStrata)
        //
        // loop through the groups in the highest hierarchy level hierArr[0] (e.g. regions)
        double NUL = 0.0;

        for (int g = 0; g < metadata->nGroups[0]; g++) {
            for (int sg = 0; sg < metadata->nGroups[1]; sg++) {
                // if subgroup is a child of group; count number of individuals in subgroup
                if (metadata->groupFromParentGroup(0, g, 1, sg) == 1) {
                    NUL += metadata->countIndsInGroup(1, sg);
                }
            }
        }

        // TODO join these two up below
        //
        //  NUR=0;
        //  for each region
        //  		Nig2_g=0;
        //		Nig_g=0;
        //  		for each population in region
        //  			Nig2_g += nIndPerStrata^2
        //  		for each population in region
        //  			Nig_g += nIndPerStrata
        //  		NUR += (Nig2_g / Nig_g)
        //
        double NUR = 0.0;
        double Nig2_g = 0.0;
        double Nig_g = 0.0;

        for (int g = 0; g < metadata->nGroups[0]; g++) {
            Nig2_g = 0.0;
            Nig_g = 0.0;
            for (int sg = 0; sg < metadata->nGroups[1]; sg++) {
                // if subgroup is a child of group; count number of individuals in subgroup
                if (metadata->groupFromParentGroup(0, g, 1, sg) == 1) {
                    Nig2_g += SQUARE(metadata->countIndsInGroup(1, sg));
                    Nig_g += metadata->countIndsInGroup(1, sg);
                }
            }
            NUR += Nig2_g / Nig_g;
        }

        //
        // NL=0;
        // for each region
        // 		NL += (number of populations in region - 1 )
        //
        //
        double NL = 0.0;
        for (int g = 0; g < metadata->nGroups[0]; g++) {
            NL += metadata->countNSubgroupAtLevel(0, g, 1) - 1;
        }

        amova->ncoef[0] = (NUL - NUR) / NL;

        // n' =
        //  (NUR) \frac{ \sum^G_{g=1} (\frac{\sum^{I_g}_{j=1} N^2_{ig} }{\sum^{I_g}_{i=1} N{ig}  } ) -
        //  (N1UR) 	\frac{\sum^G_{g=1} \sum^{I_g}_{j=1} N^2_{ig} }
        //  (N1L)	{ G - 1 }
        //
        // n' = (NUR-N1UR)/N1L
        //
        // NUR => already calculated above
        //
        // N1UR =0;
        // Nig2=0;
        // Nig=0;
        // for each region
        // 		for each population in region
        // 			Nig2 += nIndPerStrata^2
        //
        // for each region
        // 		for each population in region
        // 			Nig += nIndPerStrata
        //
        // N1UR = Nig2 / Nig
        // 				Nig => already calculated above as NUL
        //
        double N1UR = 0.0;
        double Nig2 = 0.0;

        // g=group, sg=subgroup
        for (int g = 0; g < metadata->nGroups[0]; g++) {
            for (int sg = 0; sg < metadata->nGroups[1]; sg++) {
                if (metadata->groupFromParentGroup(0, g, 1, sg) == 1) {
                    Nig2 += SQUARE(metadata->nIndPerStrata[1][sg]);
                }
            }
        }

        N1UR = Nig2 / NUL;

        // N1L=0;
        // N1L = nRegions - 1;
        //
        double N1L = 0.0;
        N1L = metadata->nGroups[0] - 1;

        amova->ncoef[1] = (NUR - N1UR) / N1L;

        // n'' =
        // 	(NUL) frac{\sum^G_{g=1} \sum^{I_g}_{i=1} N_{ig} -
        //	(N2UR) 	\frac{\sum^G_{g=1} \left( \sum^{I_g}_{i=1} N_{ig} \right)^2 }
        // 	(N1L) {\sum^G_{g=1} \sum^{I_g}_{i=1} N_{ig} }}{ G - 1 }
        //
        // n'' = (NUL-N2UR)/N1L
        //
        // NUL, N1L => already calculated above
        //
        // N2UR=0.0;
        // Nig=Nig; => already calculated above
        // Ni2g=0.0;
        // for each region
        //		Ni2g=0.0
        // 		for each population in region
        // 			Ni2g += nIndPerStrata
        //		Ni2g = Ni2g ^ 2
        //
        // N2UR = Ni2g/Nig;
        //
        double N2UR = 0.0;
        double Ni2g = 0.0;
        for (int g = 0; g < metadata->nGroups[0]; g++) {
            double xx = 0.0;
            for (int sg = 0; sg < metadata->nGroups[1]; sg++) {
                if (metadata->groupFromParentGroup(0, g, 1, sg) == 1) {
                    xx += metadata->countIndsInGroup(1, sg);
                }
            }
            Ni2g += SQUARE(xx);
        }
        N2UR = Ni2g / NUL;

        amova->ncoef[2] = (NUL - N2UR) / N1L;

    } else {
        NEVER;
    }
}

// calculate variance components (sigma squared)
void get_varianceComponents(amovaStruct *amova, metadataStruct *metadata) {
    if (metadata->nLevels == 1) {
        amova->sigmasq[0] = (amova->msd[0] - amova->msd[1]) / amova->ncoef[0];
        amova->sigmasq[1] = amova->msd[1];
        amova->sigmasq_total = amova->sigmasq[0] + amova->sigmasq[1];

        // TODO add more warnings for all levels and explain
        // if(amova->sigmasq[0] <= 0){
        //     WARNING("Sigma squared 0 is %f.");
        // }

    } else if (metadata->nLevels == 2) {
        amova->sigmasq[2] = amova->msd[2];
        amova->sigmasq[1] = (amova->msd[1] - amova->sigmasq[2]) / amova->ncoef[0];
        amova->sigmasq[0] = (amova->msd[0] - amova->msd[2] - (amova->ncoef[1] * amova->sigmasq[1])) / amova->ncoef[2];
        amova->sigmasq_total = amova->sigmasq[0] + amova->sigmasq[1] + amova->sigmasq[2];
    } else {
        NEVER;
    }

    calculate_PercentageTotalVariance(amova);
}

void calculate_PercentageTotalVariance(amovaStruct *amova) {
    for (int i = 0; i < (int)amova->_ncoef; i++) {
        amova->pct_sigmasq[i] = (amova->sigmasq[i] / amova->sigmasq_total) * 100.0;
    }
}

// TODO rm metadata and use amova->nLevels instead
void get_phiStatistics(amovaStruct *amova, metadataStruct *metadata) {
    if (metadata->nLevels == 1) {
        amova->phi[0] = amova->sigmasq[0] / (amova->sigmasq[0] + amova->sigmasq[1]);

    } else if (metadata->nLevels == 2) {
        // Individual,Region,Population,Total

        // lvl0 (i.e. reg) in TOTAL
        // Phi_CT
        if (amova->sigmasq_total <= 0) {
            fprintf(stderr, "\n[ERROR]\tTotal variance is %f (sigmasq_total is <= 0). Please check your analysis settings and make sure you have enough data.\n", amova->sigmasq_total);
            exit(1);
        }

        amova->phi[0] = amova->sigmasq[0] / (amova->sigmasq_total);

        // lvl1 in lvl0
        // Phi_SC
        amova->phi[1] = amova->sigmasq[1] / (amova->sigmasq[1] + amova->sigmasq[2]);

        // lvl1 (i.e. pop) in TOTAL
        // Phi_ST
        amova->phi[2] = (amova->sigmasq[0] + amova->sigmasq[1]) / (amova->sigmasq_total);
    }
}

int doAmova(amovaStruct *amova, distanceMatrixStruct *dm, metadataStruct *metadata) {
    get_degreesOfFreedom(amova, metadata);
    get_sumOfSquares(amova, metadata, dm);
    get_varianceCoefficients(amova, metadata, dm);
    get_varianceComponents(amova, metadata);
    get_phiStatistics(amova, metadata);
    return (0);
}

void amovaStruct::_print(FILE *fp) {
    fprintf(fp, "\nnAmovaLevels = %d\n", nAmovaLevels);
    fprintf(fp, "nLevels = %d\n", nLevels);
    fprintf(fp, "\n_ncoef = %zu\n", _ncoef);
    fprintf(fp, "_phi = %zu\n", _phi);
    for (size_t i = 0; i < (size_t)nAmovaLevels; i++) {
        fprintf(fp, "ssd[%zu] = %f\n", i, ssd[i]);
        fprintf(fp, "df[%zu] = %d\n", i, df[i]);
        fprintf(fp, "msd[%zu] = %f\n", i, msd[i]);
    }
    for (size_t i = 0; i < _ncoef; i++) {
        fprintf(fp, "ss[%zu] = %f\n", i, ss[i]);
        fprintf(fp, "ncoef[%zu] = %f\n", i, ncoef[i]);
        fprintf(fp, "sigmasq[%zu] = %f\n", i, sigmasq[i]);
        fprintf(fp, "pct_sigmasq[%zu] = %f\n", i, pct_sigmasq[i]);
    }
    for (size_t i = 0; i < _phi; i++) {
        fprintf(fp, "phi[%zu] = %f\n", i, phi[i]);
    }
}

void amovaStruct::print_as_table(FILE *fp, metadataStruct *metadata) {
    fprintf(fp, "\n");
    fprintf(fp, "\n\n");
    fprintf(fp, "==========================================  AMOVA  ==========================================");
    fprintf(fp, "\n");
    fprintf(fp, "Source of variation\t\t\t\t\td.f.\tSSD\t\tMSD");
    fprintf(fp, "\n");
    fprintf(fp, "---------------------------------------------------------------------------------------------");
    fprintf(fp, "\n\n");
    fprintf(fp, "\n");

    // TODO print formula
    int x = 1;
    fprintf(fp, "Among %-15s\t\t\t\t\t%d\t%f\t%f", metadata->levelNames[x], df[0], ssd[0], msd[0]);

    while (x < metadata->nLevels + 1) {
        if (x == metadata->nLevels) {
            fprintf(fp, "\n");
            fprintf(fp, "Among %s within %-25s\t%d\t%f\t%f", metadata->levelNames[0], metadata->levelNames[metadata->nLevels], df[x], ssd[x], msd[x]);
            fprintf(fp, "\n");
        } else {
            fprintf(fp, "\n");
            fprintf(fp, "Among %s within %-25s\t%d\t%f\t%f", metadata->levelNames[x + 1], metadata->levelNames[x], df[x], ssd[x], msd[x]);
            fprintf(fp, "\n");
        }
        x++;
    }

    fprintf(fp, "\n");
    fprintf(fp, "Total\t\t\t\t\t\t\t%d\t%f\t%f", df[nAmovaLevels - 1], ssd[nAmovaLevels - 1], msd[nAmovaLevels - 1]);
    fprintf(fp, "\n\n\n");
    fprintf(fp, "Variance components:\n\n\tVariance Component\tSigma^2\t\t%% of total");
    for (size_t i = 0; i < _ncoef - 1; i++) {
        fprintf(fp, "\n\t%-20s", metadata->levelNames[i + 1]);
        fprintf(fp, "\t%f", sigmasq[i]);
        fprintf(fp, "\t%f", pct_sigmasq[i]);
    }
    // Lowest level (i.e. Individual)
    fprintf(fp, "\n\t%-20s", metadata->levelNames[0]);
    fprintf(fp, "\t%f", sigmasq[_ncoef - 1]);
    fprintf(fp, "\t%f", pct_sigmasq[_ncoef - 1]);

    fprintf(fp, "\n\n\n");
    fprintf(fp, "\nVariance coefficients:\n\n\t");
    if (nLevels == 1) {
        fprintf(fp, "%f", ncoef[0]);
    } else {
        for (size_t i = 0; i < _ncoef; i++) {
            fprintf(fp, "%f\t", ncoef[i]);
        }
    }
    fprintf(fp, "\n\n\n");
    fprintf(fp, "Phi-statistic:\n\n");
    for (size_t i = 0; i < _phi; i++) {
        fprintf(fp, "\t%f", phi[i]);
    }
    fprintf(fp, "\n\n");
    fprintf(fp, "=============================================================================================");
    fprintf(fp, "\n\n");
}

// TODO add pct variance to output table and csv
// TODO add hdr to output
void amovaStruct::print_as_csv(metadataStruct *metadata) {
    // header
    //  type,label,value
    //  SSD,Among_region,0.1234
    //  fprintf(fp, "type,label,value\n");

    outFiles->out_amova_fs->kbuf = kbuf_init();
    kstring_t *kbuf = outFiles->out_amova_fs->kbuf;

    ksprintf(kbuf, "df,Total,%d\n", df[nAmovaLevels - 1]);
    ksprintf(kbuf, "SSD,Total,%f\n", ssd[nAmovaLevels - 1]);
    ksprintf(kbuf, "MSD,Total,%f\n", msd[nAmovaLevels - 1]);

    int x = 1;
    ksprintf(kbuf, "df,Among_%s_within_%s,%d\n", metadata->levelNames[x], "Total", df[0]);
    ksprintf(kbuf, "SSD,Among_%s_within_%s,%f\n", metadata->levelNames[x], "Total", ssd[0]);
    ksprintf(kbuf, "MSD,Among_%s_within_%s,%f\n", metadata->levelNames[x], "Total", msd[0]);
    while (x < metadata->nLevels + 1) {
        if (x == metadata->nLevels) {
            ksprintf(kbuf, "df,Among_%s_within_%s,%d\n", metadata->levelNames[0], metadata->levelNames[metadata->nLevels], df[x]);
            ksprintf(kbuf, "SSD,Among_%s_within_%s,%f\n", metadata->levelNames[0], metadata->levelNames[metadata->nLevels], ssd[x]);
            ksprintf(kbuf, "MSD,Among_%s_within_%s,%f\n", metadata->levelNames[0], metadata->levelNames[metadata->nLevels], msd[x]);
        } else {
            ksprintf(kbuf, "df,Among_%s_within_%s,%d\n", metadata->levelNames[x + 1], metadata->levelNames[x], df[x]);
            ksprintf(kbuf, "SSD,Among_%s_within_%s,%f\n", metadata->levelNames[x + 1], metadata->levelNames[x], ssd[x]);
            ksprintf(kbuf, "MSD,Among_%s_within_%s,%f\n", metadata->levelNames[x + 1], metadata->levelNames[x], msd[x]);
        }
        x++;
    }

    if (nLevels == 1) {
        ksprintf(kbuf, "Phi,%s_in_%s,%f\n", metadata->levelNames[1], "Total", phi[0]);
        ksprintf(kbuf, "Variance_coefficient,a,%f\n", ncoef[0]);
        ksprintf(kbuf, "Variance_component,%s,%f\n", metadata->levelNames[1], sigmasq[0]);
        ksprintf(kbuf, "Variance_component,%s,%f\n", metadata->levelNames[0], sigmasq[1]);
        ksprintf(kbuf, "Percentage_variance,%s,%f\n", metadata->levelNames[1], pct_sigmasq[0]);
        ksprintf(kbuf, "Percentage_variance,%s,%f\n", metadata->levelNames[0], pct_sigmasq[1]);
    } else if (nLevels == 2) {
        ksprintf(kbuf, "Phi,%s_in_%s,%f\n", metadata->levelNames[1], "Total", phi[0]);
        ksprintf(kbuf, "Phi,%s_in_%s,%f\n", metadata->levelNames[2], metadata->levelNames[1], phi[1]);
        ksprintf(kbuf, "Phi,%s_in_%s,%f\n", metadata->levelNames[2], "Total", phi[2]);
        ksprintf(kbuf, "Variance_coefficient,a,%f\n", ncoef[0]);
        ksprintf(kbuf, "Variance_coefficient,b,%f\n", ncoef[1]);
        ksprintf(kbuf, "Variance_coefficient,c,%f\n", ncoef[2]);
        ksprintf(kbuf, "Variance_component,%s,%f\n", metadata->levelNames[1], sigmasq[0]);
        ksprintf(kbuf, "Variance_component,%s,%f\n", metadata->levelNames[2], sigmasq[1]);
        ksprintf(kbuf, "Variance_component,%s,%f\n", metadata->levelNames[0], sigmasq[2]);
        ksprintf(kbuf, "Percentage_variance,%s,%f\n", metadata->levelNames[1], pct_sigmasq[0]);
        ksprintf(kbuf, "Percentage_variance,%s,%f\n", metadata->levelNames[2], pct_sigmasq[1]);
        ksprintf(kbuf, "Percentage_variance,%s,%f\n", metadata->levelNames[0], pct_sigmasq[2]);

    } else {
        fprintf(stderr, "[ERROR]: nLevels > 2 not supported yet\n");
    }

    outFiles->out_amova_fs->kbuf_write();
    ASSERT(outFiles->out_amova_fs->kbuf == NULL);  // TODO delme
}

amovaStruct::amovaStruct(metadataStruct *metadata) {
    // storing levels in order of highest to lowest level + total
    // {level1, level2, level3, ..., total}
    // where level1 is the highest level
    // (e.g. region in Individual ~ Region/Population/Subpopulation)

    nLevels = metadata->nLevels;
    ASSERT(nLevels > 0);

    // set number of amova levels
    // number of metadata levels + 2
    nAmovaLevels = nLevels + 2;
    _ncoef = nLevels + 1;
    _phi = nCk(_ncoef, 2);

    ssd = new double[nAmovaLevels];
    df = new int[nAmovaLevels];
    msd = new double[nAmovaLevels];
    for (size_t i = 0; i < (size_t)nAmovaLevels; i++) {
        ssd[i] = 0.0;
        df[i] = 0;
        msd[i] = 0.0;
    }

    ss = new double[_ncoef];
    ncoef = new double[_ncoef];
    sigmasq = new double[_ncoef];
    pct_sigmasq = new double[_ncoef];
    for (size_t i = 0; i < _ncoef; i++) {
        ncoef[i] = 0.0;
        sigmasq[i] = 0.0;
        ss[i] = 0.0;
        pct_sigmasq[i] = 0.0;
    }

    phi = new double[_phi];
    for (size_t i = 0; i < _phi; i++) {
        phi[i] = 0.0;
    }
}

amovaStruct::~amovaStruct() {
    DEL1D(df);
    DEL1D(ss);
    DEL1D(ssd);
    DEL1D(msd);
    DEL1D(ncoef);
    DEL1D(sigmasq);
    DEL1D(phi);
    DEL1D(pct_sigmasq);
}

double calculate_SumOfSquares_Total(distanceMatrixStruct *dm) {
    double ssd_TOTAL = 0.0;

    for (int px = 0; px < dm->nIndCmb; px++) {
        // the distance is already stored in squared form
        ssd_TOTAL += dm->M[px];
    }

    ssd_TOTAL = ssd_TOTAL / (double)dm->nInd;
    return ssd_TOTAL;
}

double calculate_SumOfSquares_Within(int lvl, amovaStruct *amova, distanceMatrixStruct *dm, metadataStruct *metadata) {
    ASSERT(lvl >= 0 && lvl < amova->nAmovaLevels);
    // if at the highest level, then calculate the total sum of squares
    if (lvl == amova->nAmovaLevels - 2) {
        return calculate_SumOfSquares_Total(dm);
    } else {
        double sum = 0.0;
        double ssd = 0.0;
        int n = 0;

        // for (int g = 0; g < metadata->hierArr[lvl]->nStrata; g++)
        for (int g = 0; g < metadata->nGroups[lvl]; g++) {
            sum = 0.0;
            n = metadata->countIndsInGroup(lvl, g);

            // only use the individuals that belong to this group (localGrpIdx:g at level:lvl)
            for (int i1 = 0; i1 < dm->nInd - 1; i1++) {
                for (int i2 = i1 + 1; i2 < dm->nInd; i2++) {
                    if (metadata->indFromGroup(i1, lvl, g) && metadata->indFromGroup(i2, lvl, g)) {
                        sum += dm->M[dm->inds2idx[i1][i2]] / n;
                    }
                }
            }

            ssd += sum;
        }
        return ssd;
    }
}
