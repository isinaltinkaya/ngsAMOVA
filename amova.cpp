#include "amova.h"

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
    fprintf(fp, "Variance components:\n\n");
    for (size_t i = 0; i < _ncoef - 1; i++) {
        fprintf(fp, "\n\t%-20s", metadata->levelNames[i + 1]);
        fprintf(fp, "\t%f", sigmasq[i]);
    }
    // Lowest level (i.e. Individual)
    fprintf(fp, "\n\t%-20s", metadata->levelNames[0]);
    fprintf(fp, "\t%f", sigmasq[_ncoef - 1]);

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

void amovaStruct::print_as_csv(FILE *fp, metadataStruct *metadata) {
    // TODO add percentage total?
    // header
    //  type,label,value
    //  SSD,Among_region,0.1234
    //  fprintf(fp, "type,label,value\n");
    ASSERT(fp != NULL);
    fprintf(fp, "df,Total,%d\n", df[nAmovaLevels - 1]);
    fprintf(fp, "SSD,Total,%f\n", ssd[nAmovaLevels - 1]);
    fprintf(fp, "MSD,Total,%f\n", msd[nAmovaLevels - 1]);

    int x = 1;
    fprintf(fp, "df,Among_%s_within_%s,%d\n", metadata->levelNames[x], "Total", df[0]);
    fprintf(fp, "SSD,Among_%s_within_%s,%f\n", metadata->levelNames[x], "Total", ssd[0]);
    fprintf(fp, "MSD,Among_%s_within_%s,%f\n", metadata->levelNames[x], "Total", msd[0]);
    while (x < metadata->nLevels + 1) {
        if (x == metadata->nLevels) {
            fprintf(fp, "df,Among_%s_within_%s,%d\n", metadata->levelNames[0], metadata->levelNames[metadata->nLevels], df[x]);
            fprintf(fp, "SSD,Among_%s_within_%s,%f\n", metadata->levelNames[0], metadata->levelNames[metadata->nLevels], ssd[x]);
            fprintf(fp, "MSD,Among_%s_within_%s,%f\n", metadata->levelNames[0], metadata->levelNames[metadata->nLevels], msd[x]);
        } else {
            fprintf(fp, "df,Among_%s_within_%s,%d\n", metadata->levelNames[x + 1], metadata->levelNames[x], df[x]);
            fprintf(fp, "SSD,Among_%s_within_%s,%f\n", metadata->levelNames[x + 1], metadata->levelNames[x], ssd[x]);
            fprintf(fp, "MSD,Among_%s_within_%s,%f\n", metadata->levelNames[x + 1], metadata->levelNames[x], msd[x]);
        }
        x++;
    }

    if (nLevels == 1) {
        fprintf(fp, "Phi,%s_in_%s,%f\n", metadata->levelNames[1], "Total", phi[0]);
        fprintf(fp, "Variance_coefficient,a,%f\n", ncoef[0]);
        fprintf(fp, "Variance_component,%s,%f\n", metadata->levelNames[1], sigmasq[0]);
        fprintf(fp, "Variance_component,%s,%f\n", metadata->levelNames[0], sigmasq[1]);
    } else if (nLevels == 2) {
        fprintf(fp, "Phi,%s_in_%s,%f\n", metadata->levelNames[1], "Total", phi[0]);
        fprintf(fp, "Phi,%s_in_%s,%f\n", metadata->levelNames[2], metadata->levelNames[1], phi[1]);
        fprintf(fp, "Phi,%s_in_%s,%f\n", metadata->levelNames[2], "Total", phi[2]);
        fprintf(fp, "Variance_coefficient,a,%f\n", ncoef[0]);
        fprintf(fp, "Variance_coefficient,b,%f\n", ncoef[1]);
        fprintf(fp, "Variance_coefficient,c,%f\n", ncoef[2]);
        fprintf(fp, "Variance_component,%s,%f\n", metadata->levelNames[1], sigmasq[0]);
        fprintf(fp, "Variance_component,%s,%f\n", metadata->levelNames[2], sigmasq[1]);
        fprintf(fp, "Variance_component,%s,%f\n", metadata->levelNames[0], sigmasq[2]);
    } else {
        fprintf(stderr, "[ERROR]: nLevels > 2 not supported yet\n");
    }
}

amovaStruct::amovaStruct(metadataStruct *metadata) {
    // storing levels in order of highest to lowest level + total
    // {level1, level2, level3, ..., total}
    // where level1 is the highest level
    // (e.g. region in Individual ~ Region/Population/Subpopulation)

    nLevels = metadata->nLevels;
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
    for (size_t i = 0; i < _ncoef; i++) {
        ncoef[i] = 0.0;
        sigmasq[i] = 0.0;
        ss[i] = 0.0;
    }

    phi = new double[_phi];
    for (size_t i = 0; i < _phi; i++) {
        phi[i] = 0.0;
    }
}

amovaStruct::~amovaStruct() {
    delete[] ssd;
    delete[] ss;
    delete[] df;
    delete[] msd;
    delete[] ncoef;
    delete[] phi;
    delete[] sigmasq;
}

double calculate_SumOfSquares_Total(distanceMatrixStruct *dMS) {
    double ssd_TOTAL = 0.0;

    for (int px = 0; px < dMS->nIndCmb; px++) {
        // the distance is already stored in squared form
        ssd_TOTAL += dMS->M[px];
    }

    ssd_TOTAL = ssd_TOTAL / (double)dMS->nInd;
    return ssd_TOTAL;
}

double calculate_SumOfSquares_Within(int lvl, amovaStruct *aS, distanceMatrixStruct *dMS, metadataStruct *metadata, paramStruct *pars) {
    ASSERT(lvl >= 0 && lvl < aS->nAmovaLevels);
    // if at the highest level, then calculate the total sum of squares
    if (lvl == aS->nAmovaLevels - 2) {
        return calculate_SumOfSquares_Total(dMS);
    } else {
        double sum = 0.0;
        double ssd = 0.0;
        int n = 0;

        // for (int g = 0; g < metadata->hierArr[lvl]->nStrata; g++)
        for (int g = 0; g < metadata->nGroups[lvl]; g++) {
            sum = 0.0;
            n = metadata->countIndsInGroup(lvl, g);

            // only use the individuals that belong to this group (localGrpIdx:g at level:lvl)
            for (int i1 = 0; i1 < dMS->nInd - 1; i1++) {
                for (int i2 = i1 + 1; i2 < dMS->nInd; i2++) {
                    if (metadata->indFromGroup(i1, lvl, g) && metadata->indFromGroup(i2, lvl, g)) {
                        sum += dMS->M[dMS->inds2idx[i1][i2]] / n;
                    }
                }
            }

            ssd += sum;
        }
        return ssd;
    }
}
/// @brief doAMOVA - perform AMOVA analysis
/// @param dMS  pointer to distanceMatrixStruct
/// @param metadata   pointer to metadataStruct
/// @param pars pointer to paramStruct
/// @return
amovaStruct *doAmova(distanceMatrixStruct *dMS, metadataStruct *metadata, paramStruct *pars) {
    amovaStruct *aS = new amovaStruct(metadata);

    // DEGREES OF FREEDOM
    // ----------------------------------------------- //
    // calculate the degrees of freedom for each level

    // highest level df is nStrata - 1
    // e.g. Individual ~ Region / Population / Subpopulation
    // if there are 3 regions, then nStrata is number of unique regions - 1 = 3-1 = 2
    aS->df[0] = metadata->nGroups[0] - 1;

    // total df
    aS->df[aS->nAmovaLevels - 1] = dMS->nInd - 1;

    if (metadata->nLevels == 1) {
        // df for lowest amova level
        // e.g. Individual ~ Region / Population / Subpopulation
        // nInd - sum(for each Population, nStrata (num_Subpopulation) in Population)
        aS->df[1] = dMS->nInd - metadata->nGroups[0];
    } else if (metadata->nLevels == 2) {
        int sum0 = 0;
        for (int i = 0; i < metadata->nGroups[0]; i++) {
            sum0 += metadata->countNSubgroupAtLevel(0, i, 1);
        }
        int sum1 = metadata->nGroups[0];
        aS->df[1] = sum0 - sum1;
        aS->df[2] = dMS->nInd - sum0;
    } else {
        fprintf(stderr, "[ERROR] NOT IMPLEMENTED YET");
        ASSERT(0 == 1);
    }

    // ----------------------------------------------- //
    // SUM OF SQUARES
    // ----------------------------------------------- //
    // calculate the sum of squares within for each level
    //
    for (size_t i = 0; i < aS->_ncoef; i++) {
        // within level sum of squares
        aS->ss[i] = calculate_SumOfSquares_Within(i, aS, dMS, metadata, pars);
    }

    // ----------------------------------------------- //
    // SUM OF SQUARED DEVIATIONS
    // ----------------------------------------------- //

    if (metadata->nLevels == 1) {
        aS->ssd[0] = aS->ss[1] - aS->ss[0];
        aS->ssd[1] = aS->ss[0];
        aS->ssd[2] = aS->ss[1];
    } else if (metadata->nLevels == 2) {
        aS->ssd[3] = aS->ss[2];
        aS->ssd[1] = aS->ss[0] - aS->ss[1];
        aS->ssd[2] = aS->ss[1];
        aS->ssd[0] = aS->ssd[3] - (aS->ssd[0] + aS->ssd[1] + aS->ssd[2]);
    } else {
        fprintf(stderr, "[ERROR] NOT IMPLEMENTED YET");
        ASSERT(0 == 1);
    }

    // get msd from ssd and df
    for (size_t i = 0; i < (size_t)aS->nAmovaLevels; i++) {
        aS->msd[i] = aS->ssd[i] / aS->df[i];
    }

    if (metadata->nLevels == 1) {
        // ----------------------------------------------- //
        // VARIANCE COEFFICIENTS
        // ----------------------------------------------- //
        // calculate variance coefficients (n) for 1 level AMOVA

        double n_gi = 0.0;

        // n = [ N - \sum_{g \in G} ( N^2_{g}/N) ) ]  /   G - 1
        for (int sti = 0; sti < metadata->nGroups[0]; sti++) {
            ASSERT(metadata->nIndPerStrata[0] != 0);
            n_gi += SQUARE(metadata->nIndPerStrata[0][sti]);
        }
        n_gi = n_gi / (double)dMS->nInd;
        aS->ncoef[0] = (double)((double)dMS->nInd - (double)n_gi) / (double)(metadata->nGroups[0] - 1);

        // ----------------------------------------------- //
        // VARIANCE COMPONENTS
        // ----------------------------------------------- //
        // calculate variance components (sigma squared)
        // for 1 level AMOVA
        aS->sigmasq_total = aS->sigmasq[0] + aS->sigmasq[1];

        aS->sigmasq[0] = (aS->msd[0] - aS->msd[1]) / aS->ncoef[0];
        aS->sigmasq[1] = aS->msd[1];

        aS->phi[0] = aS->sigmasq[0] / (aS->sigmasq[0] + aS->sigmasq[1]);
    } else if (metadata->nLevels == 2) {
        // ----------------------------------------------- //
        // VARIANCE COEFFICIENTS
        // ----------------------------------------------- //
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

        aS->ncoef[0] = (NUL - NUR) / NL;

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

        aS->ncoef[1] = (NUR - N1UR) / N1L;

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

        aS->ncoef[2] = (NUL - N2UR) / N1L;

        // ----------------------------------------------- //
        // VARIANCE COMPONENTS
        // ----------------------------------------------- //
        // calculate variance components (sigma squared)
        aS->sigmasq[2] = aS->msd[2];
        aS->sigmasq[1] = (aS->msd[1] - aS->sigmasq[2]) / aS->ncoef[0];
        aS->sigmasq[0] = (aS->msd[0] - aS->msd[2] - (aS->ncoef[1] * aS->sigmasq[1])) / aS->ncoef[2];
        aS->sigmasq_total = aS->sigmasq[0] + aS->sigmasq[1] + aS->sigmasq[2];

        // ----------------------------------------------- //
        // PHI STATISTIC
        // ----------------------------------------------- //
        // calculate phi statistics
        //
        // Individual,Region,Population,Total
        //

        // lvl0 (i.e. reg) in TOTAL
        // Phi_CT
        if (aS->sigmasq_total <= 0) {
            fprintf(stderr, "\n[ERROR]\tTotal variance is %f (sigmasq_total is <= 0). Please check your analysis settings and make sure you have enough data.\n", aS->sigmasq_total);
            exit(1);
        }

        aS->phi[0] = aS->sigmasq[0] / (aS->sigmasq_total);

        // lvl1 in lvl0
        // Phi_SC
        aS->phi[1] = aS->sigmasq[1] / (aS->sigmasq[1] + aS->sigmasq[2]);

        // lvl1 (i.e. pop) in TOTAL
        // Phi_ST
        aS->phi[2] = (aS->sigmasq[0] + aS->sigmasq[1]) / (aS->sigmasq_total);
    }

    IO::vprint(stderr, 1, "Finished running AMOVA\n");

    return aS;
}
