#include "vcfReader.h"

#include "bootstrap.h"

int glData::ind_data_isMissing(const int ind_i) {
    float *ind_data = ind_ptr(ind_i);
    if (ind_data == NULL) {
        NEVER;
        return 1;
    }

    if (isnan(ind_data[0])) {
        return 1;
    }

#if DEBUG
    for (int i = 1; i < n_gls; ++i) {
        if (1 == isnan(ind_data[i]))
            ERROR("nan found after a non-nan value at ind_i: %d, i: %d", ind_i, i);
    }
#endif

    int z = 0;
    for (int i = 0; i < n_gls; ++i) {
        if (0 == ind_data[i])
            z++;
    }
    if (z == n_gls)
        return 1;

    return 0;
}

glData::glData(vcfData *vcfd) {
    size_e = 0;
    n_values = bcf_get_format_float(vcfd->hdr, vcfd->bcf, "GL", &data, &size_e);
    if (n_values <= 0) {
        ERROR("Could not read GL tag from the VCF file.");
    }

    n_gls = n_values / vcfd->nInd;

    if (10 == n_gls) {
        IO::vprint(3, "GL field has 10 values per individual.");
    } else if (3 == n_gls) {
        IO::vprint(3, "GL field has 3 values per individual.");
    } else {
        ERROR("GL field has %d values per individual. Only 3 or 10 are supported.", n_gls);
    }
}

glData::~glData() {
    FREE(data);
}

float *glData::ind_ptr(const int ind_i) {
    return (data + ind_i * n_gls);
}

gtData::gtData(vcfData *vcfd) {
    size_e = 0;
    n_values = bcf_get_genotypes(vcfd->hdr, vcfd->bcf, &data, &size_e);
    if (n_values <= 0) {
        ERROR("GT not present.");
    }

    ploidy = n_values / vcfd->nInd;
    if (ploidy != 2) {
        ERROR("Ploidy %d not supported.", ploidy);
    }
}

gtData::~gtData() {
    FREE(data);
}

int *gtData::ind_ptr(const int ind_i) {
    return (data + ind_i * ploidy);
}

// from angsd analysisFunction.cpp
extern const int bcf_allele_charToInt[(int)256] = {
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 15
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 31
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 47
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 63
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,  // 79
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 95
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,  // 111
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 127
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 143
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 159
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 175
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 191
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 207
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 223
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 239
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4   // 255
};

int bcf_alleles_get_gtidx(int a1, int a2) {
    return bcf_alleles2gt(a1, a2);
}

int bcf_alleles_get_gtidx(char a1, char a2) {
    return bcf_alleles2gt(bcf_allele_charToInt[(int)(unsigned char)a1], bcf_allele_charToInt[(int)(unsigned char)a2]);
}

int bcf_alleles_get_gtidx(unsigned char a1, unsigned char a2) {
    return bcf_alleles2gt(bcf_allele_charToInt[(int)a1], bcf_allele_charToInt[(int)a2]);
}

// return 1: skip site for all individuals
int site_read_GL(const int contig_i, const int site_i, vcfData *vcfd, paramStruct *pars, pairStruct **pairs) {
    const int nInd = pars->nInd;

    // data from VCF GL field is in log10 scale thus [l]gls
    glData lgls(vcfd);

    int cmbArr[pars->nIndCmb];
    for (int i = 0; i < pars->nIndCmb; i++) {
        cmbArr[i] = 0;
    }

    if (1 == bcf_is_snp(vcfd->bcf)) {
        char a1;
        char a2;

        if (NULL != pars->ancder) {
            a1 = pars->ancder->a1[contig_i][site_i];
            a2 = pars->ancder->a2[contig_i][site_i];
        } else if (NULL != pars->majmin) {
            a1 = pars->majmin->a1[contig_i][site_i];
            a2 = pars->majmin->a2[contig_i][site_i];
        } else if (1 == args->isSim) {
            // reference allele
            a1 = bcf_allele_charToInt[(int)(unsigned char)vcfd->bcf->d.allele[0][0]];

            // alternative allele 1
            a2 = bcf_allele_charToInt[(int)(unsigned char)vcfd->bcf->d.allele[1][0]];
        } else {
            NEVER;
        }

        // ASSERT(a1 != a2); debug
        // fprintf(stderr, "a1: %c, a2: %c\n", vcfd->bcf->d.allele[0][0], vcfd->bcf->d.allele[1][0]);

        for (int indi = 0; indi < nInd; indi++) {
            const int lngls_ind_start = indi * vcfd->nGT;

            if (1 == lgls.ind_data_isMissing(indi)) {
                // if only use sites shared across all individuals (minInd 0); skip site when you first encounter nan
                if (args->minInd == 0) {
                    return 1;
                } else {
                    // if there are only 2 individuals any missing will skip the site
                    if (nInd == 2) {
                        return 1;
                    }

                    lgls.n_missing_ind++;

                    if (nInd == lgls.n_missing_ind) {
                        return 1;
                    }

                    if (args->minInd != 2) {
                        // skip site if minInd is defined and #non-missing inds=<nInd
                        if ((nInd - lgls.n_missing_ind) < args->minInd) {
                            // fprintf(stderr,"\n\nMinimum number of individuals -minInd is set to %d, but nInd-n_missing_ind==n_nonmissing_ind is %d at site %d\n\n",args->minInd,pars->nInd-n_missing_ind,site);
                            return 1;
                        }
                    }
                }
            } else {
                int pidx = -1;
                // if not first individual, check previous individuals pairing with current ind
                if (indi != 0) {
                    for (int indi2 = indi - 1; indi2 > -1; indi2--) {
                        // both inds has data
                        pidx = nCk_idx(nInd, indi, indi2);

                        if (cmbArr[pidx]) {
                            // append site to sharedSites shared sites list
                            pairs[pidx]->sharedSites_add(pars->nSites);
                        }
                    }
                }

                // if not last individual, check latter individuals pairing with current ind
                if (indi != nInd - 1) {
                    for (int indi2 = indi + 1; indi2 < nInd; indi2++) {
                        pidx = nCk_idx(nInd, indi, indi2);
                        cmbArr[pidx]++;
                    }
                }
                // TODO can this prev and latter check be checkin stuff twice?

                vcfd->lngl[pars->nSites][lngls_ind_start + 0] = (double)LOG2LN(lgls.ind_ptr(indi)[bcf_alleles_get_gtidx(a1, a1)]);
                vcfd->lngl[pars->nSites][lngls_ind_start + 1] = (double)LOG2LN(lgls.ind_ptr(indi)[bcf_alleles_get_gtidx(a1, a2)]);
                vcfd->lngl[pars->nSites][lngls_ind_start + 2] = (double)LOG2LN(lgls.ind_ptr(indi)[bcf_alleles_get_gtidx(a2, a2)]);
                // fprintf(stderr, "site: %d, indi: %d, a1: %c, a2: %c, lngl: %f, %f, %f, log10gl: %f, %f, %f\n", pars->nSites, indi, a1, a2, vcfd->lngl[pars->nSites][lngls_ind_start + 0], vcfd->lngl[pars->nSites][lngls_ind_start + 1], vcfd->lngl[pars->nSites][lngls_ind_start + 2], lgls.ind_ptr(indi)[bcf_alleles_get_gtidx(a1, a1)], lgls.ind_ptr(indi)[bcf_alleles_get_gtidx(a1, a2)], lgls.ind_ptr(indi)[bcf_alleles_get_gtidx(a2, a2)]);
            }
        }
    } else {
        // TODO skip nonsnp data?
        fprintf(stderr, "\n\nbcf_IS_SNP==0!!!\n\n");
        exit(1);
    }

    return 0;
}

int get_JointGenoDist_GT(const int contig_i, const int site_i, vcfData *vcfd, paramStruct *pars, blobStruct *blob, const int block_i) {
    const int nInd = pars->nInd;

    if (1 == bcf_is_snp(vcfd->bcf)) {
        char a1;
        char a2;

        if (NULL != pars->ancder) {
            a1 = bcf_allele_charToInt[(int)pars->ancder->a1[contig_i][site_i]];
            a2 = bcf_allele_charToInt[(int)pars->ancder->a2[contig_i][site_i]];
        } else if (NULL != pars->majmin) {
            a1 = bcf_allele_charToInt[(int)pars->majmin->a1[contig_i][site_i]];
            a2 = bcf_allele_charToInt[(int)pars->majmin->a2[contig_i][site_i]];
        } else if (1 == args->isSim) {
            // reference allele
            a1 = bcf_allele_charToInt[(int)(unsigned char)vcfd->bcf->d.allele[0][0]];
            // alternative allele 1
            a2 = bcf_allele_charToInt[(int)(unsigned char)vcfd->bcf->d.allele[1][0]];
        } else {
            NEVER;
        }

        // ACGT
        int alleles_lut[4] = {-1, -1, -1, -1};
        alleles_lut[(int)a1] = 0;  // major/ancestral
        alleles_lut[(int)a2] = 1;  // minor/derived
        // rest is -1 == not ancestral or derived

        gtData gts(vcfd);

        // --- minInd threshold check
        for (int i1 = 0; i1 < nInd; ++i1) {
            int *ptr1 = gts.ind_ptr(i1);

            if (1 == bcf_gt_is_missing(ptr1[0]) || 1 == bcf_gt_is_missing(ptr1[1])) {
                // if arg=only use sites shared across all individuals (minInd 0); skip site when you first encounter nan
                if (args->minInd == 0) {
                    return 1;
                }

                gts.n_missing_ind++;
                if (nInd == gts.n_missing_ind) {
                    return 1;
                }

                // skip site if minInd is defined and #non-missing inds=<nInd
                if (args->minInd != 2 && (nInd - gts.n_missing_ind) < args->minInd) {
                    return 1;
                }
                continue;
            }
        }

        for (int i1 = 0; i1 < nInd; ++i1) {
            int *ptr1 = gts.ind_ptr(i1);

            if (1 == bcf_gt_is_missing(ptr1[0]) || 1 == bcf_gt_is_missing(ptr1[1]))
                continue;

            // index in the bcf->d.allele array
            int i1_1 = bcf_gt_allele(ptr1[0]);
            int i1_2 = bcf_gt_allele(ptr1[1]);

            // corresponding allele char
            char *i1a1 = vcfd->bcf->d.allele[i1_1];
            char *i1a2 = vcfd->bcf->d.allele[i1_2];

            // index in ACGT
            int i1a1i = bcf_allele_charToInt[(int)*i1a1];
            int i1a2i = bcf_allele_charToInt[(int)*i1a2];

            int i1a1_state = alleles_lut[i1a1i];
            int i1a2_state = alleles_lut[i1a2i];

            // if any allele in i1's gt is not in the allelic state list, skip
            if (i1a1_state == -1 && i1a2_state == -1) {
                continue;
            }

            int i1_state_gt = i1a1_state + i1a2_state;

            for (int i2 = i1 + 1; i2 < nInd; i2++) {
                int *ptr2 = gts.ind_ptr(i2);

                if (1 == bcf_gt_is_missing(ptr2[0]) || 1 == bcf_gt_is_missing(ptr2[1]))
                    continue;

                // index in the bcf->d.allele array
                int i2_1 = bcf_gt_allele(ptr2[0]);
                int i2_2 = bcf_gt_allele(ptr2[1]);

                // corresponding allele char
                char *i2a1 = vcfd->bcf->d.allele[i2_1];
                char *i2a2 = vcfd->bcf->d.allele[i2_2];

                // index in ACGT
                int i2a1i = bcf_allele_charToInt[(int)*i2a1];
                int i2a2i = bcf_allele_charToInt[(int)*i2a2];

                int i2a1_state = alleles_lut[i2a1i];
                int i2a2_state = alleles_lut[i2a2i];

                // if any allele in i2's gt is not in the allelic state list, skip
                if (i2a1_state == -1 && i2a2_state == -1) {
                    continue;
                }

                int pair_idx = nCk_idx(nInd, i1, i2);

                int i2_state_gt = i2a1_state + i2a2_state;

                // vcfd->JointGenoCountDistGT[pair_idx][get_3x3_idx[i1_state_gt][i2_state_gt]]++;
                // vcfd->JointGenoCountDistGT[pair_idx][9]++;  // #shared sites

                vcfd->jgcd_gt[block_i][pair_idx][get_3x3_idx[i1_state_gt][i2_state_gt]]++;
                vcfd->pair_shared_nSites[block_i][pair_idx]++;
            }
        }

    } else {
        // TODO skip nonsnp data?
        fprintf(stderr, "\n\nbcf_IS_SNP==0!!!\n\n");
        exit(1);
    }

    return 0;
}

// // TODO this would only work with my simulated data, consider removing
// see the same in get_jointgenodist_gl
// a1 = bcf_allele_charToInt[(int)(unsigned char)vcfd->bcf->d.allele[0][0]];
// a2 = bcf_allele_charToInt[(int)(unsigned char)vcfd->bcf->d.allele[1][0]];
// ASSERT(a1 != a2);

int get_JointGenoDist_GT(const int contig_i, const int site_i, vcfData *vcfd, paramStruct *pars) {
    const int nInd = pars->nInd;

    ASSERT(NULL != pars->ancder->a1[contig_i]);
    // assume: ancderfile contains all contigs in vcf
    const char a1 = bcf_allele_charToInt[(int)pars->ancder->a1[contig_i][site_i]];
    const char a2 = bcf_allele_charToInt[(int)pars->ancder->a2[contig_i][site_i]];
    // ACGT
    int alleles_lut[4] = {-1, -1, -1, -1};
    alleles_lut[(int)a1] = 0;  // major/ancestral
    alleles_lut[(int)a2] = 1;  // minor/derived
    // rest is -1 == not ancestral or derived

    gtData gts(vcfd);

    for (int i1 = 0; i1 < nInd; ++i1) {
        int *ptr1 = gts.ind_ptr(i1);

        if (1 == bcf_gt_is_missing(ptr1[0]) || 1 == bcf_gt_is_missing(ptr1[1])) {
            // if arg=only use sites shared across all individuals (minInd 0); skip site when you first encounter nan
            if (args->minInd == 0) {
                return 1;
            }

            gts.n_missing_ind++;
            if (nInd == gts.n_missing_ind) {
                return 1;
            }

            // skip site if minInd is defined and #non-missing inds=<nInd
            if (args->minInd != 2 && (nInd - gts.n_missing_ind) < args->minInd) {
                return 1;
            }
            continue;
        }
    }

    for (int i1 = 0; i1 < nInd; ++i1) {
        int *ptr1 = gts.ind_ptr(i1);

        if (1 == bcf_gt_is_missing(ptr1[0]) || 1 == bcf_gt_is_missing(ptr1[1]))
            continue;

        // index in the bcf->d.allele array
        int i1_1 = bcf_gt_allele(ptr1[0]);
        int i1_2 = bcf_gt_allele(ptr1[1]);

        // corresponding allele char
        char *i1a1 = vcfd->bcf->d.allele[i1_1];
        char *i1a2 = vcfd->bcf->d.allele[i1_2];

        // index in ACGT
        int i1a1i = bcf_allele_charToInt[(int)*i1a1];
        int i1a2i = bcf_allele_charToInt[(int)*i1a2];

        int i1a1_state = alleles_lut[i1a1i];
        int i1a2_state = alleles_lut[i1a2i];

        // if any allele in i1's gt is not in the allelic state list, skip
        if (i1a1_state == -1 && i1a2_state == -1) {
            continue;
        }

        int i1_state_gt = i1a1_state + i1a2_state;

        for (int i2 = i1 + 1; i2 < nInd; i2++) {
            int *ptr2 = gts.ind_ptr(i2);

            if (1 == bcf_gt_is_missing(ptr2[0]) || 1 == bcf_gt_is_missing(ptr2[1]))
                continue;

            // index in the bcf->d.allele array
            int i2_1 = bcf_gt_allele(ptr2[0]);
            int i2_2 = bcf_gt_allele(ptr2[1]);

            // corresponding allele char
            char *i2a1 = vcfd->bcf->d.allele[i2_1];
            char *i2a2 = vcfd->bcf->d.allele[i2_2];

            // index in ACGT
            int i2a1i = bcf_allele_charToInt[(int)*i2a1];
            int i2a2i = bcf_allele_charToInt[(int)*i2a2];

            int i2a1_state = alleles_lut[i2a1i];
            int i2a2_state = alleles_lut[i2a2i];

            // if any allele in i2's gt is not in the allelic state list, skip
            if (i2a1_state == -1 && i2a2_state == -1) {
                continue;
            }

            int pair_idx = nCk_idx(nInd, i1, i2);

            int i2_state_gt = i2a1_state + i2a2_state;

            vcfd->JointGenoCountDistGT[pair_idx][get_3x3_idx[i1_state_gt][i2_state_gt]]++;
            vcfd->JointGenoCountDistGT[pair_idx][9]++;  // #shared sites
        }
    }

    return 0;
}

void vcfData::set_n_joint_categories(paramStruct *pars) {
    if (1 == args->doDist) {
        // EM using 3 GL values
        if (args->doEM == 1) {
            nGT = 3;
            nJointClasses = 9;
        }
        // EM using 10 GL values
        else if (args->doEM == 2) {
            nGT = 10;
            nJointClasses = 100;
        }
    } else if (2 == args->doDist) {
        // GT
        nGT = 3;
        nJointClasses = 9;
    } else {
        NEVER;
    }
}

vcfData *vcfData_init(paramStruct *pars) {
    vcfData *vcfd = new vcfData;

    vcfd->in_fp = bcf_open(args->in_vcf_fn, "r");
    if (vcfd->in_fp == NULL) {
        fprintf(stderr, "\n[ERROR] Could not open bcf file: %s\n", args->in_vcf_fn);
        exit(1);
    }

    vcfd->set_n_joint_categories(pars);

    vcfd->hdr = bcf_hdr_read(vcfd->in_fp);
    vcfd->bcf = bcf_init();

    if (NULL != args->in_region) {
        // vcfd->tbx = IO::load_vcf_tabix_idx(args->in_vcf_fn);
        // vcfd->itr = tbx_itr_querys(vcfd->tbx, args->in_region);

        vcfd->idx = IO::load_bcf_csi_idx(args->in_vcf_fn);

        vcfd->itr = bcf_itr_querys(vcfd->idx, vcfd->hdr, args->in_region);
        if (NULL == vcfd->itr) {
            ERROR("Could not parse region: %s. Please make sure region exists and defined in the correct format.", args->in_region);
        }
        ASSERT(vcfd->itr != NULL);
    }

    if (NULL != args->in_regions_bed_fn) {
        ERROR("Not implemented yet");
    }
    if (NULL != args->in_regions_tab_fn) {
        ERROR("Not implemented yet");
    }

    pars->nInd = bcf_hdr_nsamples(vcfd->hdr);
    // TODO unnecessary
    vcfd->nInd = pars->nInd;
    pars->nIndCmb = nC2(pars->nInd);
    vcfd->nIndCmb = pars->nIndCmb;

    if (1 == args->doDist) {
        vcfd->lngl_init(args->doEM);
        vcfd->init_JointGenoCountDistGL();
        vcfd->init_JointGenoProbDistGL();
    } else if (2 == args->doDist) {
        vcfd->init_JointGenoCountDistGT();
    }

    fprintf(stderr, "\nNumber of individual pairs: %d\n", pars->nIndCmb);

    vcfd->addIndNames();

    vcfd->nContigs = vcfd->hdr->n[BCF_DT_CTG];

    check_consistency_args_pars(args, pars);

    return vcfd;
}

void vcfData_destroy(vcfData *v) {
    bcf_hdr_destroy(v->hdr);
    bcf_destroy(v->bcf);

    int BCF_CLOSE = bcf_close(v->in_fp);
    if (0 != BCF_CLOSE) {
        fprintf(stderr, "\n[ERROR]\tbcf_close had non-zero status %d\n", BCF_CLOSE);
        exit(BCF_CLOSE);
    }

    if (v->lngl != NULL) {
        for (size_t i = 0; i < v->_lngl; i++) {
            FREE(v->lngl[i]);
        }
        FREE(v->lngl);
    }

    if (v->JointGenoCountDistGT != NULL) {
        for (size_t s = 0; s < (size_t)v->nIndCmb; s++) {
            FREE(v->JointGenoCountDistGT[s]);
        }
        FREE(v->JointGenoCountDistGT);
    }
    if (v->JointGenoCountDistGL != NULL) {
        for (size_t s = 0; s < (size_t)v->nIndCmb; s++) {
            FREE(v->JointGenoCountDistGL[s]);
        }
        FREE(v->JointGenoCountDistGL);
    }
    if (v->JointGenoProbDistGL != NULL) {
        for (size_t s = 0; s < (size_t)v->nIndCmb; s++) {
            FREE(v->JointGenoProbDistGL[s]);
        }
        FREE(v->JointGenoProbDistGL);
    }
    for (int i = 0; i < v->nInd; i++) {
        FREE(v->indNames[i]);
    }
    FREE(v->indNames);

    if (NULL != v->itr) {
        hts_itr_destroy(v->itr);
    }

    if (NULL != v->idx) {
        hts_idx_destroy(v->idx);
    }

    delete v;
}

void readSites_GL(vcfData *vcfd, paramStruct *pars, pairStruct **pairSt) {
    int skip_site = 0;
    int contig_i = 0;
    int site_i = 0;
    char prev_contig[100];

    while (vcfd->records_next()) {
        while (pars->nSites >= vcfd->_lngl) {
            vcfd->lngl_expand();
        }

        // if cleared in the previous loop to skip the site, allocate memory again
        if (NULL == vcfd->lngl[pars->nSites]) {
            vcfd->lngl[pars->nSites] = (double *)malloc(pars->nInd * vcfd->nGT * sizeof(double));
            for (int indi = 0; indi < pars->nInd; indi++) {
                const int lngls_ind_start = indi * vcfd->nGT;
                for (int j = 0; j < vcfd->nGT; j++) {
                    vcfd->lngl[pars->nSites][lngls_ind_start + j] = NEG_INF;
                }
            }
        }

        const char *contig_id = bcf_seqname(vcfd->hdr, vcfd->bcf);
        if (0 != contig_i && 0 != strcmp(contig_id, prev_contig)) {
            strcpy(prev_contig, contig_id);
            ++contig_i;
            if (contig_i >= vcfd->nContigs) {
                ERROR("Number of contigs in the VCF file header is larger than the number of contigs in the VCF file");
            }
            site_i = 0;
        }

        skip_site = site_read_GL(contig_i, site_i, vcfd, pars, pairSt);

        if (skip_site == 1) {
            fprintf(stderr, "\n->\tSkipping site %lu for all individuals\n\n", pars->totSites);
            pars->totSites++;

            //  next loop will skip this and use the same pars->nSites
            FREE(vcfd->lngl[pars->nSites]);
            continue;
        }

        // check if contig changed
        pars->nSites++;
        ++site_i;
        pars->totSites++;
    }

    vcfd->nSites = pars->nSites;
    vcfd->totSites = pars->totSites;

    fprintf(stderr, "\n\t-> Finished reading sites\n");
}

void readSites_GT(vcfData *vcfd, paramStruct *pars, pairStruct **pairSt) {
    int skip_site = 0;
    /*
     * [START] Reading sites
     *
     * nSites=0; totSites=0
     *
     * if minInd is set
     * 		nSites==where minInd threshold can be passed
     * 		totSites==all sites processed
     *
     */

    ASSERT(0 == pars->nSites);
    ASSERT(0 == pars->nContigs);

    int contig_i = 0;
    int site_i = 0;

    while (vcfd->records_next()) {
        skip_site = get_JointGenoDist_GT(contig_i, site_i, vcfd, pars);

        if (skip_site == 1) {
            fprintf(stderr, "\n->\tSkipping site %lu for all individuals\n\n", pars->totSites);
            pars->totSites++;

            continue;
        }

        ++site_i;
        pars->nSites++;
        pars->totSites++;
    }

    vcfd->nSites = pars->nSites;
    vcfd->totSites = pars->totSites;
#if 0
	fprintf(stderr,"\n\n\t-> Printing at site %d\n\n",pars->nSites);
	fprintf(stderr,"\nPrinting at (idx: %lu, pos: %lu 1based %lu) totSites:%d\n\n",pars->nSites,vcfd->bcf->pos,vcfd->bcf->pos+1,pars->totSites);
#endif

    fprintf(stderr, "\n\t-> Finished reading sites\n");
}

void readSites_GT(vcfData *vcfd, paramStruct *pars, pairStruct **pairSt, blobStruct *blob) {
    const int nBlocks = blob->bootstraps->nBlocks;

    int contig_i = 0;
    int site_i = 0;
    int skip_site = 0;

    // read vcf records block by block
    // assume blocks are sorted in the same order as the vcf file
    for (int block_i = 0; block_i < nBlocks; block_i++) {
        int block_start = blob->blocks[block_i]->start;
        int block_end = blob->blocks[block_i]->end;

        int nSites_block_i = 0;

        // read vcf records until the end of the block (or the end of the vcf file)
        int ret = vcfd->records_next();
        while (1 == ret) {
            ret = vcfd->records_next();

            // if site is in block_i
            if (0 == strcmp(blob->blocks[block_i]->chr, bcf_seqname(vcfd->hdr, vcfd->bcf)) && block_start <= vcfd->bcf->pos && block_end > vcfd->bcf->pos) {
                nSites_block_i++;

            } else {
                // dragon
                //  if site is not in the block and the block is empty
                if (nSites_block_i == 0) {
                    ERROR("Block %d is empty", block_i);
                }

                // block change
                block_i++;
                block_start = blob->blocks[block_i]->start;
                block_end = blob->blocks[block_i]->end;
                nSites_block_i = 0;

                if (0 == strcmp(blob->blocks[block_i]->chr, bcf_seqname(vcfd->hdr, vcfd->bcf)) && block_start <= vcfd->bcf->pos && block_end > vcfd->bcf->pos) {
                    nSites_block_i++;
                } else {
                    // skip the site
                    continue;
                }
            }

            skip_site = get_JointGenoDist_GT(contig_i, site_i, vcfd, pars, blob, block_i);

            if (skip_site == 1) {
                fprintf(stderr, "\n->\tSkipping site %lu for all individuals\n\n", pars->totSites);
                pars->totSites++;

                continue;
            }

            ++site_i;
            pars->nSites++;
            pars->totSites++;
        }
    }
}

// dragon
void vcfData::init_JointGenoCountDistGL() {
    ASSERT(nJointClasses > 0);
    ASSERT(nIndCmb > 0);
    JointGenoCountDistGL = (double **)malloc(nIndCmb * sizeof(double *));
    for (int i = 0; i < nIndCmb; i++) {
        JointGenoCountDistGL[i] = (double *)malloc((nJointClasses + 1) * sizeof(double));
        for (int j = 0; j < nJointClasses + 1; j++) {
            JointGenoCountDistGL[i][j] = 0.0;
        }
    }
}

void vcfData::init_JointGenoProbDistGL() {
    ASSERT(nJointClasses > 0);
    ASSERT(nIndCmb > 0);
    JointGenoProbDistGL = (double **)malloc(nIndCmb * sizeof(double *));
    for (int i = 0; i < nIndCmb; i++) {
        JointGenoProbDistGL[i] = (double *)malloc((nJointClasses + 1) * sizeof(double));
        for (int j = 0; j < nJointClasses + 1; j++) {
            JointGenoProbDistGL[i][j] = 0.0;
        }
    }
}

void vcfData::init_JointGenoCountDistGT() {
    ASSERT(nIndCmb > 0);
    JointGenoCountDistGT = (int **)malloc(nIndCmb * sizeof(int *));
    for (int i = 0; i < nIndCmb; i++) {
        JointGenoCountDistGT[i] = (int *)malloc((nJointClasses + 1) * sizeof(int));
        for (int j = 0; j < nJointClasses + 1; j++) {
            JointGenoCountDistGT[i][j] = 0;
        }
    }
}

void vcfData::print_JointGenoCountDist() {
    if (outFiles->out_jgcd_fs != NULL) {
        kstring_t *kbuf = kbuf_init();

        if (args->doDist == 1) {
            for (int i = 0; i < nIndCmb; i++) {
                ksprintf(kbuf, "%i,", i);
                for (int j = 0; j < nJointClasses + 1; j++) {
                    ksprintf(kbuf, "%f", JointGenoCountDistGL[i][j]);
                    if (j == nJointClasses) {
                        ksprintf(kbuf, "\n");
                    } else {
                        ksprintf(kbuf, ",");
                    }
                }
            }
        } else if (args->doDist == 2) {
            for (int i = 0; i < nIndCmb; i++) {
                ksprintf(kbuf, "%i,", i);
                for (int j = 0; j < nJointClasses + 1; j++) {
                    ksprintf(kbuf, "%i", JointGenoCountDistGT[i][j]);
                    if (j == nJointClasses) {
                        ksprintf(kbuf, "\n");
                    } else {
                        ksprintf(kbuf, ",");
                    }
                }
            }
        }
        outFiles->out_jgcd_fs->write(kbuf);
        kbuf_destroy(kbuf);
    }
}

void vcfData::print_JointGenoProbDist() {
    if (args->printJointGenoProbDist != 0) {
        kstring_t *kbuf = kbuf_init();
        if (args->doDist == 1) {
            for (int i = 0; i < nIndCmb; i++) {
                ksprintf(kbuf, "%i,", i);
                for (int j = 0; j < nJointClasses + 1; j++) {
                    ksprintf(kbuf, "%f", JointGenoProbDistGL[i][j]);
                    if (j == nJointClasses) {
                        ksprintf(kbuf, "\n");
                    } else {
                        ksprintf(kbuf, ",");
                    }
                }
            }
        } else if (args->doDist == 2) {
            ASSERT(0 == 1);
        }
        ASSERT(kbuf != NULL);
        outFiles->out_jgpd_fs->write(kbuf);
        kbuf_destroy(kbuf);
    }
}

void vcfData::lngl_init(int doEM) {
    lngl = (double **)malloc(_lngl * sizeof(double *));

    ASSERT(nInd > 0);
    ASSERT(nGT > 0);
    ASSERT(_lngl > 0);

    for (size_t i = 0; i < _lngl; i++) {
        lngl[i] = (double *)malloc(nInd * nGT * sizeof(double));
        for (int indi = 0; indi < nInd; indi++) {
            const int lngls_ind_start = indi * nGT;
            for (int j = 0; j < nGT; j++) {
                lngl[i][lngls_ind_start + j] = NEG_INF;
            }
        }
    }
}

void vcfData::lngl_expand() {
    ASSERT(nInd > 0);
    ASSERT(nGT > 0);
    ASSERT(lngl != NULL);

    int prev_lngl = _lngl;
    _lngl = _lngl * 2;
    lngl = (double **)realloc(lngl, _lngl * sizeof(double *));
    ASSERT(lngl != NULL);
    for (int i = prev_lngl; i < (int)_lngl; i++) {
        ASSERT(lngl[i] == NULL);
        lngl[i] = (double *)malloc(nInd * nGT * sizeof(double));
        for (int indi = 0; indi < nInd; indi++) {
            const int lngls_ind_start = indi * nGT;
            for (int j = 0; j < nGT; j++) {
                lngl[i][lngls_ind_start + j] = NEG_INF;
            }
        }
    }
}

void vcfData::_print(FILE *fp) {
    fprintf(stderr, "\nNumber of samples: %i", nInd);
    fprintf(stderr, "\nNumber of contigs: %d", nContigs);
}

void vcfData::_print() {
    _print(stderr);
}

int vcfData::records_next() {
    int ret = -42;

    // at least BCF_UN_STR unpack is needed to use bcf->d.allele with bcf_get_genotypes

    if (NULL == itr) {
        // no region specified
        ret = bcf_read(in_fp, hdr, bcf);
        bcf_unpack(bcf, BCF_UN_ALL);
    } else {
        // region reading using iterator
        ret = bcf_itr_next(in_fp, itr, bcf);

        // TODO is this check necessary? htslib might be doing it already
        // if ( -1 > ret ){
        // 	// error
        // 	fprintf(stderr, "\n[ERROR:%d]\tError reading from VCF file",ret);
        // 	exit(1);
        // }

        bcf_unpack(bcf, BCF_UN_ALL);
    }

    if (-1 == ret) {
        // no more data
        return 0;
    }

    if (bcf->rlen != 1) {
        ERROR("VCF file has a REF allele of length %ld. This is currently not allowed.", bcf->rlen);
    }

    return 1;
}
