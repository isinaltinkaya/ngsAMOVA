#include "vcfReader.h"

#include "bootstrap.h"

bool gtData::pass_minInd_threshold(const int nInd) {
    for (int i = 0; i < nInd; ++i) {
        int *p = this->ind_ptr(i);

        if (1 == bcf_gt_is_missing(p[0]) || 1 == bcf_gt_is_missing(p[1])) {
            // minInd 0
            // only use sites shared across all individuals
            // so skip site when you first encounter nan
            if (0 == args->minInd) {
                return (false);
            }

            this->n_missing_ind++;
            // all individuals are missing
            if (nInd == this->n_missing_ind) {
                return (false);
            }

            // threshold requires: #non-missing inds >= minInd
            if ((nInd - this->n_missing_ind) < args->minInd) {
                return (false);
            }
        }
    }
    return (true);
}

int gtData::get_alleleidx2state(const int alleleidx) {
    if (alleleidx > 3) {
        ERROR("Allele not found in the alleleidx2state list");
    } else if (alleleidx < 0) {
        ERROR("Allele index is negative");
    }
    int state = alleleidx2state[alleleidx];
    if (state == -1) {
        ERROR("Allele not found in the alleleidx2state list");
    }
    return (state);
}

int gtData::get_n_derived_alleles_ind(const int ind_i) {
    int *p = this->ind_ptr(ind_i);

    if (1 == bcf_gt_is_missing(p[0]) || 1 == bcf_gt_is_missing(p[1]))
        return (-1);

    int ind_a1 = bcf_gt_allele(p[0]);
    int ind_a2 = bcf_gt_allele(p[1]);

    // homozygous for ancestral/major allele
    if (ind_a1 == ind_a2) {
        if (ind_a1 == 0) {
            return (0);
        } else if (ind_a1 == 1) {
            return (2);
        } else {
            NEVER;
        }
    }

    int s1 = get_alleleidx2state(ind_a1);
    int s2 = get_alleleidx2state(ind_a2);
    ASSERT(-1 != s1 && -1 != s2);

    return (s1 + s2);
}

// TODO .vcf vcf.gz .bcf .bcf.gz .*.csi .*.tbi
int require_index(paramStruct *pars) {
    if (pars->in_ft & IN_VCF) {
        if (NULL != args->in_region) {
            return (IDX_CSI);
        }
    }
    if (pars->in_ft & IN_DM) {
        return (IDX_NONE);
    }
    return (IDX_NONE);
}

int require_unpack() {
    // use GT
    if (2 == args->doDist) {
        // BCF_UN_STR unpack is needed to use bcf->d.allele with bcf_get_genotypes
        return (BCF_UN_STR);
    }

    return (0);  // no unpacking needed
}

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

void vcfData::site_gts_get(const int a1, const int a2) {
    if (a1 == a2) {
        ERROR("Ancestral/Major allele cannot be the same as the derived/minor allele.");
    }

    if (NULL != this->gts) {
        this->gts->site_gts_clear();
    } else {
        this->gts = new gtData;
    }
    gts->n_values = bcf_get_genotypes(hdr, bcf, &gts->data, &gts->size_e);

    if (gts->n_values <= 0) {
        ERROR("GT not present.");
    }

    gts->ploidy = gts->n_values / nInd;
    if (gts->ploidy != 2) {
        ERROR("Ploidy %d not supported.", gts->ploidy);
    }

    this->gts->intbase2state[a1] = 0;  // major/ancestral
    this->gts->intbase2state[a2] = 1;  // minor/derived

    for (int a = 0; a < bcf->n_allele; ++a) {
        // this->gts->acgt2alleles[acgt_charToInt[(int)*bcf->d.allele[a]]] = a;
        // this->gts->alleleidx2intbase[a] = acgt_charToInt[(int)*bcf->d.allele[a]];
        this->gts->alleleidx2state[a] = gts->intbase2state[acgt_charToInt[(int)*bcf->d.allele[a]]];
    }
    // TODO allow multiallelic/or has multiple alleles in alt column but actually all individuals are biallelic?
    // TODO should we only skip individuals with genotype alleles not in {0,1}?
    //  if (2 < bcf->n_allele) {
    //      NEVER;
    //  }
}

void gtData::site_gts_clear() {
    FREE(data);
    size_e = 0;
    n_values = 0;
    n_missing_ind = 0;
    ploidy = 0;

    intbase2state[0] = -1;
    intbase2state[1] = -1;
    intbase2state[2] = -1;
    intbase2state[3] = -1;

    alleleidx2state[0] = -1;
    alleleidx2state[1] = -1;
    alleleidx2state[2] = -1;
    alleleidx2state[3] = -1;
}

void vcfData::site_get_gls(void) {
    vcfData *v = this;
    this->gls = new glData(v);
}

gtData::gtData() {
    size_e = 0;

    intbase2state = new int[4];
    intbase2state[0] = -1;
    intbase2state[1] = -1;
    intbase2state[2] = -1;
    intbase2state[3] = -1;

    alleleidx2state = new int[4];
    alleleidx2state[0] = -1;
    alleleidx2state[1] = -1;
    alleleidx2state[2] = -1;
    alleleidx2state[3] = -1;
}

gtData::~gtData() {
    DEL1D(intbase2state);
    DEL1D(alleleidx2state);

    FREE(data);
}

int *gtData::ind_ptr(const int ind_i) {
    return (data + ind_i * ploidy);
}

extern const int acgt_charToInt[256] = {
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
    return bcf_alleles2gt(acgt_charToInt[(int)(unsigned char)a1], acgt_charToInt[(int)(unsigned char)a2]);
}

int bcf_alleles_get_gtidx(unsigned char a1, unsigned char a2) {
    return bcf_alleles2gt(acgt_charToInt[(int)a1], acgt_charToInt[(int)a2]);
}

// return 1: skip site for all individuals
int site_read_GL(const int contig_i, const int site_i, vcfData *vcfd, paramStruct *pars, pairStruct **pairs) {
    const int nInd = pars->nInd;

    glData lgls(vcfd);

    int cmbArr[pars->nIndCmb];
    for (int i = 0; i < pars->nIndCmb; i++) {
        cmbArr[i] = 0;
    }

    if (1 == bcf_is_snp(vcfd->bcf)) {
        int a1 = -1;
        int a2 = -1;

        if (NULL != pars->ancder) {
            a1 = acgt_charToInt[(int)pars->ancder->a1[contig_i][site_i]];
            a2 = acgt_charToInt[(int)pars->ancder->a2[contig_i][site_i]];

            if (pars->ancder->nContigs < contig_i + 1) {
                ERROR("Could not find contig %s (index: %d) in ancderFile.", bcf_seqname_safe(vcfd->hdr, vcfd->bcf), contig_i);
            }
            if (strcmp(pars->ancder->contigNames[contig_i], bcf_seqname_safe(vcfd->hdr, vcfd->bcf)) != 0) {
                ERROR("Contig name mismatch: %s != %s", pars->ancder->contigNames[contig_i], bcf_seqname_safe(vcfd->hdr, vcfd->bcf));
            }

        } else if (NULL != pars->majmin) {
            a1 = acgt_charToInt[(int)pars->majmin->a1[contig_i][site_i]];
            a2 = acgt_charToInt[(int)pars->majmin->a2[contig_i][site_i]];

            if (pars->majmin->nContigs < contig_i + 1) {
                ERROR("Could not find contig %s (index: %d) in majminFile.", bcf_seqname_safe(vcfd->hdr, vcfd->bcf), contig_i);
            }
            if (strcmp(pars->majmin->contigNames[contig_i], bcf_seqname_safe(vcfd->hdr, vcfd->bcf)) != 0) {
                ERROR("Contig name mismatch: %s != %s", pars->majmin->contigNames[contig_i], bcf_seqname_safe(vcfd->hdr, vcfd->bcf));  // dragon
            }

        } else if (1 == args->isSim) {
            // reference allele == 0 in d.allele
            ASSERT(1 == strlen(vcfd->bcf->d.allele[0]));
            a1 = acgt_charToInt[(int)*vcfd->bcf->d.allele[0]];
            // first alternative allele == 1 in d.allele
            ASSERT(1 == strlen(vcfd->bcf->d.allele[1]));
            a2 = acgt_charToInt[(int)*vcfd->bcf->d.allele[1]];

        } else {
            NEVER;
        }

        ASSERT(-1 != a1);
        ASSERT(-1 != a2);

        // debug
        ASSERT(a1 != a2);

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
            }
        }
    } else {
        ERROR("Only SNPs are supported at the moment");
        // TODO skip non-SNPs
    }

    return 0;
}

int get_JointGenoDist_GT(const int contig_i, const int site_i, vcfData *vcfd, paramStruct *pars, const int block_i) {
    const int nInd = pars->nInd;

    if (1 == bcf_is_snp(vcfd->bcf)) {
        int a1 = -1;
        int a2 = -1;

        if (NULL != pars->ancder) {
            a1 = acgt_charToInt[(int)pars->ancder->a1[contig_i][site_i]];
            a2 = acgt_charToInt[(int)pars->ancder->a2[contig_i][site_i]];

            if (pars->ancder->nContigs < contig_i + 1) {
                ERROR("Could not find contig %s (index: %d) in ancderFile.", pars->ancder->contigNames[contig_i], contig_i);
            }
            if (strcmp(pars->ancder->contigNames[contig_i], bcf_seqname_safe(vcfd->hdr, vcfd->bcf)) != 0) {
                ERROR("Contig name mismatch: %s != %s", pars->ancder->contigNames[contig_i], bcf_seqname_safe(vcfd->hdr, vcfd->bcf));
            }

        } else if (NULL != pars->majmin) {
            a1 = acgt_charToInt[(int)pars->majmin->a1[contig_i][site_i]];
            a2 = acgt_charToInt[(int)pars->majmin->a2[contig_i][site_i]];

            if (pars->majmin->nContigs < contig_i + 1) {
                ERROR("Could not find contig %s (index: %d) in majminFile.", pars->majmin->contigNames[contig_i], contig_i);  // dragon
            }
            if (strcmp(pars->majmin->contigNames[contig_i], bcf_seqname_safe(vcfd->hdr, vcfd->bcf)) != 0) {
                ERROR("Contig name mismatch: %s != %s", pars->majmin->contigNames[contig_i], bcf_seqname_safe(vcfd->hdr, vcfd->bcf));  // dragon
            }

        } else if (1 == args->isSim) {
            // reference allele == 0 in d.allele
            ASSERT(1 == strlen(vcfd->bcf->d.allele[0]));
            a1 = acgt_charToInt[(int)*vcfd->bcf->d.allele[0]];
            // first alternative allele == 1 in d.allele
            ASSERT(1 == strlen(vcfd->bcf->d.allele[1]));
            a2 = acgt_charToInt[(int)*vcfd->bcf->d.allele[1]];

        } else {
            NEVER;
        }

        ASSERT(-1 != a1);
        ASSERT(-1 != a2);

        vcfd->site_gts_get(a1, a2);

        if (!vcfd->gts->pass_minInd_threshold(nInd)) {
            return 1;
        }

        for (int i1 = 0; i1 < nInd; ++i1) {
            int i1_nder = vcfd->gts->get_n_derived_alleles_ind(i1);
            if (-1 == i1_nder) {
                continue;
            }

            for (int i2 = i1 + 1; i2 < nInd; i2++) {
                int i2_nder = vcfd->gts->get_n_derived_alleles_ind(i2);
                if (-1 == i2_nder) {
                    continue;
                }

                int pair_idx = nCk_idx(nInd, i1, i2);

                // no block bootstrapping
                if (-1 == block_i) {
                    vcfd->JointGenoCountDistGT[pair_idx][nDerToM33Idx[i1_nder][i2_nder]]++;
                    vcfd->JointGenoCountDistGT[pair_idx][9]++;  // #shared sites
                } else {
                    vcfd->jgcd_gt[block_i][pair_idx][nDerToM33Idx[i1_nder][i2_nder]]++;
                    vcfd->pair_shared_nSites[block_i][pair_idx]++;
                }
            }
        }
    } else {
        LOG("Skipping non-SNP position at site %d in contig %d", site_i, contig_i);
        return 1;
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

    if (require_index(pars) & IDX_CSI) {
        vcfd->idx = IO::load_bcf_csi_idx(args->in_vcf_fn);
    }

    // if (require_index(pars) & IDX_TBI) {
    // vcfd->tbx = IO::load_vcf_tabix_idx(args->in_vcf_fn);
    // vcfd->itr = tbx_itr_querys(vcfd->tbx, args->in_region);
    // }

    if (NULL != args->in_region) {
        vcfd->itr = bcf_itr_querys(vcfd->idx, vcfd->hdr, args->in_region);
        if (NULL == vcfd->itr) {
            ERROR("Could not parse region: %s. Please make sure region exists and defined in the correct format.", args->in_region);
        }
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

    if (NULL != v->gts) {
        DEL(v->gts);
    }

    if (NULL != v->gls) {
        DEL(v->gls);
    }

    DEL(v);
}

void readSites_GL(vcfData *vcfd, paramStruct *pars, pairStruct **pairSt, blobStruct *blob) {
    NEVER;
    // int skip_site = 0;
    // int contig_i = 0;
    // int site_i = 0;
    // char prev_contig[100];

    // while (vcfd->records_next()) {
    //     while (pars->nSites >= vcfd->_lngl) {
    //         vcfd->lngl_expand();
    //     }

    //     // if cleared in the previous loop to skip the site, allocate memory again
    //     if (NULL == vcfd->lngl[pars->nSites]) {
    //         vcfd->lngl[pars->nSites] = (double *)malloc(pars->nInd * vcfd->nGT * sizeof(double));
    //         for (int indi = 0; indi < pars->nInd; indi++) {
    //             const int lngls_ind_start = indi * vcfd->nGT;
    //             for (int j = 0; j < vcfd->nGT; j++) {
    //                 vcfd->lngl[pars->nSites][lngls_ind_start + j] = NEG_INF;
    //             }
    //         }
    //     }

    //     const char *contig_id = bcf_seqname(vcfd->hdr, vcfd->bcf);
    //     if (0 != contig_i && 0 != strcmp(contig_id, prev_contig)) {
    //         strcpy(prev_contig, contig_id);
    //         ++contig_i;
    //         if (contig_i >= vcfd->nContigs) {
    //             ERROR("Number of contigs in the VCF file header is larger than the number of contigs in the VCF file");
    //         }
    //         site_i = 0;
    //     }

    //     skip_site = site_read_GL(contig_i, site_i, vcfd, pars, pairSt);

    //     if (skip_site == 1) {
    //         fprintf(stderr, "\n->\tSkipping site %lu for all individuals\n\n", pars->totSites);
    //         pars->totSites++;

    //         //  next loop will skip this and use the same pars->nSites
    //         FREE(vcfd->lngl[pars->nSites]);
    //         continue;
    //     }

    //     // check if contig changed
    //     pars->nSites++;
    //     ++site_i;
    //     pars->totSites++;
    // }

    // vcfd->nSites = pars->nSites;
    // vcfd->totSites = pars->totSites;

    // fprintf(stderr, "\n\t-> Finished reading sites\n");
}

void readSites_GL(vcfData *vcfd, paramStruct *pars, pairStruct **pairSt) {
    int skip_site = 0;
    int contig_i = -1;
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
        if (0 != strcmp(contig_id, prev_contig)) {
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
        skip_site = get_JointGenoDist_GT(contig_i, site_i, vcfd, pars, -1);

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
        // int ret = vcfd->records_next();
        while (1 == (vcfd->records_next())) {
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

            skip_site = get_JointGenoDist_GT(contig_i, site_i, vcfd, pars, block_i);

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
        outFiles->out_jgcd_fs->kbuf = kbuf_init();
        kstring_t *kbuf = outFiles->out_jgcd_fs->kbuf;
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
        outFiles->out_jgcd_fs->kbuf_write();
    }
}

void vcfData::print_JointGenoProbDist() {
    if (args->printJointGenoProbDist != 0) {
        outFiles->out_jgpd_fs->kbuf = kbuf_init();
        kstring_t *kbuf = outFiles->out_jgpd_fs->kbuf;
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
        outFiles->out_jgpd_fs->kbuf_write();
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

void vcfData::unpack() {
    int unpack = require_unpack();
    if (0 != unpack) {
        bcf_unpack(bcf, unpack);
    }
}

int vcfData::records_next() {
    int ret = -42;

    if (NULL == itr) {
        // no region specified
        ret = bcf_read(in_fp, hdr, bcf);
    } else {
        // region reading using iterator
        ret = bcf_itr_next(in_fp, itr, bcf);

        if (-1 > ret) {
            ERROR("An error occurred while reading the VCF file.");
        }
    }

    unpack();

    ASSERT(0 == bcf->errcode);

    if (-1 == ret) {
        // no more data
        return 0;
    }

    return 1;
}
