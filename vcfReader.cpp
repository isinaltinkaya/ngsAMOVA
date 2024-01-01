#include "vcfReader.h"

#include "bootstrap.h"

void change_contig_reset_nSites() {

}

// NB> in full genome arrays use pars->nSites (for total)
// in contig arrays use site_i (per contig)
void arrays_nSites_expand(paramStruct* pars, lnglStruct* lngl, int** a1a2) {
    const int oldSize = pars->nSites_arrays_size;
    ASSERT(pars->nSites_arrays_size == lngl->size1);
    const int newSize = oldSize + 4096;
    pars->nSites_arrays_size = newSize;
    lngl->size1 = newSize;

    if (lngl != NULL) {
        lngl->size1 = newSize;
        REALLOC(lngl->d, lngl->size1, double**);

        for (size_t i = oldSize; i < lngl->size1;++i) {
            lngl->d[i] = (double*)malloc(lngl->size2 * sizeof(double));
            ASSERT(lngl->d[i] != NULL);
            for (size_t j = 0;j < lngl->size2;++j) {
                lngl->d[i][j] = NEG_INF;
            }
        }
    }

    if (a1a2 != NULL) {
        REALLOC(pars->a1a2, pars->nSites_arrays_size, int**);

        for (size_t i = oldSize; i < newSize;++i) {
            pars->a1a2[i] = (int*)malloc(2 * sizeof(int));
            pars->a1a2[i][0] = -1;
            pars->a1a2[i][1] = -1;
        }
    }

}


const char* vcfData::get_contig_name(void) {
    const char* ret = bcf_hdr_id2name(this->hdr, this->rec ? this->rec->rid : -1);
    if (NULL == ret) {
        ERROR("Could not retrieve contig name for record %d.", this->rec->rid);
    }
    return (ret);
}

const char* vcfData::get_contig_name(const int32_t i) {
    const char* ret = bcf_hdr_id2name(this->hdr, this->rec ? i : -1);
    if (NULL == ret) {
        ERROR("Could not retrieve contig name for record %d.", i);
    }
    return (ret);
}

bool gtData::pass_minInd_threshold(const int nInd) {
    for (int i = 0; i < nInd; ++i) {
        int* p = this->ind_ptr(i);

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

    return (state);
}

// -1   missing
// -2
int gtData::get_n_derived_alleles_ind(const int ind_i) {
    int* p = this->ind_ptr(ind_i);

    if (1 == bcf_gt_is_missing(p[0]) || 1 == bcf_gt_is_missing(p[1]))
        return (-1);

    // 0, 1, 2, 3 (e.g. 0|0)
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
    if (-1 == s1 || -1 == s2) {
        return (-2);
    }

    return (s1 + s2);
}

int require_unpack() {

    int val = 0;

    // BCF_UN_STR: if we need to fetch alleles from bcf->d.allele 
    if (2 == args->doDist) {
        if (val < BCF_UN_STR) {
            val = BCF_UN_STR;
        }
    } else if (1 == args->doDist) {
        if (val < BCF_UN_STR) {
            val = BCF_UN_STR;
        }
    }

    // TODO set to all for now
    val = BCF_UN_ALL;
    return(val);
}

// TODO .vcf vcf.gz .bcf .bcf.gz .*.csi .*.tbi
int require_index(paramStruct* pars) {
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



glData::glData(vcfData* vcfd) {
    size_e = 0;
    n_values = bcf_get_format_float(vcfd->hdr, vcfd->rec, "GL", &data, &size_e);
    if (n_values <= 0) {
        ERROR("Could not read GL tag from the VCF file.");
    }

    n_gls = n_values / vcfd->nInd;

    if (-1 == vcfd->allele_unobserved) {
        if (10 == n_gls) {
            return;
        } else if (3 == n_gls) {
            return;
        } else {
            ERROR("Found %d values (%d GLs) in the GL field at %s:%ld.", n_values, n_gls, vcfd->get_contig_name(), vcfd->rec->pos + 1);
        }
    } else {
        // e.g. A, UNOBSERVED -> 3 GLs
        if (3 == n_gls) {
            NEVER; // site shouldve been already skipped
        }
        // e.g. A,C, UNOBSERVED -> 6 GLs
        if (6 == n_gls) {
            return;
        }
        // e.g. A,C,G, UNOBSERVED -> 10 GLs
        if (10 == n_gls) {
            return;
        }
        // e.g. A,C,G,T, UNOBSERVED -> 15 GLs
        if (15 == n_gls) {
            return;
        }
        ERROR("Found %d values in the GL field at %s:%ld.", n_gls, vcfd->get_contig_name(), vcfd->rec->pos + 1);
    }
}


glData::~glData() {
    FREE(data);
}


//TODO deprec
float* glData::ind_ptr(const int ind_i) {
    return (data + ind_i * n_gls);
}


void vcfData::site_gts_get(const int a1, const int a2) {
    // if (a1 == a2) {  // debug
    //     ERROR("Ancestral/Major allele cannot be the same as the derived/minor allele.");
    // }

    // if (NULL != this->gts) {
    //     this->gts->site_gts_clear();
    // } else {
    //     this->gts = new gtData;
    // }
    // gts->n_values = bcf_get_genotypes(hdr, rec, &gts->data, &gts->size_e);

    // if (gts->n_values <= 0) {
    //     ERROR("GT not present.");
    // }

    // gts->ploidy = gts->n_values / nInd;
    // if (gts->ploidy != 2) {
    //     ERROR("Ploidy %d not supported.", gts->ploidy);
    // }

    // // intbase = int representation of internal acgt base
    // this->gts->intbase2state[a1] = 0;  // major/ancestral
    // this->gts->intbase2state[a2] = 1;  // minor/derived

    // for (int a = 0; a < rec->n_allele; ++a) {
    //     // this->gts->acgt2alleles[acgt_charToInt[(int)*bcf->d.allele[a]]] = a;
    //     // this->gts->alleleidx2intbase[a] = acgt_charToInt[(int)*bcf->d.allele[a]];
    //     this->gts->alleleidx2state[a] = gts->intbase2state[acgt_charToInt[(int)*rec->d.allele[a]]];
    //     //TODO fix
    // }
}


void gtData::site_gts_clear() {
    FREE(data);
    size_e = 0;
    n_values = 0;
    n_missing_ind = 0;

    intbase2state[0] = -1;
    intbase2state[1] = -1;
    intbase2state[2] = -1;
    intbase2state[3] = -1;

    alleleidx2state[0] = -1;
    alleleidx2state[1] = -1;
    alleleidx2state[2] = -1;
    alleleidx2state[3] = -1;
}


gtData::gtData(vcfData* vcfd) {
    size_e = 0;

    n_values = bcf_get_genotypes(vcfd->hdr, vcfd->rec, &data, &size_e);
    if (n_values <= 0) {
        ERROR("Could not read GT tag from the VCF file.");
    }

    if ((n_values / vcfd->nInd) != PROGRAM_PLOIDY) {
        ERROR("Ploidy %d not supported.", n_values / vcfd->nInd);
    }

    //TODO DEPREC
    // intbase2state = new int[4];
    // intbase2state[0] = -1;
    // intbase2state[1] = -1;
    // intbase2state[2] = -1;
    // intbase2state[3] = -1;

    //TODO DEPREC
    // alleleidx2state = new int[4];
    // alleleidx2state[0] = -1;
    // alleleidx2state[1] = -1;
    // alleleidx2state[2] = -1;
    // alleleidx2state[3] = -1;
}

gtData::~gtData() {
    // delete(intbase2state);
    // delete(alleleidx2state);

    FREE(data);
}

int* gtData::ind_ptr(const int ind_i) {
    return (data + ind_i * PROGRAM_PLOIDY);
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


// return 1: skip site for all individuals
int read_site_with_alleles(const int site_i, vcfData* vcfd, paramStruct* pars) {

    const int nAlleles = vcfd->rec->n_allele;

    vcfd->allele_unobserved = -1;
    for (int i = 0; i < nAlleles; ++i) {

        // if allele is single char
        if ('\0' == vcfd->rec->d.allele[i][1]) {
            if ('*' == vcfd->rec->d.allele[i][0]) {
                // '*' symbol = allele missing due to overlapping deletion
                // source: VCF specs v4.3 
                IO::vprint(2, "Skipping site at %s:%ld. Reason: REF allele is set to '*' (allele missing due to overlapping deletion).\n", vcfd->get_contig_name(), vcfd->rec->pos + 1);
                return(1);
            }
            if ('.' == vcfd->rec->d.allele[i][0]) {
                // '.' symbol = MISSING value (no variant)
                // source: VCF specs v4.3 
                IO::vprint(2, "Skipping site at %s:%ld. Reason: REF allele is set to '*' (allele missing due to overlapping deletion).\n", vcfd->get_contig_name(), vcfd->rec->pos + 1);
                return(1);
            }

            continue;
        }

        // -> allele is not a single char

        // if not a single char, it can be one of these:
        // 1) An unobserved allele notation (<*> or <NON_REF>)
        // 2) Insertion/deletion 
        // otherwise we do not support such a site
        if (vcfd->rec->d.allele[i][0] != '<') {

            if (0 == i) {
                // if at allele REF and all chars are one of ACGT, it is indel
                int j = 0;
                while (vcfd->rec->d.allele[i][j] != '\0') {
                    if (vcfd->rec->d.allele[i][j] != 'A' && vcfd->rec->d.allele[i][j] != 'C' && vcfd->rec->d.allele[i][j] != 'G' && vcfd->rec->d.allele[i][j] != 'T') {
                        ERROR("Found bad allele '%s' at %s:%ld. If you believe this is a valid allele, please contact the developers.\n", vcfd->rec->d.allele[i], vcfd->get_contig_name(), vcfd->rec->pos + 1);
                        break;
                    }
                    j++;
                }
                IO::vprint(2, "Skipping site at %s:%ld. Reason: REF allele is an INDEL.\n", vcfd->get_contig_name(), vcfd->rec->pos + 1);
                return(1);
            }

            ERROR("Found bad allele '%s' at %s:%ld. If you believe this is a valid allele, please contact the developers.\n", vcfd->rec->d.allele[i], vcfd->get_contig_name(), vcfd->rec->pos + 1);
        }

        // -> allele starts with '<'

        // if alelle is 3 chars long
        if (0 == vcfd->rec->d.allele[i][3]) {
            // if unobserved allele <*>
            if ('*' == vcfd->rec->d.allele[i][1] && '>' == vcfd->rec->d.allele[i][2]) {
                vcfd->allele_unobserved = i;
                continue;
            } else {
                ERROR("Found bad allele '%s' at %s:%ld. If you believe this is a valid allele, please contact the developers.\n", vcfd->rec->d.allele[i], vcfd->get_contig_name(), vcfd->rec->pos + 1);
            }


            // if allele is 9 chars long
        } else if (0 == vcfd->rec->d.allele[i][9]) {
            if ('N' == vcfd->rec->d.allele[i][1] && 'O' == vcfd->rec->d.allele[i][2] && 'N' == vcfd->rec->d.allele[i][3] && '_' == vcfd->rec->d.allele[i][4] && 'R' == vcfd->rec->d.allele[i][5] && 'E' == vcfd->rec->d.allele[i][6] && 'F' == vcfd->rec->d.allele[i][7] && '>' == vcfd->rec->d.allele[i][8]) {
                // if unobserved allele <NON_REF>
                vcfd->allele_unobserved = i;
                continue;
            } else {
                ERROR("Found bad allele '%s' at %s:%ld. If you believe this is a valid allele, please contact the developers.\n", vcfd->rec->d.allele[i], vcfd->get_contig_name(), vcfd->rec->pos + 1);
            }
        }
        NEVER;
    }

    int b1 = -1;
    int b2 = -1;

    ASSERT(nAlleles != 0);

    if (1 == nAlleles) {
        IO::vprint(2, "Skipping site at %s:%ld. Reason: No data found at site.\n", vcfd->get_contig_name(), vcfd->rec->pos + 1);
        return(1);
    }

    // if REF==ALLELE_UNOBSERVED
    if (0 == vcfd->allele_unobserved) {


        // NEVER; //TODO this is only supported if major minor file or refalt file or an estimation method for these is defined; add proper check and handling

        // check if ALT is set to any allele, maybe site has data but no ref allele specified

        // if 1==nAlleles ALT is not set to any allele, so no data at site
        // this is already skipped above

        if (2 == nAlleles) {
            // ALT is set but REF is not, so probably an invariable site with no REF allele specified

            if (args->rmInvarSites) {
                IO::vprint(2, "Skipping site at %s:%ld. Reason: Site is not variable.\n", vcfd->get_contig_name(), vcfd->rec->pos + 1);
                return(1);
            }
        }

        b1 = BASE_UNOBSERVED;
        b2 = acgt_charToInt[(int)vcfd->rec->d.allele[1][0]];


    } else if (1 == vcfd->allele_unobserved) {

        // NEVER;//TODO

        // if ALT==ALLELE_UNOBSERVED

        if (2 == nAlleles) {
            // invariable site
            if (args->rmInvarSites) {
                IO::vprint(2, "Skipping site at %s:%ld. Reason: Site is not variable.\n", vcfd->get_contig_name(), vcfd->rec->pos + 1);
                return(1);
            }
        }

        b1 = acgt_charToInt[(int)vcfd->rec->d.allele[0][0]];
        b2 = BASE_UNOBSERVED;

    } else if (-1 == vcfd->allele_unobserved) {
        // NEVER;//TODO
        // if no unobserved alt allele notation is found

        if (1 == nAlleles) {
            // Invariable
            if (args->rmInvarSites) {
                IO::vprint(2, "Skipping site at %s:%ld. Reason: Site is not variable.\n", vcfd->get_contig_name(), vcfd->rec->pos + 1);
                return(1);
            }
        } else {
            b2 = acgt_charToInt[(int)vcfd->rec->d.allele[1][0]];
        }

        b1 = acgt_charToInt[(int)vcfd->rec->d.allele[0][0]];

    } else {

        // unobserved allele notation is found but is not a1 or a2
        b1 = acgt_charToInt[(int)vcfd->rec->d.allele[0][0]];
        b2 = acgt_charToInt[(int)vcfd->rec->d.allele[1][0]];

    }


    if (args->rmInvarSites && args->doDist != 2) {
        int32_t* ads = NULL;
        int32_t n_ads = 0;
        int32_t n = 0;
        n = bcf_get_format_int32(vcfd->hdr, vcfd->rec, "AD", &ads, &n_ads);
        if (n <= 0) {
            ERROR("Could not read AD tag from the VCF file. AD tag is required for removing invariable sites via --rm-invar-sites 1.");
        }


        // make sure that AD is non-0 for at least 2 alleles 
        int n_non0[nAlleles];
        for (int i = 0; i < nAlleles; ++i) {
            n_non0[i] = 0;
        }
        // loop through samples 
        for (int i = 0; i < vcfd->nInd; ++i) {
            // loop through alleles
            for (int j = 0; j < nAlleles; ++j) {
                if (ads[i * nAlleles + j] > 0) {
                    n_non0[j]++;
                }
            }
        }
        int tot_non0 = 0;
        for (int i = 0; i < nAlleles; ++i) {
            if (n_non0[i] > 0) {
                tot_non0++;
            }
        }
        FREE(ads);
        if (tot_non0 < 2) {
            IO::vprint(2, "Skipping site at %s:%ld. Reason: Site is not variable.\n", vcfd->get_contig_name(), vcfd->rec->pos + 1);
            return(1);
        }



    }

    pars->a1a2[pars->nSites][0] = b1;
    pars->a1a2[pars->nSites][1] = b2;


    // DEVASSERT(b1 != -1);
    // DEVASSERT(b2 != -1);
    DEVASSERT(b1 != b2);

    //TODO
    // DEVPRINT("b1: %d, b2: %d\n", b1, b2);

    return(0);

}

// return 1: skip site for all individuals
int site_read_GL(const int contig_i, const int site_i, vcfData* vcfd, paramStruct* pars) {

    const int nAlleles = vcfd->rec->n_allele;

    int ret;
    if (0 != (ret = read_site_with_alleles(site_i, vcfd, pars))) {
        return(ret);
    }

    //TODO
    int b1 = pars->a1a2[pars->nSites][0];
    int b2 = pars->a1a2[pars->nSites][1];

    const int nInd = pars->nInd;

    glData lgls(vcfd);

    //TODO make sure these are compatible!!!
    // DEVPRINT("n_gls: %d nGT:%d", lgls.n_gls, vcfd->nGT);

    float* sample_lgls = NULL;
    double* sample_lngls = NULL;
    bool isMissing;
    bool includeInds[nInd];

    for (int indi = 0; indi < nInd; indi++) {

        includeInds[indi] = false;
        isMissing = false;

        sample_lgls = lgls.data + indi * lgls.n_gls;

        if (bcf_float_is_missing(sample_lgls[0])) {
            // missing data check 1
            isMissing = true;

        } else {

            // missing data check 2
            int z = 1;
            for (int i = 1; i < lgls.n_gls; ++i) {
                if (sample_lgls[0] == sample_lgls[i]) {
                    ++z;
                }
            }
            if (z == lgls.n_gls) {
                // missing data (all gl values are the same)
                isMissing = true;
            }
        }

        if (isMissing) {

            // if minInd 0 (the program will only use sites shared across all individuals) 
            // skip site when you encounter any missing data
            if (0 == args->minInd) {
                return(1);
            }

            // skip site if there are 2 inds and at least one of them is missing
            if (2 == nInd) {
                return(1);
            }

            lgls.n_missing_ind++;

            // all individuals are missing
            if (nInd == lgls.n_missing_ind) {
                return(1);
            }

            // skip site if minInd is defined and #non-missing inds=<nInd
            if (args->minInd != 2) {
                if ((nInd - lgls.n_missing_ind) < args->minInd) {
                    IO::vprint(2, "Skipping site at %s:%ld. Reason: Number of non-missing individuals is less than the minimum number of individuals -minInd.\n", vcfd->get_contig_name(), vcfd->rec->pos + 1);
                    return(1);
                }
            }
            continue;
        }
        includeInds[indi] = true;
    }

    // -> site passed all missingness filters



//TODO
    int a1a1 = bcf_alleles2gt(b1, b1);
    int a1a2 = bcf_alleles2gt(b1, b2);
    int a2a2 = bcf_alleles2gt(b2, b2);



    for (int indi = 0; indi < nInd; indi++) {

        if (includeInds[indi]) {

            //TODO currently this only works if a1a2 is defined using REF and the first ALT

            sample_lgls = lgls.data + indi * lgls.n_gls;

            sample_lngls = vcfd->lngl->d[pars->nSites] + indi * vcfd->nGT;

            // sample_lngls[0] = (double)LOG2LN(sample_lgls[0]);
            // sample_lngls[1] = (double)LOG2LN(sample_lgls[1]);
            // sample_lngls[2] = (double)LOG2LN(sample_lgls[2]);

    //TODO
            sample_lngls[0] = (double)LOG2LN(sample_lgls[a1a1]);
            sample_lngls[1] = (double)LOG2LN(sample_lgls[a1a2]);
            sample_lngls[2] = (double)LOG2LN(sample_lgls[a2a2]);

            // fprintf(stdout, "%f,%f,%f\n", sample_lngls[0], sample_lngls[1], sample_lngls[2]);


            for (int indi2 = indi + 1; indi2 < nInd; indi2++) {
                if (includeInds[indi2]) {
                    vcfd->snSites[nCk_idx(nInd, indi, indi2)]++;
                }
            }

        }
    }

    return(0);

}


int site_read_GT(const int contig_i, const int site_i, vcfData* vcfd, paramStruct* pars) {

    const int nAlleles = vcfd->rec->n_allele;

    int ret;
    if (0 != (ret = read_site_with_alleles(site_i, vcfd, pars))) {
        return(ret);
    }

    const int nInd = pars->nInd;

    bool isMissing;

    gtData gts(vcfd);


    int32_t* sample_gts = NULL;

    bool includeInds[nInd];

    for (int indi = 0; indi < nInd; ++indi) {

        includeInds[indi] = false;
        isMissing = false;

        sample_gts = gts.data + indi * PROGRAM_PLOIDY;

        if (1 == bcf_gt_is_missing(sample_gts[0]) || 1 == bcf_gt_is_missing(sample_gts[1])) {
            isMissing = true;
        }

        if (isMissing) {
            // minInd 0
            // only use sites shared across all individuals
            // so skip site when you first encounter nan
            if (0 == args->minInd) {
                return(1);
            }

            // skip site if there are 2 inds and at least one of them is missing
            if (2 == nInd) {
                return(1);
            }

            gts.n_missing_ind++;

            // all individuals are missing
            if (nInd == gts.n_missing_ind) {
                return(1);
            }

            // skip site if minInd is defined and #non-missing inds=<nInd
            if (args->minInd != 2) {
                if ((nInd - gts.n_missing_ind) < args->minInd) {
                    IO::vprint(2, "Skipping site at %s:%ld. Reason: Number of non-missing individuals is less than the minimum number of individuals -minInd.\n", vcfd->get_contig_name(), vcfd->rec->pos + 1);
                    return(1);
                }
            }
            continue;
        }
        includeInds[indi] = true;
    }

    // -> site passed all missingness filters


    int32_t* i1_gts = NULL;
    int32_t* i2_gts = NULL;


    int a1_base = pars->a1a2[pars->nSites][0];
    int a2_base = pars->a1a2[pars->nSites][1];

    int ind_a1_base, ind_a2_base;


    int recAlleles2a1a2[nAlleles];
    int recAlleles2Base[nAlleles];
    int nMatch = 0;
    for (int a = 0;a < nAlleles;++a) {
        if (a == vcfd->allele_unobserved) {
            recAlleles2Base[a] = BASE_UNOBSERVED;
            recAlleles2a1a2[a] = -1;
        } else {
            recAlleles2Base[a] = acgt_charToInt[(int)vcfd->rec->d.allele[a][0]];
            if (pars->a1a2[pars->nSites][0] == recAlleles2Base[a]) {
                recAlleles2a1a2[a] = 0;
                nMatch++;
            } else if (pars->a1a2[pars->nSites][1] == recAlleles2Base[a]) {
                recAlleles2a1a2[a] = 1;
                nMatch++;
            }
        }
    }
    if (nMatch != 2) {
        ERROR("At site %s:%ld, the alleles in the VCF record (%s) does not match the alleles defined based on user input (a1:%c,a2:%c).", vcfd->get_contig_name(), vcfd->rec->pos + 1, vcfd->rec->d.allele[0], "ACGT"[pars->a1a2[pars->nSites][0]], "ACGT"[pars->a1a2[pars->nSites][1]]);
    }

    int i1_nder = -1;
    int i2_nder = -1;
    int gta1 = -1;
    int gta2 = -1;
    int pair_idx = -1;

    if (args->rmInvarSites) {

        int tot_nder = 0;

        for (int i = 0; i < nInd; ++i) {

            if (includeInds[i]) {

                i1_nder = -1;
                gta1 = -1;
                gta2 = -1;

                i1_gts = gts.data + i * PROGRAM_PLOIDY;

                gta1 = bcf_gt_allele(i1_gts[0]);
                gta2 = bcf_gt_allele(i1_gts[1]);

                i1_nder = recAlleles2a1a2[gta1] + recAlleles2a1a2[gta2];
                tot_nder += i1_nder;

            }

        }


        if (tot_nder == 0) {
            IO::vprint(2, "Skipping site at %s:%ld. Reason: Site is not variable.\n", vcfd->get_contig_name(), vcfd->rec->pos + 1);
            return(1);
        } else if (tot_nder == nInd * PROGRAM_PLOIDY) {
            IO::vprint(2, "Skipping site at %s:%ld. Reason: Site is not variable.\n", vcfd->get_contig_name(), vcfd->rec->pos + 1);
            return(1);
        }
    }

    for (int i1 = 0; i1 < nInd - 1; ++i1) {

        if (includeInds[i1]) {

            i1_nder = -1;
            gta1 = -1;
            gta2 = -1;

            i1_gts = gts.data + i1 * PROGRAM_PLOIDY;

            gta1 = bcf_gt_allele(i1_gts[0]);
            gta2 = bcf_gt_allele(i1_gts[1]);
            // DEVPRINT("%d %d", gta1, gta2);
            // DEVPRINT("%d %d", recAlleles2Base[gta1], recAlleles2Base[gta2]);
            // DEVPRINT("%c %c", "ACGT"[recAlleles2Base[gta1]], "ACGT"[recAlleles2Base[gta2]]);
            // DEVPRINT("%d %d", recAlleles2a1a2[gta1], recAlleles2a1a2[gta2]);
            // DEVPRINT("%c %c", "Mm"[recAlleles2a1a2[gta1]], "Mm"[recAlleles2a1a2[gta2]]);


            DEVASSERT(recAlleles2a1a2[gta1] != -1);
            DEVASSERT(recAlleles2a1a2[gta2] != -1);

            i1_nder = recAlleles2a1a2[gta1] + recAlleles2a1a2[gta2];


            for (int i2 = i1 + 1; i2 < nInd; i2++) {


                if (includeInds[i2]) {

                    pair_idx = -1;
                    i2_nder = -1;
                    gta1 = -1;
                    gta2 = -1;

                    i2_gts = gts.data + i2 * PROGRAM_PLOIDY;

                    gta1 = bcf_gt_allele(i2_gts[0]);
                    gta2 = bcf_gt_allele(i2_gts[1]);

                    i2_nder = recAlleles2a1a2[gta1] + recAlleles2a1a2[gta2];


                    pair_idx = nCk_idx(nInd, i1, i2);

                    // no block bootstrapping
                    //TODO add block thingy here?
                    // if (-1 == block_i) {
                    //     vcfd->jointGenotypeMatrixGT[pair_idx][nDerToM33Idx[i1_nder][i2_nder]]++;
                    //     vcfd->snSites[pair_idx]++;// #shared sites
                    // } else {
                        // if no block bootstrapping, block_i is -1
                    // vcfd->jgcd_gt[0][pair_idx][nDerToM33Idx[i1_nder][i2_nder]]++;
                    // vcfd->pair_shared_nSites[0][pair_idx]++;
                    // }

                    vcfd->jointGenotypeMatrixGT[pair_idx][nDerToM33Idx[i1_nder][i2_nder]]++;
                    vcfd->snSites[pair_idx]++;


                }


            }
        }

    }


    return(0);
}



vcfData* vcfData_init(paramStruct* pars) {

    vcfData* vcfd = new vcfData;

    vcfd->DO_BCF_UNPACK = require_unpack();

    vcfd->in_fp = bcf_open(args->in_vcf_fn, "r");
    if (vcfd->in_fp == NULL) {
        ERROR("Could not open bcf file: %s\n", args->in_vcf_fn);
    }

    vcfd->hdr = bcf_hdr_read(vcfd->in_fp);
    vcfd->rec = bcf_init();

    if (require_index(pars) & IDX_CSI) {
        vcfd->idx = IO::load_bcf_csi_idx(args->in_vcf_fn);
    }

    if (require_index(pars) & IDX_TBI) {
        vcfd->tbx = IO::load_vcf_tabix_idx(args->in_vcf_fn);
        vcfd->itr = tbx_itr_querys(vcfd->tbx, args->in_region);
    }

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
    pars->nIndCmb = nC2(pars->nInd);

    vcfd->nInd = pars->nInd;
    vcfd->nIndCmb = pars->nIndCmb;

    pars->pidx2inds = (int**)malloc(pars->nIndCmb * sizeof(int*));
    ASSERT(pars->pidx2inds != NULL);
    for (int i = 0;i < pars->nIndCmb;++i) {
        pars->pidx2inds[i] = (int*)malloc(2 * sizeof(int));
    }

    int pidx = 0;

    //TODO check for loops
    for (int i1 = 0; i1 < pars->nInd - 1;++i1) {
        for (int i2 = i1 + 1; i2 < pars->nInd;++i2) {
            pidx = nCk_idx(pars->nInd, i1, i2);
            pars->pidx2inds[pidx][0] = i1;
            pars->pidx2inds[pidx][1] = i2;
        }
    }

    if (1 == args->doDist) {
        // EM using 3 GL values
        if (args->doEM == 1) {
            vcfd->nGT = 3;
            vcfd->nJointClasses = 9;
        }
        // EM using 10 GL values
        else if (args->doEM == 2) {
            NEVER;//TODO
            // vcfd->nGT = 10;
            // vcfd->nJointClasses = 100;
        }
    } else if (2 == args->doDist) {
        // GT
        vcfd->nGT = 3;
        vcfd->nJointClasses = 9;
    } else if (1 == args->doIbd) {
        // using 3 GL values
        vcfd->nGT = 3;
        vcfd->nJointClasses = 9;
    } else if (2 == args->doIbd) {
        // GT
        vcfd->nGT = 3;
        vcfd->nJointClasses = 9;
    } else {
        NEVER;
    }


    if (1 == args->doDist) {

        vcfd->lngl = new lnglStruct(vcfd->nGT, pars->nInd);

        // \def jointGenotypeMatrix[nIndCmb][nJointClasses]
        // jointGenotypeMatrix[i][j] == Value for pair i at joint class j
        vcfd->jointGenotypeMatrixGL = (double**)malloc(pars->nIndCmb * sizeof(double*));
        vcfd->snSites = (int*)malloc(pars->nIndCmb * sizeof(int));
        for (int i = 0; i < pars->nIndCmb; i++) {
            vcfd->snSites[i] = 0;

            vcfd->jointGenotypeMatrixGL[i] = (double*)malloc((vcfd->nJointClasses) * sizeof(double));
            for (int j = 0; j < vcfd->nJointClasses; j++) {
                // set initial guess: 1/9 flat prior
                vcfd->jointGenotypeMatrixGL[i][j] = (double)FRAC_1_9;
            }
        }

    } else if (2 == args->doDist) {

        // \def jointGenotypeMatrix[nIndCmb][nJointClasses]
        // jointGenotypeMatrix[i][j] == Value for pair j at joint class i
        vcfd->jointGenotypeMatrixGT = (int**)malloc(pars->nIndCmb * sizeof(int*));
        vcfd->snSites = (int*)malloc(pars->nIndCmb * sizeof(int));
        for (int i = 0; i < pars->nIndCmb; i++) {
            vcfd->snSites[i] = 0;
            vcfd->jointGenotypeMatrixGT[i] = (int*)malloc((vcfd->nJointClasses) * sizeof(int));
            for (int j = 0; j < vcfd->nJointClasses; j++) {
                vcfd->jointGenotypeMatrixGT[i][j] = 0;
            }
        }

    }
    //////

    BEGIN_LOGSECTION;
    LOG("Found %d individuals in the VCF file", pars->nInd);
    LOG("Number of individual pairs: %d", pars->nIndCmb);
    END_LOGSECTION;


    pars->indNames = (char**)malloc(pars->nInd * sizeof(char*));
    for (int i = 0;i < pars->nInd;++i) {
        pars->indNames[i] = (char*)strdup(vcfd->hdr->samples[i]);

    }

    vcfd->nContigs = vcfd->hdr->n[BCF_DT_CTG];
    ASSERT(vcfd->nContigs > 0);

    check_consistency_args_pars(pars);

    return vcfd;
}

void vcfData_destroy(vcfData* v) {
    bcf_hdr_destroy(v->hdr);
    bcf_destroy(v->rec);


    int BCF_CLOSE = bcf_close(v->in_fp);
    if (0 != BCF_CLOSE) {
        fprintf(stderr, "\n[ERROR]\tbcf_close had non-zero status %d\n", BCF_CLOSE);
        exit(BCF_CLOSE);
    }

    if (v->lngl != NULL) {
        delete v->lngl;
    }

    if (NULL != v->jointGenotypeMatrixGL) {
        v->print_JointGenotypeCountMatrix();
        for (int i = 0;i < v->nIndCmb;++i) {
            FREE(v->jointGenotypeMatrixGL[i]);
        }
        FREE(v->jointGenotypeMatrixGL);
    }

    if (NULL != v->jointGenotypeMatrixGT) {
        v->print_JointGenotypeCountMatrix();
        for (int i = 0;i < v->nIndCmb;++i) {
            FREE(v->jointGenotypeMatrixGT[i]);
        }
        FREE(v->jointGenotypeMatrixGT);
    }

    FREE(v->snSites);

    if (NULL != v->itr) {
        hts_itr_destroy(v->itr);
    }

    if (NULL != v->idx) {
        hts_idx_destroy(v->idx);
    }

    if (NULL != v->tbx) {
        tbx_destroy(v->tbx);
    }

    if (NULL != v->gts) {
        delete(v->gts);
    }

    if (NULL != v->gls) {
        delete(v->gls);
    }

    delete v;
}

void readSites(vcfData* vcfd, paramStruct* pars, blobStruct* blob) {


    int skip_site = 0;

    int contig_i = -1;
    int site_i = 0;

    char prev_contig[1024];

    int nBlocks = 1;
    if (NULL != blob) {
        ASSERT(blob->bootstraps->nBlocks > 1);
        nBlocks = blob->bootstraps->nBlocks;

    } else {
        // if block reading is disabled, read everything into one single block
        nBlocks = 1;
    }

    int nSites_contig_i = 0;

    // read vcf records block by block
    // assume blocks are sorted in the same order as the vcf file
    for (int block_i = 0; block_i < nBlocks; block_i++) {

        int block_start = 0;
        int block_end = 0;
        int nSites_block_i = 0;

        if (nBlocks > 1) {
            block_start = blob->blocks[block_i]->start;
            block_end = blob->blocks[block_i]->end;
            nSites_block_i = 0;
        }

        // read vcf records until the end of the block (or the end of the vcf file)

        // if region is specified, records will not contain any sites outside the region
        // but the number of contigs in the header still includes all contigs in the vcf file

        // TODO
        //  but the block may contain sites outside the region
        while (1 == (vcfd->records_next())) {
            if (NULL != blob) {
                // if site is in block_i
                if (0 == strcmp(blob->blocks[block_i]->chr, bcf_seqname(vcfd->hdr, vcfd->rec)) && block_start <= vcfd->rec->pos && block_end > vcfd->rec->pos) {
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

                    if (0 == strcmp(blob->blocks[block_i]->chr, bcf_seqname(vcfd->hdr, vcfd->rec)) && block_start <= vcfd->rec->pos && block_end > vcfd->rec->pos) {
                        nSites_block_i++;
                    } else {
                        // skip the site
                        continue;
                    }
                }
            }

            if (pars->nSites >= pars->nSites_arrays_size) {
                arrays_nSites_expand(pars, vcfd->lngl, pars->a1a2);
            }

            // if at the very beginning
            if (-1 == contig_i) {
                ++contig_i;
                strncpy(prev_contig, vcfd->get_contig_name(), sizeof(prev_contig));
                if (prev_contig[sizeof(prev_contig) - 1] != '\0') {
                    ERROR("Contig name is too long");
                }
                site_i = 0;
            }

            // ** contig change **
            if (0 != strcmp(vcfd->get_contig_name(), prev_contig)) {
                if (contig_i >= vcfd->nContigs) {
                    ERROR("Number of contigs in the VCF file header is larger than the number of contigs in the VCF file");
                }

                // if we skipped all sites in the previous contig, reuse contig_i index for the new contig
                if (0 == nSites_contig_i) {
                } else {
                    ++contig_i;
                    nSites_contig_i = 0;
                }

                strncpy(prev_contig, vcfd->get_contig_name(), sizeof(prev_contig));
                if (prev_contig[sizeof(prev_contig) - 1] != '\0') {
                    ERROR("Contig name is too long");
                }
                site_i = 0;
                //TODO empty the pars->a1a2
            }

            if (vcfd->lngl != NULL) {
                skip_site = site_read_GL(contig_i, site_i, vcfd, pars);
            } else {
                // TODOdragon
                // skip_site = get_jointGenotypeMatrixGT(contig_i, site_i, vcfd, pars, -1);
                skip_site = site_read_GT(contig_i, site_i, vcfd, pars);
            }

            if (skip_site != 0) {
                IO::vprint(1, "Skipped site at %s:%ld", vcfd->get_contig_name(), vcfd->rec->pos + 1);
                pars->totSites++;

                continue;
            }

            ++site_i;
            pars->nSites++;
            pars->totSites++;
            nSites_contig_i++;
        }
    }

    BEGIN_LOGSECTION;
    LOG("Finished reading sites.");
    LOG("Read %d (out of %d) contigs from the VCF file.", contig_i + 1, vcfd->nContigs);
    LOG("Number of contigs retained: %d", contig_i + 1);
    LOG("Number of contigs skipped: %d", vcfd->nContigs - contig_i - 1);
    LOG("Total number of sites processed: %lu", pars->totSites);
    LOG("Total number of sites skipped for all individual pairs: %lu", pars->totSites - pars->nSites);
    LOG("Total number of sites retained: %lu", pars->nSites);
    END_LOGSECTION;
}

//TODO
//TODO  231231 add condensed form (order G1 G2 vs G2 G1 does not matter; sum these)
void vcfData::print_JointGenotypeCountMatrix() {
    if (outFiles->out_jgcd_fs != NULL) {
        outFiles->out_jgcd_fs->kbuf = kbuf_init();
        kstring_t* kbuf = outFiles->out_jgcd_fs->kbuf;
        if (args->doDist == 1) {
            for (int i = 0; i < nIndCmb; i++) {
                ksprintf(kbuf, "%i,", i);
                for (int j = 0; j < nJointClasses; j++) {
                    ksprintf(kbuf, "%f", jointGenotypeMatrixGL[i][j]);
                    if (j == nJointClasses - 1) {
                        ksprintf(kbuf, "\n");
                    } else {
                        ksprintf(kbuf, ",");
                    }
                }
            }
            // print expectations 
            for (int i = 0; i < nIndCmb; i++) {
                ksprintf(kbuf, "%i,", i);
                for (int j = 0; j < nJointClasses; j++) {
                    double exp = jointGenotypeMatrixGL[i][j] * ((double)snSites[i]);
                    ksprintf(kbuf, "%f", exp);
                    if (j == nJointClasses - 1) {
                        ksprintf(kbuf, "\n");
                    } else {
                        ksprintf(kbuf, ",");
                    }
                }
            }
        } else if (args->doDist == 2) {
            for (int i = 0; i < nIndCmb; i++) {
                ksprintf(kbuf, "%i,", i);
                for (int j = 0; j < nJointClasses; j++) {
                    ksprintf(kbuf, "%i", jointGenotypeMatrixGT[i][j]);
                    if (j == nJointClasses - 1) {
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


void vcfData::_print(FILE* fp) {
    fprintf(stderr, "\nNumber of samples: %i", nInd);
    fprintf(stderr, "\nNumber of contigs: %d", nContigs);
}

void vcfData::_print() {
    _print(stderr);
}

int vcfData::records_next() {
    int ret = -42;

    if (NULL == itr) {
        // no region specified
        ret = bcf_read(in_fp, hdr, rec);
    } else {
        // region reading using iterator
        ret = bcf_itr_next(in_fp, itr, rec);

        if (-1 > ret) {
            ERROR("An error occurred while reading the VCF file.");
        }
    }

    bcf_unpack(rec, DO_BCF_UNPACK);

    ASSERT(0 == rec->errcode);

    if (-1 == ret) {
        // no more data
        return 0;
    }

    return 1;
}

