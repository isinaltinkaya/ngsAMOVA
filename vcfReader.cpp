#include "vcfReader.h"

#include "bootstrap.h"

void bblocks_match_contig(bblocks_t* bblocks, const char* this_contigName, const size_t vcf_contig_i) {
    size_t bb_contig_i;
    if (bblocks->contig_names->find(this_contigName, &bb_contig_i)) {
        if (bb_contig_i != vcf_contig_i) {
            if (PROGRAM_HAS_INPUT_BLOCKS) {
                ERROR("Block contigs do not have the same order as the VCF contigs. Please make sure the contigs in the block file are in the same order as the contigs in the VCF file.");
            } else {
                NEVER; // we generated the blocks; something went wrong
            }
        } // else OK
    } else {
        if (PROGRAM_HAS_INPUT_BLOCKS) {
            ERROR("Contig %s not found in block file", this_contigName);
        } else {
            NEVER; // we generated the blocks; something went wrong
        }
    }
    return;
}

gldata_t* gldata_init(void) {
    gldata_t* gldata = (gldata_t*)malloc(sizeof(gldata_t));
    ASSERT(gldata != NULL);

    gldata->n_gls = 0;
    gldata->n_ind = 0;
    gldata->size = 0;
    gldata->step_size = 0;
    gldata->type = ARG_GLDATA_TYPE_UNMOD;

    gldata->d = NULL;
    gldata->mis = NULL;

    return(gldata);
}

static gldata_t* gldata_alloc(const uint8_t n_gls, const size_t n_ind, const size_t init_n_sites, const uint8_t type) {
    gldata_t* gldata = gldata_init();
    gldata->n_gls = n_gls;
    gldata->n_ind = n_ind;
    gldata->size = init_n_sites;
    gldata->step_size = init_n_sites;
    gldata->type = type;

    if (n_gls != 3) {
        ERROR("Not implemented yet.");
    }

    gldata->d = (double***)malloc(n_ind * sizeof(double**));
    ASSERT(gldata->d != NULL);

    gldata->mis = (bool**)malloc(n_ind * sizeof(bool*));
    ASSERT(gldata->mis != NULL);

    for (size_t i = 0;i < n_ind;++i) {

        gldata->d[i] = (double**)malloc(init_n_sites * sizeof(double*));
        ASSERT(gldata->d[i] != NULL);
        for (size_t s = 0;s < init_n_sites;++s) {
            gldata->d[i][s] = NULL;
            // alloc in-place if !mis[i][s] 
        }

        gldata->mis[i] = NULL;
        gldata->mis[i] = (bool*)malloc(init_n_sites * sizeof(bool));
        ASSERT(gldata->mis[i] != NULL);

        for (size_t j = 0;j < init_n_sites;++j) {
            gldata->mis[i][j] = true; // init to true to exclude the unused sites at the end if array size is larger than the actual number of sites
        }

    }

    return(gldata);
}

static void gldata_realloc(gldata_t* gldata) {

    const size_t new_size = gldata->size + gldata->step_size;

    for (size_t i = 0; i < gldata->n_ind; ++i) {
        gldata->d[i] = (double**)realloc(gldata->d[i], new_size * sizeof(double*));
        ASSERT(gldata->d[i] != NULL);
        for (size_t j = gldata->size; j < new_size; ++j) {
            gldata->d[i][j] = NULL;
        }

        gldata->mis[i] = (bool*)realloc(gldata->mis[i], new_size * sizeof(bool));
        ASSERT(gldata->mis[i] != NULL);
        for (size_t j = gldata->size; j < new_size; ++j) {
            gldata->mis[i][j] = true;
        }
    }

    gldata->size = new_size;
    return;
}


/// @param indends - array of pointers to individual ends (ptr to the last allocated site)
void gldata_destroy(gldata_t* gldata) {
    for (size_t i = 0; i < gldata->n_ind; ++i) {
        for (size_t j = 0; j < gldata->size; ++j) {
            if (NULL != gldata->d[i][j]) {
                // isnull if mis[i][j]
                FREE(gldata->d[i][j]);
            }
        }
        FREE(gldata->d[i]);
        FREE(gldata->mis[i]);
    }
    FREE(gldata->d);
    FREE(gldata->mis);
    FREE(gldata);
    return;
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


// -1   missing
// -2
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

int require_index(paramStruct* pars) {
    if (PROGRAM_HAS_INPUT_VCF) {
        if (NULL != args->in_region) {
            return (IDX_CSI);
        }
    }
    if (PROGRAM_HAS_INPUT_DM) {
        return (IDX_NONE);
    }
    return (IDX_NONE);
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


// return 1: skip site (for all individuals)
// assume ref=allele1 (major/anc) and alt=allele2 (minor/der)
int read_site_with_alleles_ref_alt(vcfData* vcfd, paramStruct* pars, const bool contig_changed) {

    const int64_t this_pos = vcfd->rec->pos;

    const int nAlleles = vcfd->rec->n_allele;

    bool allow_invar_site = false; //TODO determine based on analysistype and args
    if (0 == nAlleles) {
        IO::vprint(2, "Skipping site at %s:%ld. Reason: No data found at site.\n", vcfd->get_contig_name(), this_pos + 1);
        return(1);
    }
    if (!allow_invar_site) {
        if (1 == nAlleles) {
            IO::vprint(2, "Skipping site at %s:%ld. Reason: Site is not variable.\n", vcfd->get_contig_name(), this_pos + 1);
            return(1);
        }
    }


    vcfd->allele_unobserved = -1;
    for (int i = 0; i < nAlleles; ++i) {

        // if allele is single char
        if ('\0' == vcfd->rec->d.allele[i][1]) {
            if ('*' == vcfd->rec->d.allele[i][0]) {
                // '*' symbol = allele missing due to overlapping deletion
                // source: VCF specs v4.3 
                IO::vprint(2, "Skipping site at %s:%ld. Reason: REF allele is set to '*' (allele missing due to overlapping deletion).\n", vcfd->get_contig_name(), this_pos + 1);
                return(1);
            }
            if ('.' == vcfd->rec->d.allele[i][0]) {
                // '.' symbol = MISSING value (no variant)
                // source: VCF specs v4.3 
                IO::vprint(2, "Skipping site at %s:%ld. Reason: REF allele is set to '*' (allele missing due to overlapping deletion).\n", vcfd->get_contig_name(), this_pos + 1);
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
                        ERROR("Found bad allele '%s' at %s:%ld. If you believe this is a valid allele, please contact the developers.\n", vcfd->rec->d.allele[i], vcfd->get_contig_name(), this_pos + 1);
                        break;
                    }
                    j++;
                }
                IO::vprint(2, "Skipping site at %s:%ld. Reason: REF allele is an INDEL.\n", vcfd->get_contig_name(), this_pos + 1);
                return(1);
            }

            ERROR("Found bad allele '%s' at %s:%ld. If you believe this is a valid allele, please contact the developers.\n", vcfd->rec->d.allele[i], vcfd->get_contig_name(), this_pos + 1);
        }

        // -> allele starts with '<'

        if (3 == strlen(vcfd->rec->d.allele[i])) {
            // if unobserved allele <*>
            if ('*' == vcfd->rec->d.allele[i][1] && '>' == vcfd->rec->d.allele[i][2]) {
                vcfd->allele_unobserved = i;
                continue;
            } else {
                ERROR("Found bad allele '%s' at %s:%ld. If you believe this is a valid allele, please contact the developers.\n", vcfd->rec->d.allele[i], vcfd->get_contig_name(), this_pos + 1);
            }


        } else if (9 == strlen(vcfd->rec->d.allele[i])) {
            if ('N' == vcfd->rec->d.allele[i][1] && 'O' == vcfd->rec->d.allele[i][2] && 'N' == vcfd->rec->d.allele[i][3] && '_' == vcfd->rec->d.allele[i][4] && 'R' == vcfd->rec->d.allele[i][5] && 'E' == vcfd->rec->d.allele[i][6] && 'F' == vcfd->rec->d.allele[i][7] && '>' == vcfd->rec->d.allele[i][8]) {
                // if unobserved allele <NON_REF>
                vcfd->allele_unobserved = i;
                continue;
            } else {
                ERROR("Found bad allele '%s' at %s:%ld. If you believe this is a valid allele, please contact the developers.\n", vcfd->rec->d.allele[i], vcfd->get_contig_name(), this_pos + 1);
            }
        }
        NEVER;
    }


    // if REF==ALLELE_UNOBSERVED
    if (0 == vcfd->allele_unobserved) {

        if (PROGRAM_WILL_USE_ALLELES_REF_ALT1) {
            IO::vprint(2, "Skipping site at %s:%ld. Reason: Program will use REF allele as major allele, but REF allele is set to <NON_REF> or <*>.\n", vcfd->get_contig_name(), this_pos + 1);
            return(1);
        }


        if (nAlleles > 1) {
            ERROR("Sites with REF allele set to <NON_REF> or <*> are not supported in the current version of the program.");
            //TODO this is only supported if major minor file or refalt file or an estimation method for these is defined; add proper check and handling
        }

        // check if ALT is set to any allele, maybe site has data but no ref allele specified

        // if 1==nAlleles ALT is not set to any allele, so no data at site
        // this is already skipped above
        if (2 == nAlleles) {
            // ALT is set but REF is not, so probably an invariable site with no REF allele specified
            if (args->rmInvarSites) {
                IO::vprint(2, "Skipping site at %s:%ld. Reason: Site is not variable.\n", vcfd->get_contig_name(), this_pos + 1);
                return(1);
            }
        }

        // b1 = BASE_UNOBSERVED;
        // b2 = acgt_charToInt[(int)vcfd->rec->d.allele[1][0]];

    } else if (1 == vcfd->allele_unobserved) {

        if (PROGRAM_WILL_USE_ALLELES_REF_ALT1) {
            IO::vprint(2, "Skipping site at %s:%ld. Reason: Program will use REF allele as major allele, but ALT allele is set to <NON_REF> or <*>.\n", vcfd->get_contig_name(), this_pos + 1);
            return(1);
        }

        // if ALT==ALLELE_UNOBSERVED

        if (2 == nAlleles) {
            // invariable site, only REF and UNOBSERVED alleles are found
            // if (args->rmInvarSites) {
            IO::vprint(2, "Skipping site at %s:%ld. Reason: Site is not variable.\n", vcfd->get_contig_name(), this_pos + 1);
            return(1);
            // }
        }

        // b1 = acgt_charToInt[(int)vcfd->rec->d.allele[0][0]];
        // b2 = BASE_UNOBSERVED;

    } else if (-1 == vcfd->allele_unobserved) {
        // -> no unobserved alt allele is found

        if (1 == nAlleles) {
            if (args->rmInvarSites) {
                IO::vprint(2, "Skipping site at %s:%ld. Reason: Site is not variable.\n", vcfd->get_contig_name(), this_pos + 1);
                return(1);
            }
        } else {
            // b2 = acgt_charToInt[(int)vcfd->rec->d.allele[1][0]];
            // a2 = 1;
        }

        // b1 = acgt_charToInt[(int)vcfd->rec->d.allele[0][0]];
        // a1 = 0;

    } else {
        // unobserved allele notation is found but is not a1 or a2, so we have at least 3 alleles (0,1, unobserved)
        // b1 = acgt_charToInt[(int)vcfd->rec->d.allele[0][0]];
        // a1 = 0;
        // b2 = acgt_charToInt[(int)vcfd->rec->d.allele[1][0]];
        // a2 = 1;
    }


    if (args->rmInvarSites && !(PROGRAM_WILL_USE_BCF_FMT_GT)) {
        int32_t* ads = NULL;
        int32_t n_ads = 0;
        int32_t n = 0;
        n = bcf_get_format_int32(vcfd->hdr, vcfd->rec, "AD", &ads, &n_ads);
        if (n <= 0) {
            ERROR("Could not read AD tag from the VCF file.");
        }

        // make sure that AD is non-0 for at least 2 alleles 
        int n_non0[nAlleles];
        for (int i = 0; i < nAlleles; ++i) {
            n_non0[i] = 0;
        }
        for (int i = 0; i < pars->names->len; ++i) {
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
            IO::vprint(2, "Skipping site at %s:%ld. Reason: Site is not variable.\n", vcfd->get_contig_name(), this_pos + 1);
            return(1);
        }
    }

    if (PROGRAM_WILL_USE_ALLELES_REF_ALT1) {
        return(0); // pars->a1a2 already initted to {0,1}
    }


    alleles_t* alleles = NULL;

    if (pars->majorminor != NULL) {
        alleles = pars->majorminor;
    } else if (pars->ancder != NULL) {
        alleles = pars->ancder;
    }

    size_t pos_search_start = pars->alleles_posidx + 1; // start from the next position
    if (-1 == pars->alleles_contigidx) {
        // -> first contig, first time we are using alleles
        const char* contigname = bcf_hdr_id2name(vcfd->hdr, vcfd->rec ? vcfd->rec->rid : -1);
        if (NULL == contigname) {
            ERROR("Could not retrieve contig name for record %d.", vcfd->rec->rid);
        }
        size_t cnamesidx;
        if (alleles->cnames->find(contigname, &cnamesidx)) {
            pars->alleles_contigidx = cnamesidx;
            pars->alleles_posidx = alleles->cposidx->d[cnamesidx];
            pos_search_start = pars->alleles_posidx; // start from the first pos of the contig instead
        } else {
            IO::vprint(2, "Skipping site at %s:%ld. Reason: Contig %s not found in alleles file.\n", contigname, this_pos + 1, contigname);
            return(1);
        }
    } else {

        if (contig_changed) {
            const char* contigname = bcf_hdr_id2name(vcfd->hdr, vcfd->rec ? vcfd->rec->rid : -1);
            if (NULL == contigname) {
                ERROR("Could not retrieve contig name for record %d.", vcfd->rec->rid);
            }
            // find the new contig, start with the last contig idx+1
            size_t search_start = pars->alleles_contigidx + 1;
            size_t cnamesidx;
            if (alleles->cnames->find_from(contigname, &cnamesidx, pars->alleles_contigidx)) {
                pars->alleles_contigidx = cnamesidx;
                pars->alleles_posidx = alleles->cposidx->d[cnamesidx];
                pos_search_start = pars->alleles_posidx; // start from the first pos of the contig instead
            } else {
                IO::vprint(2, "Skipping site at %s:%ld. Reason: Contig %s not found in alleles file.\n", contigname, this_pos + 1, contigname);
                return(1);
            }
        }
    }



    // if not the last contig, only search until the start of next contig
    size_t pos_search_end = (pars->alleles_contigidx < alleles->cnames->len - 1) ? alleles->cposidx->d[pars->alleles_contigidx + 1] : alleles->pos->len;

    bool foundSite;

    foundSite = false;
    for (size_t posidx = pos_search_start; posidx < pos_search_end; posidx++) {
        if (alleles->pos->d[posidx] == (size_t)this_pos) {
            pars->alleles_posidx = posidx;
            foundSite = true;
            break;
        }
    }
    if (!foundSite) {
        IO::vprint(2, "Skipping site at %s:%ld. Reason: Site not found in alleles file.\n", vcfd->get_contig_name(), this_pos + 1);
        return(1);
    }

    char allele1[2];
    char allele2[2];

    // the 64bit block containing the alleles for the position
    uint64_t tmp = alleles->d[(pars->alleles_posidx / 16)];
    // the index of the position in the 64bit block
    uint64_t j = (pars->alleles_posidx % 16) * 4;
    int baseint1 = (tmp >> j) & 3;
    int baseint2 = (tmp >> (j + 2)) & 3;
    *allele1 = "ACGT"[baseint1];
    *allele2 = "ACGT"[baseint2];


    // -> reset
    pars->a1a2[0] = -1;
    pars->a1a2[1] = -1;

    // index of vcf allele in vcf->rec->d.allele that is ref according to alleles_t
    for (size_t vcfal = 0; vcfal < vcfd->rec->n_allele; ++vcfal) {
        // no need to check if allele is not single char
        if ('\0' == vcfd->rec->d.allele[vcfal][1]) {
            if (vcfd->rec->d.allele[vcfal][0] == *allele1) {
                pars->a1a2[0] = vcfal;
                // DEVPRINT("The allele1 (%c)'s position in vcf rec->d.allele (REF:%s ALT:%s) is %d", *allele1, vcfd->rec->d.allele[0], vcfd->rec->d.allele[1], pars->a1a2[0]);
            }
            if (vcfd->rec->d.allele[vcfal][0] == *allele2) {
                pars->a1a2[1] = vcfal;
                // DEVPRINT("The allele2 (%c)'s position in vcf rec->d.allele (REF:%s ALT:%s) is %d", *allele2, vcfd->rec->d.allele[0], vcfd->rec->d.allele[1], pars->a1a2[1]);
            }
        }
    }

    // -> make sure allele1 allele2 exists in vcfd->rec->d.allele
    if (pars->a1a2[0] == -1) {
        IO::vprint(2, "Skipping site at %s:%ld. Reason: Allele1 (%c) not found in vcf rec->d.allele (REF:%s ALT:%s).\n", vcfd->get_contig_name(), this_pos + 1, *allele1, vcfd->rec->d.allele[0], vcfd->rec->d.allele[1]);
        return(1);
    }
    if (pars->a1a2[1] == -1) {
        IO::vprint(2, "Skipping site at %s:%ld. Reason: Allele2 (%c) not found in vcf rec->d.allele (REF:%s ALT:%s).\n", vcfd->get_contig_name(), this_pos + 1, *allele2, vcfd->rec->d.allele[0], vcfd->rec->d.allele[1]);
        return(1);
    }

    ASSERT(pars->a1a2[0] != pars->a1a2[1]);

    return(0);
}

// return 1: skip site for all individuals
static int site_read_3GL(vcfData* vcfd, paramStruct* pars, const size_t block_i, const bool contig_changed) {




    int ret;
    if (0 != (ret = read_site_with_alleles_ref_alt(vcfd, pars, contig_changed))) {
        return(ret);
    }

    const size_t site_idx = pars->nSites;
    const size_t gldata_n_gls = vcfd->gldata->n_gls;


    const int ngt_to_use = vcfd->nGT;
    // DEVPRINT("n_gls: %d nGT:%d", lgls.n_gls, vcfd->nGT);
    // ngls is observed number of genotypes as they appear in gl data
    // ngt_to_use is what we need to construct the pairwise genotype counts matrix (i.e. 3 or 10)

    float* lgls = NULL;
    int size_e = 0;
    int n_gls = 0;
    int n_missing_ind = 0;
    int n_values = bcf_get_format_float(vcfd->hdr, vcfd->rec, "GL", &lgls, &size_e);
    if (n_values <= 0) {
        ERROR("Could not read GL tag from the VCF file.");
    }
    n_gls = n_values / pars->names->len;


    if (n_gls < ngt_to_use) {
        IO::vprint(2, "Skipping site at %s:%ld. Reason: Number of GLs is less than the number of GLs needed to perform the specified analysis (%d).\n", vcfd->get_contig_name(), vcfd->rec->pos + 1, ngt_to_use);
        return(1);
    }

    float* sample_lgls = NULL;


    int skipSite = 0;

    bool** isMissing = vcfd->gldata->mis;

    uint64_t ordered_gt_fetch[3];
    ordered_gt_fetch[0] = bcf_alleles2gt(pars->a1a2[0], pars->a1a2[0]);
    ordered_gt_fetch[1] = bcf_alleles2gt(pars->a1a2[0], pars->a1a2[1]);
    ordered_gt_fetch[2] = bcf_alleles2gt(pars->a1a2[1], pars->a1a2[1]);

    double* iptr;

    double max;
    const size_t nInd = pars->names->len;

    for (size_t i = 0; i < nInd; i++) {

        isMissing[i][site_idx] = false;

        sample_lgls = lgls + i * n_gls;

        // ----------------------------------------------------------------- //
        // -> check for missing data
        if (bcf_float_is_missing(sample_lgls[0])) {
            // missing data check 1
            isMissing[i][site_idx] = true;

        } else {
            // missing data check 2
            int z = 1;
            for (int j = 1; j < n_gls; ++j) {
                if (sample_lgls[0] == sample_lgls[j]) {
                    ++z;
                }
            }
            if (z == n_gls) {
                // missing data (all gl values are the same)
                isMissing[i][site_idx] = true;
            }
        }

        if (isMissing[i][site_idx]) {
            // if minInd 0 (the program will only use sites shared across all ividuals) 
            // skip site when you encounter any missing data
            if (0 == args->minInd) {
                skipSite = 1;
                break;
            }

            // skip site if there are 2 inds and at least one of them is missing
            if (2 == nInd) {
                skipSite = 1;
                break;
            }

            n_missing_ind++;

            if (i == nInd - 1 && nInd == n_missing_ind) {
                // check if all ividuals are missing
                skipSite = 1;
                break;
            }

            // skip site if minInd is defined and #non-missing inds=<nInd
            if (args->minInd != 2) {
                if ((nInd - n_missing_ind) < args->minInd) {
                    IO::vprint(2, "Skipping site at %s:%ld. Reason: Number of non-missing ividuals is less than the minimum number of ividuals -minInd.\n", vcfd->get_contig_name(), vcfd->rec->pos + 1);
                    skipSite = 1;
                    break;
                }
            }
            continue;
        }


        // -> not missing
        DEVASSERT(!isMissing[i][site_idx]);

        sample_lgls = lgls + (i * n_gls);


        vcfd->gldata->d[i][site_idx] = (double*)malloc(n_gls * sizeof(double));
        ASSERT(vcfd->gldata->d[i][site_idx] != NULL);

        iptr = vcfd->gldata->d[i][site_idx];

        if (vcfd->gldata->type == ARG_GLDATA_TYPE_UNMOD) {
            iptr[0] = (double)sample_lgls[ordered_gt_fetch[0]];
            iptr[1] = (double)sample_lgls[ordered_gt_fetch[1]];
            iptr[2] = (double)sample_lgls[ordered_gt_fetch[2]];
        } else {

            if (vcfd->gldata->type == ARG_GLDATA_TYPE_LOG10GL) {
                iptr[0] = (double)sample_lgls[ordered_gt_fetch[0]];
                iptr[1] = (double)sample_lgls[ordered_gt_fetch[1]];
                iptr[2] = (double)sample_lgls[ordered_gt_fetch[2]];

            } else if (vcfd->gldata->type == ARG_GLDATA_TYPE_LNGL) {
                iptr[0] = (double)LOG2LN(sample_lgls[ordered_gt_fetch[0]]);
                iptr[1] = (double)LOG2LN(sample_lgls[ordered_gt_fetch[1]]);
                iptr[2] = (double)LOG2LN(sample_lgls[ordered_gt_fetch[2]]);

            } else {
                NEVER;
            }

            max = iptr[0];
            for (size_t i = 1; i < 3; i++) {
                max = (max < iptr[i]) ? iptr[i] : max;
            }
            for (size_t i = 0; i < 3; i++) {
                iptr[i] -= max;
            }
        }


    }


    FREE(lgls);
    return(skipSite);
}


// return 1: skip site for all individuals
int site_read_GT(jgtmat_t* jgtm, vcfData* vcfd, paramStruct* pars, const bool contig_changed) {

    int ret;
    if (0 != (ret = read_site_with_alleles_ref_alt(vcfd, pars, contig_changed))) {
        return(ret);
    }

    const int nInd = pars->names->len;
    int n_missing_ind = 0;
    int32_t* gts = NULL;
    int size_e = 0;

    int n_values = bcf_get_genotypes(vcfd->hdr, vcfd->rec, &gts, &size_e);
    if (n_values <= 0) {
        ERROR("Could not read GT tag from the VCF file.");
    }
    if ((n_values / pars->names->len) != PROGRAM_PLOIDY) {
        ERROR("Ploidy %d not supported.", n_values / pars->names->len);
    }

    int32_t* sample_gts = NULL;
    bool hasData[nInd];

    int skipSite = 0;


    int32_t* i1_gts = NULL;
    int32_t* i2_gts = NULL;
    int8_t i1_nder = -1;
    int8_t i2_nder = -1;


    for (int indi = 0; indi < nInd; ++indi) {

        hasData[indi] = false;

        sample_gts = gts + indi * PROGRAM_PLOIDY;

        if (1 == bcf_gt_is_missing(sample_gts[0]) || 1 == bcf_gt_is_missing(sample_gts[1])) {
            // minInd 0
            // only use sites shared across all individuals
            // so skip site when you first encounter nan
            if (0 == args->minInd) {
                skipSite = 1;
                break;
            }

            // skip site if there are 2 inds and at least one of them is missing
            if (2 == nInd) {
                skipSite = 1;
                break;
            }

            n_missing_ind++;

            // all individuals are missing
            if (nInd == n_missing_ind) {
                skipSite = 1;
                break;
            }

            // skip site if minInd is defined and #non-missing inds=<nInd
            if (args->minInd != 2) {
                if ((nInd - n_missing_ind) < args->minInd) {
                    IO::vprint(2, "Skipping site at %s:%ld. Reason: Number of non-missing individuals is less than the minimum number of individuals -minInd.\n", vcfd->get_contig_name(), vcfd->rec->pos + 1);
                    skipSite = 1;
                    break;
                }
            }
            continue;
        } else {
            // -> not missing
            hasData[indi] = true;
        }
    }

    bool skipInd;

    if (!skipSite && args->rmInvarSites != 0) {
        // -> rm invar sites based on GT
        // @dragon
        int tot_nder = 0;
        for (int i = 0; i < nInd; ++i) {
            if (!hasData[i]) {
                continue;
            }
            skipInd = false;
            i1_gts = gts + i * PROGRAM_PLOIDY;
            if (bcf_gt_allele(i1_gts[0]) == pars->a1a2[1]) {
                tot_nder++;
            } else if (bcf_gt_allele(i1_gts[0]) == pars->a1a2[0]) {
                //
            } else {
                skipInd = true;
                continue;
            }
            if (bcf_gt_allele(i1_gts[1]) == pars->a1a2[1]) {
                tot_nder++;
            } else if (bcf_gt_allele(i1_gts[1]) == pars->a1a2[0]) {
                //
            } else {
                skipInd = true;
                continue;
            }
        }
        if (tot_nder == 0) {
            IO::vprint(2, "Skipping site at %s:%ld. Reason: Site is not variable.\n", vcfd->get_contig_name(), vcfd->rec->pos + 1);
            skipSite = 1;
        } else if (tot_nder == nInd * PROGRAM_PLOIDY) {
            IO::vprint(2, "Skipping site at %s:%ld. Reason: Site is not variable.\n", vcfd->get_contig_name(), vcfd->rec->pos + 1);
            skipSite = 1;
        }
    }


    do {
        if (skipSite) {
            break;
        }

        int gt1 = -1;
        int gt2 = -1;
        int pidx = 0;
        for (int i1 = 1; i1 < nInd; ++i1) {

            if (!hasData[i1]) {
                ++pidx;
                continue;
            }

            skipInd = false;
            i1_nder = 0;
            i1_gts = gts + i1 * PROGRAM_PLOIDY;
            gt1 = bcf_gt_allele(i1_gts[0]);
            gt2 = bcf_gt_allele(i1_gts[1]);


            if (gt1 == pars->a1a2[0]) {
                // first base in genotype is ref allele
            } else if (gt1 == pars->a1a2[1]) {
                i1_nder++;
            } else {
                skipInd = true;
            }

            if (gt2 == pars->a1a2[0]) {
                // second base in genotype is ref allele
            } else if (gt2 == pars->a1a2[1]) {
                i1_nder++;
            } else {
                skipInd = true;
            }

            if (skipInd) {
                IO::vprint(2, "Skipping individual %d (%s) at site %s:%ld. Reason: Allele in genotype (%c) is not one of alleles to use (%c %c).", i1, vcfd->hdr->samples[i1], vcfd->get_contig_name(), vcfd->rec->pos + 1, vcfd->rec->d.allele[gt1][0], vcfd->rec->d.allele[0][0], vcfd->rec->d.allele[1][0]);
                hasData[i1] = false;
                ++pidx;
                continue;
            }

            for (int i2 = 0;i2 < i1; ++i2) {

                if (!hasData[i2]) {
                    ++pidx;
                    continue;
                }

                skipInd = false;
                i2_nder = 0;
                i2_gts = gts + i2 * PROGRAM_PLOIDY;
                gt1 = bcf_gt_allele(i2_gts[0]);
                gt2 = bcf_gt_allele(i2_gts[1]);

                if (gt1 == pars->a1a2[0]) {
                    // first base in genotype is ref allele
                } else if (gt1 == pars->a1a2[1]) {
                    i2_nder++;
                } else {
                    skipInd = true;
                }

                if (gt2 == pars->a1a2[0]) {
                    // second base in genotype is ref allele
                } else if (gt2 == pars->a1a2[1]) {
                    i2_nder++;
                } else {
                    skipInd = true;
                }

                if (skipInd) {
                    IO::vprint(2, "Skipping individual %d (%s) at site %s:%ld. Reason: Allele in genotype (%c) is not one of alleles to use (%c %c).", i2, vcfd->hdr->samples[i2], vcfd->get_contig_name(), vcfd->rec->pos + 1, vcfd->rec->d.allele[gt1][0], vcfd->rec->d.allele[0][0], vcfd->rec->d.allele[1][0]);
                    hasData[i2] = false;
                    ++pidx;
                    continue;
                }

                jgtm->m[pidx][(i1_nder * 3) + i2_nder]++;
                ++pidx;
            }

        }

    } while (0);

    FREE(gts);
    return(skipSite);
}

vcfData* vcfData_init(paramStruct* pars, metadataStruct* metadata) {

    BEGIN_LOGSECTION;

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


    const size_t in_vcf_nind = bcf_hdr_nsamples(vcfd->hdr);

    if (metadata == NULL) {
        pars->names = strArray_alloc(in_vcf_nind);
        for (size_t i = 0;i < in_vcf_nind;++i) {
            pars->names->add(vcfd->hdr->samples[i]);
        }
        DEVASSERT(pars->names->len == in_vcf_nind);
        LOGADD("Found %d individuals in the VCF file", pars->names->len);
    } else {
        // -> set pars with vcfdata
        const size_t mtd_nind = pars->metadata->indNames->len;

        size_t midx;

        bool badorder = false;
        uint64_t samples_metad2vcf[mtd_nind];

        size_t newvidx = 0;
        for (size_t vidx = 0;vidx < in_vcf_nind;++vidx) {
            if (pars->metadata->indNames->find(vcfd->hdr->samples[vidx], &midx)) {
                samples_metad2vcf[midx] = vidx;
                if (newvidx != midx) {
                    badorder = true;
                    break;
                }
                ++newvidx;
            } else {
                LOGADD("Skipping individual %s in the VCF file. Reason: Individual not present in metadata file.", vcfd->hdr->samples[vidx]);
            }
        }
        if (newvidx != mtd_nind) {
            ERROR("All individuals listed in metadata file must be present in the VCF file. Program could only find %ld out of %ld individuals.", newvidx, mtd_nind);
        }

        if (badorder) {
            ERROR("The order of individuals in the VCF file does not match the order of individuals in the metadata file. Please sort the individuals in the metadata file to match the order of individuals in the VCF file.");
        }

        pars->names = pars->metadata->indNames;

        kstring_t tmp = KS_INITIALIZE;
        for (size_t i = 0;i < mtd_nind;++i) {
            ksprintf(&tmp, "%s", vcfd->hdr->samples[samples_metad2vcf[i]]);
            if (i != mtd_nind - 1) {
                kputc(',', &tmp);
            }
        }
        if (bcf_hdr_set_samples(vcfd->hdr, tmp.s, 0)) {
            ERROR("Could not set samples in bcf header.");
        }
        FREE(tmp.s);
        int nsamples = bcf_hdr_nsamples(vcfd->hdr); // sanity check
        ASSERT(nsamples == pars->metadata->indNames->len);

        if (nsamples != in_vcf_nind) {
            LOGADD("Will use %d individuals (out of %d) found in the VCF file.", nsamples, in_vcf_nind);
        } else {
            LOGADD("Found %d individuals in the VCF file", nsamples);
        }

    }





    vcfd->nContigs = vcfd->hdr->n[BCF_DT_CTG];
    ASSERT(vcfd->nContigs > 0);

    if (args->minInd == pars->names->len) {
        DEVPRINT("-minInd %d is equal to the number of individuals found in file: %d. Setting -minInd to 0 (all).\n", args->minInd, pars->names->len);
        args->minInd = 0;
    }

    if (pars->names->len == 1) {
        ERROR("Only one individual found in the VCF file. At least two individuals are required to perform the analysis.");
    }

    if (pars->names->len < args->minInd) {
        ERROR("Number of individuals in the VCF file is less than the minimum number of individuals provided by -minInd (%d).", args->minInd);
    }


    if (PROGRAM_WILL_USE_BCF_FMT_GL) {
        int tag_id = bcf_hdr_id2int(vcfd->hdr, BCF_DT_ID, "GL");
        if (!bcf_hdr_idinfo_exists(vcfd->hdr, BCF_HL_FMT, tag_id)) {
            ERROR("GL tag not found in the VCF file.");
        }
    }

    if (PROGRAM_WILL_USE_BCF_FMT_GT) {
        int tag_id = bcf_hdr_id2int(vcfd->hdr, BCF_DT_ID, "GT");
        if (!bcf_hdr_idinfo_exists(vcfd->hdr, BCF_HL_FMT, tag_id)) {
            ERROR("GT tag not found in the VCF file.");
        }
    }
    if (args->rmInvarSites && !(PROGRAM_WILL_USE_BCF_FMT_GT)) {
        int tag_id = bcf_hdr_id2int(vcfd->hdr, BCF_DT_ID, "AD");
        if (!bcf_hdr_idinfo_exists(vcfd->hdr, BCF_HL_FMT, tag_id)) {
            ERROR("AD tag not found in the VCF file.");
        }
    }

    END_LOGSECTION;

    return(vcfd);
}

void vcfData_destroy(vcfData* v) {
    bcf_hdr_destroy(v->hdr);
    bcf_destroy(v->rec);


    int BCF_CLOSE = bcf_close(v->in_fp);
    if (0 != BCF_CLOSE) {
        fprintf(stderr, "\n[ERROR]\tbcf_close had non-zero status %d\n", BCF_CLOSE);
        exit(BCF_CLOSE);
    }

    if (v->gldata != NULL) {
        gldata_destroy(v->gldata);
    }

    if (NULL != v->itr) {
        hts_itr_destroy(v->itr);
    }

    if (NULL != v->idx) {
        hts_idx_destroy(v->idx);
    }

    if (NULL != v->tbx) {
        tbx_destroy(v->tbx);
    }

    delete v;
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




void readSites(jgtmat_t* jgtm, bblocks_t* bblocks, vcfData* vcfd, paramStruct* pars) {

    BEGIN_LOGSECTION;

    int skip_site = 0;


    char prev_contig[1024];

    const size_t nInd = pars->names->len;



    // while we read the data into memory, all we do for bootstrap is to keep track of the number of sites in each block


    const bool doBlocks = (bblocks != NULL);
    size_t nBlocks = (doBlocks) ? bblocks->n_blocks : 1;
    uint64_t* nsites_per_block = NULL;
    if (doBlocks) {
        nsites_per_block = bblocks->nsites_per_block;
    }

    if (args->bcfSrc & ARG_INTPLUS_BCFSRC_FMT_GL) {
        vcfd->gldata = gldata_alloc(vcfd->nGT, pars->names->len, 4096, ARG_GLDATA_TYPE_LNGL);
    }



    const char* this_contigName = NULL;

    int contig_i = -1;
    size_t block_i = 0;

    bool isFirstRecord = true;
    bool newBlock;

    bool contig_changed;

    size_t block_start_pos;
    size_t block_end_pos;
    size_t this_pos;



    while (1 == (vcfd->records_next())) {

        if (vcfd->rec->pos >= 0) {
            this_pos = vcfd->rec->pos;
        } else {
            NEVER;
        }

        if (vcfd->gldata != NULL && pars->nSites == vcfd->gldata->size) {
            gldata_realloc(vcfd->gldata);
        } else if (vcfd->gldata != NULL && pars->nSites > vcfd->gldata->size) {
            NEVER;
        }


        newBlock = false;
        contig_changed = false;

        ASSERT(NULL != (this_contigName = bcf_hdr_id2name(vcfd->hdr, vcfd->rec ? vcfd->rec->rid : -1)));

        if (!isFirstRecord) {
            if (0 != strcmp(prev_contig, this_contigName)) {
                // ** contig change **
                contig_changed = true;
                ++contig_i;
                strncpy(prev_contig, this_contigName, sizeof(prev_contig));
                if (prev_contig[sizeof(prev_contig) - 1] != '\0') {
                    ERROR("Contig name is too long");
                }

                if (doBlocks) {
                    bblocks_match_contig(bblocks, this_contigName, contig_i);

                    newBlock = true;
                }
            }

        } else {
            newBlock = true;
            // if at the very beginning; first ever pos
            contig_changed = false;
            contig_i = 0;
            strncpy(prev_contig, this_contigName, sizeof(prev_contig));
            ASSERT(prev_contig[sizeof(prev_contig) - 1] == '\0');
            isFirstRecord = false;


            if (doBlocks) {
                bblocks_match_contig(bblocks, this_contigName, contig_i);
            }

        }

        // ------------------------------------------------- //
        // -> read data into memory


        // -> first, if doBlocks, find which block we are in
        if (doBlocks) {




            while (1) {
                block_start_pos = bblocks->block_start_pos[block_i];
                block_end_pos = bblocks->block_end_pos[block_i];

                if (block_i < nBlocks) {
                    if (this_pos > block_start_pos) {
                        if (this_pos >= block_end_pos) {
                            // inc '=' since end is exclusive
                            ++block_i;
                            newBlock = true;
                            continue;
                        } else if (this_pos < block_end_pos) {
                            break;
                        }

                    } else if (this_pos == block_start_pos) {
                        break;
                    } else {
                        NEVER;
                    }
                } else if (block_i == nBlocks) {
                    // out of blocks but still not found
                    ERROR("Could not find the block for the current position %ld in contig %s", this_pos + 1, this_contigName);
                } else {
                    NEVER;
                }

            }

            ASSERT(bblocks->block_contig[block_i] == contig_i);

            // DEVPRINT("Position %ld is located at block %ld (%s:%ld-%ld)", this_pos + 1, block_i, this_contigName, block_start_pos, block_end_pos);

            if (newBlock) {
                bblocks->block_start_siteidx[block_i] = pars->nSites;
                bblocks->block_contig[block_i] = contig_i;
            }

        }
        // <- block bootstrapping prep


        if (doBlocks && !newBlock && 0 == nsites_per_block[block_i]) {
            //TODO add testcase forthis
            // a previous site that belongs to this block was found before but was skipped for all indvs
            // so treat it as a new block
            newBlock = true;
        }

        // -> then, read the needed bcf format/data fields 
        if (PROGRAM_WILL_USE_BCF_FMT_GL) {
            skip_site = site_read_3GL(vcfd, pars, block_i, contig_changed);
        } else if (PROGRAM_WILL_USE_BCF_FMT_GT) {
            if (args->doJGTM) {
                skip_site = site_read_GT(jgtm, vcfd, pars, contig_changed);
            } else {
                NEVER;
            }
        } else {
            NEVER;
        }

        // -> skip site if needed
        if (skip_site != 0) {
            IO::vprint(1, "Skipped site at %s:%ld", vcfd->get_contig_name(), this_pos + 1);
            pars->totSites++;
            continue;
        }

        // -> if site is NOT skipped
        if (doBlocks) {
            nsites_per_block[block_i]++;
        }
        pars->nSites++;
        pars->totSites++;

    } // end sites loop




    if (doBlocks) {
        for (size_t i = 0;i < nBlocks;++i) {
            if (nsites_per_block[i] == 0) {
                ERROR("Block %ld has no sites", i);
            }
            if (PROGRAM_VERBOSITY_LEVEL >= 1) {
                LOG("Block %ld has %ld sites", i, nsites_per_block[i]);

            }
        }
    }


    if (args->bcfSrc & ARG_INTPLUS_BCFSRC_FMT_GT) {
        // jgtm->size
        size_t pair_shared_nSites;
        for (size_t pairidx = 0;pairidx < jgtm->n;++pairidx) {
            pair_shared_nSites = 0;
            for (size_t i = 0;i < jgtm->size;++i) {
                pair_shared_nSites += jgtm->m[pairidx][i];
            }
            if (0 == pair_shared_nSites) {
                ERROR("Pair %ld has no shared sites", pairidx);
            }
        }
    }

    LOG("Finished reading sites.");
    LOG("Read %d (out of %d) contigs from the VCF file.", contig_i + 1, vcfd->nContigs);
    LOG("Number of contigs retained: %d", contig_i + 1);
    LOG("Number of contigs skipped: %d", vcfd->nContigs - contig_i - 1);
    LOG("Total number of sites processed: %lu", pars->totSites);
    LOG("Total number of sites skipped for all individual pairs: %lu", pars->totSites - pars->nSites);
    LOG("Total number of sites retained: %lu", pars->nSites);
    END_LOGSECTION;
}
