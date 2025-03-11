#include "vcfReader.h"
#include <math.h>
#include "jgtmat.h"
#include "ibd.h"

#include "bootstrap.h"
#include "metadata.h"


static void launch_percontig_jobs(ibds_t* ibds, paramStruct* pars) {
    if (args->doIbd) {

        // -> analyse ibd results for the data we have in total for the last contig before the current pos (pos where contig change has occured)
        if (args->ibd_alpha < 0.0 && args->ibd_beta < 0.0) {
            identify_ibd_segments_ibdseq(ibds, pars);
        } else {
            identify_ibd_segments(ibds, pars);
        }

        // -> reset and prepare for the new contig
        ibds_reset_percontig(ibds);
    }
    return;
}

typedef enum {
    SKIP_SITE_REASON_MININD = 0,

    SKIP_SITE_REASON_MIN_A2C,

    SKIP_SITE_REASON_MIN_A2F,

    SKIP_SITE_REASON_INVAR,

    SKIP_SITE_REASON_OVERLAP_DEL,

    SKIP_SITE_REASON_INDEL,

    SKIP_SITE_REASON_DOMAJORMINOR_REF,
    // program will use REF as major but REF is set to <NON_REF> or <*>

    SKIP_SITE_REASON_DOMAJORMINOR_ALT,
    // program will use ALT as major but ALT is set to <NON_REF> or <*>

    SKIP_SITE_REASON_ALLELES_CONTIG_NOMATCH,
    // contig name in vcf is not present in alleles file

    SKIP_SITE_REASON_ALLELES_SITE_NOMATCH,
    // site in vcf is not present in alleles file

    SKIP_SITE_REASON_ALLELES_REF_BASE_NOMATCH,
    // ref base specified for the site in alleles file is not one of the alleles observed in the vcf site

    SKIP_SITE_REASON_ALLELES_ALT_BASE_NOMATCH,
    // alt base specified for the site in alleles file is not one of the alleles observed in the vcf site

    SKIP_SITE_REASON_N_GLS,

    SKIP_SITE_REASON_MULTIALLELIC,

    SKIP_SITE_REASON_NOSKIP
    // signals that there is no reason to skip the site, so the site will NOT be skipped!

} skip_site_reason;

const char* skip_site_reason_strs[SKIP_SITE_REASON_NOSKIP] = {
    // SKIP_SITE_REASON_MININD
    "(--min-ind) Minimum number of individuals",

    // SKIP_SITE_REASON_MIN_A2C
    "(--min-a2c) Minimum total a2 allele count for a site to be included in analyses",

    // SKIP_SITE_REASON_MIN_A2F
    "(--min-a2f) Minimum a2 allele frequency for a site to be included in analyses",

    // SKIP_SITE_REASON_INVAR
    "(--rm-invar-sites) Invariable site",

    // SKIP_SITE_REASON_OVERLAP_DEL
    "Overlapping deletion",

    // SKIP_SITE_REASON_INDEL
    "Indel",

    // SKIP_SITE_REASON_DOMAJORMINOR_REF
    "Ref allele set to <NON_REF> or <*>",

    // SKIP_SITE_REASON_DOMAJORMINOR_ALT
    "Alt allele set to <NON_REF> or <*>",

    // SKIP_SITE_REASON_ALLELES_CONTIG_NOMATCH
    "Contig name not found in alleles file",

    // SKIP_SITE_REASON_ALLELES_SITE_NOMATCH
    "Site not found in alleles file",

    // SKIP_SITE_REASON_ALLELES_REF_BASE_NOMATCH
    "Ref base not found in alleles file",

    // SKIP_SITE_REASON_ALLELES_ALT_BASE_NOMATCH
    "Alt base not found in alleles file",

    // SKIP_SITE_REASON_MULTIALLELIC
    "(--rm-multiallelic-sites) Multiallelic site"

    // SKIP_SITE_REASON_N_GLS
    "Number of GLs",
};

static const char* get_skip_site_reason_str(skip_site_reason reason) {
    if (reason >= 0 && reason < SKIP_SITE_REASON_NOSKIP) {
        return (skip_site_reason_strs[reason]);
    }
    ERROR("Unknown skip site reason (%d).", reason);
}

static void bblocks_match_contig(bblocks_t* bblocks, const char* this_contigName, const size_t vcf_contig_i) {
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

    val = BCF_UN_ALL;
    return(val);
}

int require_index(void) {
    if (PROGRAM_HAS_INPUT_VCF) {
        if (args->doBlockBootstrap) {
            return (IDX_CSI);
        }
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
// TODO!! : check evverything about me
static skip_site_reason read_site_a1a2(vcfData* vcfd, paramStruct* pars, const bool contig_changed) {
    //TODO N.B. inform user that POS are 0-based like bed files
    // e.g. bcftools query -f "%CHROM\t%POS0\t%REF\t%ALT"

    static const size_t nInd = pars->nInd;

    const int64_t this_pos = vcfd->rec->pos;
    const int nAlleles = vcfd->rec->n_allele;

    if (nAlleles > 5) {
        ERROR("Found %d alleles at %s:%ld. Program supports only up to 5 alleles.", nAlleles, vcfd->get_contig_name(), this_pos + 1);
    }

    // ------------------------------------------------
    // FILTER: --rm-multiallelic-sites 
    if (args->rmMultiallelicSites) {
        if (nAlleles > 2) {
            return(SKIP_SITE_REASON_MULTIALLELIC);
        }
    }

    // ------------------------------------------------
    // FILTER: --rm-invar-sites (FIRST PASS; SECOND AFTER DETERMINING MAJOR/MINOR)
    if (args->rmInvarSites) {
        if (nAlleles < 2) {
            return(SKIP_SITE_REASON_INVAR);
        }
    }

    // static const bool allow_invar_site = false; //TODO determine based on analysistype and args
    // if (0 == nAlleles) {
    //     NEVER;
    //     // IO::vprint(2, "Skipping site at %s:%ld. Reason: No data found at site.\n", vcfd->get_contig_name(), this_pos + 1);
    // }
    // if (!allow_invar_site) {
    //     if (1 == nAlleles) {
    //         return(SKIP_SITE_REASON_INVAR);
    //     }
    // }



    vcfd->allele_unobserved = -1;
    for (int i = 0; i < nAlleles; ++i) {

        // if allele is single char
        if ('\0' == vcfd->rec->d.allele[i][1]) {
            if ('*' == vcfd->rec->d.allele[i][0]) {
                // '*' symbol = allele missing due to overlapping deletion
                // source: VCF specs v4.3 
                return(SKIP_SITE_REASON_OVERLAP_DEL);
            }
            if ('.' == vcfd->rec->d.allele[i][0]) {
                // '.' symbol = MISSING value (no variant)
                // source: VCF specs v4.3 
                if (args->rmInvarSites) {
                    return(SKIP_SITE_REASON_INVAR);
                }
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
                    }
                    j++;
                }
                return(SKIP_SITE_REASON_INDEL);
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
            return(SKIP_SITE_REASON_DOMAJORMINOR_REF);
        }

        if (nAlleles > 1) {
            ERROR("Sites with REF allele set to <NON_REF> or <*> are not supported.");
            //TODO this is only supported if major minor file or refalt file or an estimation method for these is defined; add proper check and handling
        }

        // -> check if ALT is set to any allele, maybe site has data but no ref allele specified
        // if 1==nAlleles ALT is not set to any allele, so no data at site
        // this is already skipped above
        if (2 == nAlleles) {
            // ALT is set but REF is not, so probably an invariable site with no REF allele specified
            if (args->rmInvarSites) {
                return(SKIP_SITE_REASON_INVAR);
            }
        }

        // b1 = BASE_UNOBSERVED;
        // b2 = acgt_charToInt[(int)vcfd->rec->d.allele[1][0]];

    } else if (1 == vcfd->allele_unobserved) {
        // if ALT[0] is ALLELE_UNOBSERVED

        if (PROGRAM_WILL_USE_ALLELES_REF_ALT1) {
            if (PROGRAM_WILL_USE_ALLELES_REF_ALT1) {
                return(SKIP_SITE_REASON_DOMAJORMINOR_ALT);
            }
        }

        if (2 == nAlleles) {
            if (args->rmInvarSites) {
                return(SKIP_SITE_REASON_INVAR);
            }
        }


    } else if (-1 == vcfd->allele_unobserved) {
        // -> no unobserved alt allele is found

        if (1 == nAlleles) {
            if (args->rmInvarSites) {
                return(SKIP_SITE_REASON_INVAR);
            }
            //} else {
                // b2 = acgt_charToInt[(int)vcfd->rec->d.allele[1][0]];
                // a2 = 1;
        }

        // b1 = acgt_charToInt[(int)vcfd->rec->d.allele[0][0]];
        // a1 = 0;

    //} else {
        // unobserved allele notation is found but is not a1 or a2, so we have at least 3 alleles (0,1, unobserved)
        // b1 = acgt_charToInt[(int)vcfd->rec->d.allele[0][0]];
        // a1 = 0;
        // b2 = acgt_charToInt[(int)vcfd->rec->d.allele[1][0]];
        // a2 = 1;
    }

    int32_t* ac_arr = NULL, ac_n = 0, ac_ndst = 0;
    int32_t* an_arr = NULL, an_ndst = 0, an_n = 0, an_val = 0;


    if (args->doMajorMinor == ARG_DOMAJORMINOR_BCF_REFALT1) {

        DEVASSERT(pars->a1a2[0] == 0);
        DEVASSERT(pars->a1a2[1] == 1);


    } else if (args->doMajorMinor == ARG_DOMAJORMINOR_INFILE) {

        static alleles_t* alleles = NULL;

        if (pars->majorminor != NULL) {
            alleles = pars->majorminor;
        } else {
            NEVER;
        }

        //TODO possible improvement: instead use lastfoundposidx to start search from there
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
                return(SKIP_SITE_REASON_ALLELES_CONTIG_NOMATCH);
            }
        } else {

            if (contig_changed) {
                const char* contigname = bcf_hdr_id2name(vcfd->hdr, vcfd->rec ? vcfd->rec->rid : -1);
                if (NULL == contigname) {
                    ERROR("Could not retrieve contig name for record %d.", vcfd->rec->rid);
                }
                // find the new contig, start with the last contig idx+1
                size_t cnamesidx;
                if (alleles->cnames->find_from(contigname, &cnamesidx, pars->alleles_contigidx)) {
                    pars->alleles_contigidx = cnamesidx;
                    pars->alleles_posidx = alleles->cposidx->d[cnamesidx];
                    pos_search_start = pars->alleles_posidx; // start from the first pos of the contig instead
                } else {
                    return(SKIP_SITE_REASON_ALLELES_CONTIG_NOMATCH);
                }
            }
        }

        // if not the last contig, only search until the start of next contig
        size_t pos_search_end = (pars->alleles_contigidx < (int)alleles->cnames->len - 1) ? alleles->cposidx->d[pars->alleles_contigidx + 1] : alleles->pos->len;

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
            return(SKIP_SITE_REASON_ALLELES_SITE_NOMATCH);
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
            return(SKIP_SITE_REASON_ALLELES_REF_BASE_NOMATCH);
        }
        if (pars->a1a2[1] == -1) {
            return(SKIP_SITE_REASON_ALLELES_ALT_BASE_NOMATCH);
        }

        ASSERT(pars->a1a2[0] != pars->a1a2[1]);


    } else if (args->doMajorMinor == ARG_DOMAJORMINOR_BCF_AC) {
        // -> set major/minor based on allele counts 

        an_n = bcf_get_info_int32(vcfd->hdr, vcfd->rec, "AN", &an_arr, &an_ndst);

        if (an_n <= 0) {
            ERROR("INFO/AN tag is needed for -doMajorMinor %d", args->doMajorMinor);
        } else if (an_n == 1) {
            an_val = an_arr[0];
        } else {
            ERROR("Multiple AN values found in the INFO/AN tag. Only one value is expected.");
        }

        int ref_ac_val = an_val; // init

        // -> apply filters
        ac_n = bcf_get_info_int32(vcfd->hdr, vcfd->rec, "AC", &ac_arr, &ac_ndst);
        if (ac_n <= 0) {
            ERROR("INFO/AC tag is needed for -doMajorMinor %d", args->doMajorMinor);
        }


        // determine major/minor based on ac arr 
        int max_ac_idx = 0;
        for (int i = 0; i < ac_n; ++i) {
            ref_ac_val -= ac_arr[i];
            if (ac_arr[max_ac_idx] < ac_arr[i]) {
                max_ac_idx = i;
            }
        }

        if (ref_ac_val >= ac_arr[max_ac_idx]) {
            pars->a1a2[0] = 0;
            pars->a1a2[1] = max_ac_idx + 1; // +1 to move to indexing for all alleles including ref 
        } else {
            // ref allele is not the most frequent (major) allele 

            if (ac_n == 1) {
                // site is diallelic
                pars->a1a2[0] = 1;
                pars->a1a2[1] = 0;
            } else {
                // search for the second most frequent allele among the alts
                int max2_ac_idx = max_ac_idx == 0 ? 1 : 0;
                for (int i = 0; i < ac_n; ++i) {
                    if (i == max_ac_idx) {
                        continue;
                    }
                    if (ac_arr[max2_ac_idx] < ac_arr[i]) {
                        max2_ac_idx = i;
                    }
                }
                if (ref_ac_val >= ac_arr[max2_ac_idx]) {
                    pars->a1a2[0] = max_ac_idx + 1; // +1 to move to indexing for all alleles including ref
                    pars->a1a2[1] = 0;
                } else {
                    // ref allele is not the second most frequent allele either
                    // so we need to set both to the alts
                    pars->a1a2[0] = max_ac_idx + 1; // +1 to move to indexing for all alleles including ref
                    pars->a1a2[1] = max2_ac_idx + 1; // +1 to move to indexing for all alleles including ref
                }
            }
        }


        ASSERT(pars->a1a2[0] != pars->a1a2[1]);
        ASSERT(pars->a1a2[0] >= 0);
        ASSERT(pars->a1a2[1] >= 0);
        ASSERT(pars->a1a2[0] < nAlleles);
        ASSERT(pars->a1a2[1] < nAlleles);

        DEVASSERT(pars->a1a2[0] != pars->a1a2[1]);

    }

    // ------------------------------------------------
    // FILTER: --min-a2c
    if (0 != args->min_a2c) {

        if (an_arr == NULL) {
            an_n = bcf_get_info_int32(vcfd->hdr, vcfd->rec, "AN", &an_arr, &an_ndst);
            if (an_n <= 0) {
                ERROR("INFO/AN tag is needed for --min-a2c %d", args->min_a2c);
            } else if (an_n == 1) {
                an_val = an_arr[0];
            } else {
                ERROR("Multiple AN values found in the INFO/AN tag. Only one value is expected.");
            }
        }

        if (ac_arr == NULL) {
            ac_n = bcf_get_info_int32(vcfd->hdr, vcfd->rec, "AC", &ac_arr, &ac_ndst);
            if (ac_n <= 0) {
                ERROR("INFO/AC tag is needed for --min-a2c %d", args->min_a2c);
            }
        }

        // calculate a2count based on pars->a1a2[1]
        int a2count = 0;
        if (pars->a1a2[1] == 0) {
            a2count = an_val;
            for (int i = 0; i < ac_n; ++i) {
                a2count -= ac_arr[i];
            }
        } else {
            ASSERT(pars->a1a2[1] <= ac_n);
            a2count = ac_arr[pars->a1a2[1] - 1];
        }


        if (a2count < args->min_a2c) {
            FREE(an_arr);
            FREE(ac_arr);
            return(SKIP_SITE_REASON_MIN_A2C);
        }
    }

    // ------------------------------------------------
    // FILTER: --rm-invar-sites
    if (args->rmInvarSites) {

        bool has_an_tag = true;
        if (an_arr == NULL) {
            an_n = bcf_get_info_int32(vcfd->hdr, vcfd->rec, "AN", &an_arr, &an_ndst);
            if (an_n <= 0) {
                has_an_tag = false;
            } else if (an_n == 1) {
                an_val = an_arr[0];
            } else {
                ERROR("Multiple AN values found in the INFO/AN tag. Only one value is expected.");
            }
        }
        // if not NULL, avoid reading AN tag again if it is already read in --min-a2c or doMajorMinor
        if (an_arr != NULL) {
            FREE(an_arr);
        }

        if (has_an_tag) {

            if (ac_arr == NULL) {
                ac_n = bcf_get_info_int32(vcfd->hdr, vcfd->rec, "AC", &ac_arr, &ac_ndst);
                if (ac_n <= 0) {
                    ERROR("INFO/AC tag is needed for --rm-invar-sites %d", args->rmInvarSites);
                }
            }
            // if not NULL, avoid reading AC tag again if it is already read in --min-a2c or doMajorMinor

            // calculate a2count based on pars->a1a2[1]
            int a2count = 0;
            if (pars->a1a2[1] == 0) {
                a2count = an_val;
                for (int i = 0; i < ac_n; ++i) {
                    a2count -= ac_arr[i];
                }
            } else {
                a2count = ac_arr[pars->a1a2[1] - 1];
            }

            FREE(ac_arr);

            if (a2count == 0) {
                return(SKIP_SITE_REASON_INVAR);
            } else if (a2count == (int)nInd * PROGRAM_PLOIDY) {
                return(SKIP_SITE_REASON_INVAR);
            }

        } else {
            // try AD tag

            int32_t* ads = NULL;
            int32_t n_ads = 0;
            int32_t n = 0;
            n = bcf_get_format_int32(vcfd->hdr, vcfd->rec, "AD", &ads, &n_ads);
            if (n <= 0) {
                ERROR("Program could not find INFO/AN or FORMAT/AD tags to use for --rm-invar-sites %d", args->rmInvarSites);
            }

            // make sure that AD is non-0 for at least 2 alleles 
            int n_non0[10] = { 0 };
            for (size_t i = 0; i < nInd; ++i) {
                for (size_t j = 0; j < (size_t)nAlleles; ++j) {
                    if (ads[i * nAlleles + j] > 0) {
                        n_non0[j]++;
                    }
                }
            }
            int tot_non0 = 0;
            for (size_t i = 0; i < (size_t)nAlleles; ++i) {
                if (n_non0[i] > 0) {
                    tot_non0++;
                }
            }
            FREE(ads);
            if (tot_non0 < 2) {
                return(SKIP_SITE_REASON_INVAR);
            }

        }

    }

    if (NULL != an_arr) {
        FREE(an_arr);
    }
    if (NULL != ac_arr) {
        FREE(ac_arr);
    }
    return(SKIP_SITE_REASON_NOSKIP);
}


// return 1: skip site for all individuals
//TODO check the use of block_i here, why is it not needed?
//static skip_site_reason site_read_3GL(vcfData* vcfd, paramStruct* pars, const size_t block_i, const bool contig_changed, const size_t currcontig_posidx, ibds_t* ibds, gldata_t* gldata) {
static skip_site_reason site_read_3GL(vcfData* vcfd, paramStruct* pars, const bool contig_changed, const size_t currcontig_posidx, ibds_t* ibds, gldata_t* gldata) {

    static const bool use_errorprop = (args->ibd_errorprop != 0.0);

    float* data_info_af = NULL;
    int r_info_af = 0, n_info_af = 0;

    float* lgls = NULL;
    int size_e = 0;
    int n_gls = 0;
    const size_t nInd = pars->nInd;
    float* sample_lgls = NULL;
    double fa2 = 0.0;

    bool* skipInd = (bool*)malloc(nInd * sizeof(bool));
    ASSERT(skipInd != NULL);
    // N.B. ibdseq accepts missing gts and sets the score for the pair with at least one individuals with missing gts to 0.0
    static const bool accept_missing_gts = (args->doIbd);

    int n_missing_ind = 0;
    int n_values;
    const size_t site_idx = pars->nSites;

    skip_site_reason skip_site = SKIP_SITE_REASON_NOSKIP; // init
    skip_site = read_site_a1a2(vcfd, pars, contig_changed);
    const int recd_a1_idx = pars->a1a2[0]; // filled based on alleles file if exists in read_site_a1a2
    const int recd_a2_idx = pars->a1a2[1]; // filled based on alleles file if exists in read_site_a1a2
    if (SITE_WILL_BE_SKIPPED(skip_site)) {
        goto exit;
    }

    // ------------------------------------------------
    // FILTER: --min-a2f
    /// get data: INFO_AF
    /// require: if one of the following conditions is met:
    /// -> --min-a2f is set
    /// -> a2f needed for performing doIbd
    if ((0.0 != args->min_info_a2f) || (args->doIbd)) {

        if (NULL != args->a2freqs) {

            fa2 = args->a2freqs[pars->totSites];

        } else {

            r_info_af = bcf_get_info_float(vcfd->hdr, vcfd->rec, "AF", &data_info_af, &n_info_af);
            const bool use_af_tag_for_afs = (r_info_af <= 0) ? false : true;

            if (use_af_tag_for_afs) {

                if (0 == recd_a2_idx) {
                    // if recd_a2_idx=0 (i.e. REF), then we need to set a2f to 1-sum(data_info_af) to calculate a2f
                    fa2 = 1.0;
                    for (size_t i = 0;i < (size_t)n_info_af;++i) {
                        fa2 -= data_info_af[i];
                    }
                } else {
                    // if recd_a2_idx=x, then a2f is data_info_af[x-1] (shift indexing by -1 since info_af starts with ALT)
                    fa2 = data_info_af[recd_a2_idx - 1];
                }

            } else {

                // TRY using GT tag for calculating AFs
                int32_t* gts = NULL, size_e = 0, n_values = 0;
                n_values = bcf_get_genotypes(vcfd->hdr, vcfd->rec, &gts, &size_e);
                if (n_values <= 0) {
                    ERROR("Program requires either AF tag or GT tag for calculating AFs.");
                }
                if (n_values % nInd != 0) {
                    NEVER;
                }
                if ((n_values / nInd) != PROGRAM_PLOIDY) {
                    ERROR("Ploidy %ld not supported.", n_values / nInd);
                }

                int a2count = 0, n_missing_alleles = 0;
                int32_t* sample_gts = NULL;
                for (int indi = 0; indi < (int)nInd; ++indi) {
                    sample_gts = gts + indi * PROGRAM_PLOIDY;

                    if (bcf_gt_is_missing(sample_gts[0])) {
                        n_missing_alleles++;
                        if (bcf_gt_is_missing(sample_gts[1])) {
                            n_missing_alleles++;
                            continue;
                        } else {
                            NEVER;
                        }
                    }
                    int gt1 = bcf_gt_allele(sample_gts[0]);
                    int gt2 = bcf_gt_allele(sample_gts[1]);

                    if (gt1 == recd_a1_idx) {
                        // first base in genotype is ref allele

                        if (gt2 == recd_a1_idx) {
                            // second base in genotype is ref allele
                        } else if (gt2 == recd_a2_idx) {
                            a2count++;
                        }
                    } else if (gt1 == recd_a2_idx) {
                        a2count++;
                        if (gt2 == recd_a1_idx) {
                            // second base in genotype 2 is ref allele
                        } else if (gt2 == recd_a2_idx) {
                            // second base in genotype is alt allele
                            a2count++;
                        }
                    } else {
                        // gt1 is not one of the two alleles to use
                        if (gt2 == recd_a1_idx) {
                            // second base in genotype 2 is ref allele
                        } else if (gt2 == recd_a2_idx) {
                            // second base in genotype 2 is alt allele
                            a2count++;
                        }
                    }
                }
                fa2 = ((double)a2count) / ((double)((nInd * PROGRAM_PLOIDY) - n_missing_alleles));
                FREE(gts);
            }

            if (fa2 == 0.0 && args->doIbd) {
                ERROR("At %s:%ld, allele 2 frequency is 0.0. This is not allowed for -doIbd %d. Please use --rm-invar-sites to remove invariant sites.", vcfd->get_contig_name(), vcfd->rec->pos + 1, args->doIbd);
            }
            ASSERT(fa2 >= 0.0);

        }
        if (fa2 < args->min_info_a2f) {
            skip_site = SKIP_SITE_REASON_MIN_A2F;
            goto exit;
        }
    }

    n_values = bcf_get_format_float(vcfd->hdr, vcfd->rec, "GL", &lgls, &size_e);
    if (n_values <= 0) {
        ERROR("Could not read GL tag from the VCF file.");
    }
    if (n_values % nInd != 0) {
        NEVER;
    }
    // ngls is observed number of genotypes as they appear in gl data
    n_gls = n_values / nInd;

    if (n_gls < pars->n_gc) {
        // IO::vprint(2, "Skipping site at %s:%ld. Reason: Number of GLs is less than the number of GLs needed to perform the specified analysis (%d).\n", vcfd->get_contig_name(), vcfd->rec->pos + 1, pars->n_gc);
        // return(1);
        //TODO instead of this, determine if we allow invar etc and check beforehand
        NEVER;
    }

    uint64_t ordered_gt_fetch[3];
    ordered_gt_fetch[0] = bcf_alleles2gt(recd_a1_idx, recd_a1_idx);
    ordered_gt_fetch[1] = bcf_alleles2gt(recd_a1_idx, recd_a2_idx);
    ordered_gt_fetch[2] = bcf_alleles2gt(recd_a2_idx, recd_a2_idx);

    for (size_t i = 0; i < nInd; i++) {

        skipInd[i] = false; // init

        sample_lgls = lgls + (i * n_gls);

        // ----------------------------------------------------------------- //
        // -> check for missing data
        if (bcf_float_is_missing(sample_lgls[0])) {
            // missing data check 1
            skipInd[i] = true;

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
                skipInd[i] = true;
            }
        }

        if (skipInd[i]) {
            // if minInd 0 (the program will only use sites shared across all ividuals) 
            // skip site when you encounter any missing data
            if (0 == args->minInd) {
                skip_site = SKIP_SITE_REASON_MININD;
                goto exit;
            }

            n_missing_ind++;

            // all individuals are missing 
            if ((int)nInd == n_missing_ind) {
                if (accept_missing_gts) {
                    // if we accept missing gts, we do not skip this site
                    //skipInd[indi] = true; // already initted to true
                } else {
                    skip_site = SKIP_SITE_REASON_MININD;
                    goto exit;
                }
            }

            // skip site if minInd is defined and #non-missing inds=<nInd
            if (args->minInd != -1) {
                if (((int)nInd - n_missing_ind) < args->minInd) {
                    skip_site = SKIP_SITE_REASON_MININD;
                    goto exit;
                }
            }
        }
    }


    // ------------------------------------------------
    // -> finished all filters, no exit after this point
    if (gldata != NULL) {
        for (size_t i = 0; i < nInd; i++) {

            if (skipInd[i]) {
                continue;
            }

            float* indptr = GLDATA_GET_INDPTR_AT_SITE(gldata, i, site_idx);
            sample_lgls = lgls + (i * n_gls);

            indptr[0] = sample_lgls[ordered_gt_fetch[0]];
            indptr[1] = sample_lgls[ordered_gt_fetch[1]];
            indptr[2] = sample_lgls[ordered_gt_fetch[2]];

            // infinite sensitive version:
            if ((isinf(indptr[0]) && isinf(indptr[1]) && isinf(indptr[2]))) {
                //ERROR("All 3 gls are -inf at %s:%ld for individual %ld.", vcfd->get_contig_name(), vcfd->rec->pos + 1, i);
            } else {
                // -> normalize
                float max = indptr[0];
                if (indptr[1] > max) max = indptr[1];
                if (indptr[2] > max) max = indptr[2];
                for (size_t j = 0; j < 3; ++j) {
                    indptr[j] -= max;
                }

            }

        }
    }

    //fprintf(stdout, "Reading site %s:%ld with major allele %d (%s) and minor allele %d (%s), fa2: %f \n", vcfd->get_contig_name(), vcfd->rec->pos + 1, recd_a1_idx, vcfd->rec->d.allele[recd_a1_idx], recd_a2_idx, vcfd->rec->d.allele[recd_a2_idx], fa2);

    if (ARG_DOIBD_GL_METHOD == args->doIbd) {

        /// -> calculate IBD LOD score
        ///
        /// L(IBD=k) = \sum_{G_i,G_j} P(G_i,G_j|IBD=k) P(D|G_i) P(D|G_j)
        ///
        /// score = log10(  L(IBD=1) / L(IBD=0)  )
        ///
        ///  
        // G_j:  0     1     2
        // G_i:
        //  0    0     3     6
        //  1    1     4     7
        //  2    2     5     8
        //
        // linear idx = (G_j * 3) + G_i
        /// ibdLikes[(G_j * 3) + G_i] = P(G_i,G_j|IBD=1)
        /// nullLikes[(G_j * 3) + G_i] = P(G_i,G_j|IBD=0)

        if (fa2 > 0.5) {
            ERROR("Minor allele frequency at site %s:%ld is greater than 0.5.", vcfd->get_contig_name(), vcfd->rec->pos + 1);
        }

        // NEW METHOD: use error-added allele frequencies
        double pfa2, e;
        if (use_errorprop) {
            e = MIN(args->ibd_errorprop * fa2, args->ibd_errormax);
            pfa2 = (fa2 - e) / (1.0 - 2 * e);
        } else {
            e = 0.0;
            pfa2 = fa2;
        }
        double pfa1 = 1.0 - pfa2;

        // OLD METHOD
        //double pfa2=fa2;
        //double pfa1 = 1.0 - pfa2;

        //  -> ibd likelihoods L(IBD=1)
        //  0: pfa1 * pfa1 * pfa1;
        //  1: pfa1 * pfa1 * pfa2;
        //  2: 0.0;
        //  3: pfa1 * pfa1 * pfa2;
        //  4: pfa1 * pfa2;
        //  5: pfa1 * pfa2 * pfa2;
        //  6: 0.0;
        //  7: pfa1 * pfa2 * pfa2;
        //  8: pfa2 * pfa2 * pfa2;

#if 0
        double ibdLikes4 = pfa1 * pfa2;
        double ibdLikes5 = pfa2 * ibdLikes4;
        double ibdLikes1 = pfa1 * ibdLikes4;
        double ibdLikes3 = ibdLikes1;
        double ibdLikes7 = ibdLikes5;
        double ibdLikes0 = pfa1 * pfa1 * pfa1;
        double ibdLikes8 = pfa2 * pfa2 * pfa2;
        //double ibdLikes2 = 0.0;
        //double ibdLikes6 = 0.0;

#endif

#if 1
        double e0 = 0.0, e1 = 0.0, e2 = 0.0, e3 = 0.0, e4 = 0.0;
        if (e == args->ibd_errormax) {
            e0 = args->ibd_max_error_array[0];
            e1 = args->ibd_max_error_array[1];
            e2 = args->ibd_max_error_array[2];
            e3 = args->ibd_max_error_array[3];
            e4 = args->ibd_max_error_array[4];
        } else {
            double x = 1.0 - e;
            e0 = x * x * x * x;
            e1 = (e) * (x * x * x);
            e2 = (e * e) * (x * x);
            e3 = (e * e * e) * (x);
            e4 = (e * e * e * e);
        }

        double ibdLikes0 = e0 * (pfa1 * pfa1 * pfa1) + 2 * e1 * pfa1 * pfa1 * pfa2 + e2 * pfa1 * pfa2;
        double ibdLikes1 = e0 * pfa1 * pfa1 * pfa2 + e1 * (pfa1 * pfa2 + 2 * (pfa1 * pfa1 * pfa1)) + 3 * e2 * pfa1 * pfa2;
        double ibdLikes2 = (e1 + e2 + e3) * pfa1 * pfa2 + e2 * ((pfa1 * pfa1 * pfa1) + (pfa2 * pfa2 * pfa2));
        double ibdLikes3 = ibdLikes1;
        double ibdLikes4 = (e0 + 4 * e1 + 2 * e2) * pfa1 * pfa2 + 4 * e2 * ((pfa1 * pfa1 * pfa1) + (pfa2 * pfa2 * pfa2));
        double ibdLikes5 = e0 * pfa1 * pfa2 * pfa2 + 2 * e1 * (pfa2 * pfa2 * pfa2) + (e1 + 3 * e2 + e3) * pfa1 * pfa2 + 2 * e3 * (pfa1 * pfa1 * pfa1);
        double ibdLikes6 = ibdLikes2;
        double ibdLikes7 = ibdLikes5;
        double ibdLikes8 = e0 * (pfa2 * pfa2 * pfa2) + 2 * e1 * pfa1 * pfa2 * pfa2 + e2 * pfa1 * pfa2 + 2 * e3 * pfa1 * pfa1 * pfa2 + e4 * (pfa1 * pfa1 * pfa1);
#endif


        // -> null likelihoods L(IBD=0)
        //  0: pfa1 * pfa1 * pfa1 * pfa1;
        //  1: 2 * pfa1 * pfa1 * pfa1 * pfa2;
        //  2: pfa1 * pfa1 * pfa2 * pfa2;
        //  3: 2 * pfa1 * pfa1 * pfa1 * pfa2;
        //  4: 4 * pfa1 * pfa1 * pfa2 * pfa2;
        //  5: 2 * pfa1 * pfa2 * pfa2 * pfa2;
        //  6: pfa1 * pfa1 * pfa2 * pfa2;
        //  7: 2 * pfa1 * pfa2 * pfa2 * pfa2;
        //  8: pfa2 * pfa2 * pfa2 * pfa2;

#if 0 
        double nullLikes0 = ibdLikes0 * pfa1;
        double nullLikes1 = ibdLikes1 * pfa1 * 2;
        double nullLikes2 = ibdLikes4 * ibdLikes4;
        double nullLikes3 = ibdLikes0 * pfa2 * 2;
        double nullLikes4 = nullLikes2 * 4;
        double nullLikes5 = ibdLikes8 * pfa1 * 2;
        double nullLikes6 = nullLikes2;
        double nullLikes7 = ibdLikes8 * pfa1 * 2;
        double nullLikes8 = ibdLikes8 * pfa2;
#endif 


#if 1 
        double nullLikes0 = pfa1 * pfa1 * pfa1 * pfa1;
        double nullLikes1 = 2 * pfa1 * pfa1 * pfa1 * pfa2;
        double nullLikes2 = pfa1 * pfa1 * pfa2 * pfa2;
        double nullLikes3 = nullLikes1;
        //double nullLikes4 = nullLikes2 * 4;
        double nullLikes4 = 4 * pfa1 * pfa1 * pfa2 * pfa2;
        double nullLikes5 = 2 * pfa1 * pfa2 * pfa2 * pfa2;
        double nullLikes6 = nullLikes2;
        double nullLikes7 = 2 * pfa1 * pfa2 * pfa2 * pfa2;
        double nullLikes8 = pfa2 * pfa2 * pfa2 * pfa2;
#endif 

        float* i1_lgls;
        float* i2_lgls;
        double i1lgl1, i1lgl2, i1lgl3, i2lgl1, i2lgl2, i2lgl3;
        double gls0, gls1, gls2, gls3, gls4, gls5, gls6, gls7, gls8;
        double LIBD1, LIBD0;
        size_t pidx = 0;
        for (size_t i1 = 1;i1 < nInd;++i1) {

            if (!skipInd[i1]) {
                i1_lgls = lgls + (i1 * n_gls);
                //i1lgl1 = (double)i1_lgls[ordered_gt_fetch[0]];
                //i1lgl2 = (double)i1_lgls[ordered_gt_fetch[1]];
                //i1lgl3 = (double)i1_lgls[ordered_gt_fetch[2]];
                i1_lgls = GLDATA_GET_INDPTR_AT_SITE(gldata, i1, site_idx);
                i1lgl1 = (double)i1_lgls[0];
                i1lgl2 = (double)i1_lgls[1];
                i1lgl3 = (double)i1_lgls[2];


                //if(i1lgl1==i1lgl2 && i1lgl2==i1lgl3){
                //    // site is multiallelic 

                //    // only expected with multiallelic sites
                //    ASSERT(vcfd->rec->n_allele>2); 

                //    i1_lgls = lgls + (i1 * n_gls);

                //    int tmp_recd_a3_idx=-1;
                //    if(vcfd->rec->n_allele>=3){

                //        for(int i=0; i<vcfd->rec->n_allele; ++i){
                //            if(i==recd_a1_idx || i==recd_a2_idx){
                //                continue;
                //            }
                //            if(tmp_recd_a3_idx==-1){
                //                tmp_recd_a3_idx=i;
                //                break;
                //            }
                //        }
                //        ASSERT(tmp_recd_a3_idx!=-1);

                //        // 0/0 = combined likelihoods for 0/0, 0/2, 2/2
                //        i1lgl1 = log10(pow(10.0,(double)i1_lgls[bcf_alleles2gt(recd_a1_idx, recd_a1_idx)])+pow(10.0,(double)i1_lgls[bcf_alleles2gt(recd_a1_idx, tmp_recd_a3_idx)]) + pow(10.0,(double)i1_lgls[bcf_alleles2gt(tmp_recd_a3_idx, tmp_recd_a3_idx)]));
                //        // 0/1 = combined likelihoods for 0/1, 1/2
                //        i1lgl2 = log10(pow(10.0,(double)i1_lgls[bcf_alleles2gt(recd_a1_idx, recd_a2_idx)])+pow(10.0,(double)i1_lgls[bcf_alleles2gt(recd_a2_idx, tmp_recd_a3_idx)]));
                //        // 1/1 
                //        i1lgl3 = (double)i1_lgls[bcf_alleles2gt(recd_a2_idx, recd_a2_idx)];
                //    }

                //    int tmp_recd_a4_idx=-1;
                //    if(vcfd->rec->n_allele>=4){

                //        for(int i=0; i<vcfd->rec->n_allele; ++i){
                //            if(i==recd_a1_idx || i==recd_a2_idx || i==tmp_recd_a3_idx){
                //                continue;
                //            }
                //            if(tmp_recd_a4_idx==-1){
                //                tmp_recd_a4_idx=i;
                //                break;
                //            }
                //        }
                //        ASSERT(tmp_recd_a4_idx!=-1);

                //        // 0/0 = combined likelihoods for 0/0, 0/2, 2/2 + 0/3, 2/3, 3/3
                //        i1lgl1 += log10(pow(10.0,(double)i1_lgls[bcf_alleles2gt(recd_a1_idx, tmp_recd_a4_idx)]) + pow(10.0,(double)i1_lgls[bcf_alleles2gt(tmp_recd_a3_idx, tmp_recd_a4_idx)]) + pow(10.0,(double)i1_lgls[bcf_alleles2gt(tmp_recd_a4_idx, tmp_recd_a4_idx)]));
                //        // 0/1 = combined likelihoods for 0/1, 1/2 + 1/3
                //        i1lgl2 += log10(pow(10.0,(double)i1_lgls[bcf_alleles2gt(recd_a2_idx,tmp_recd_a4_idx)]));

                //    }

                //    if(vcfd->rec->n_allele==5){

                //        int tmp_recd_a5_idx=-1;
                //        for(int i=0; i<vcfd->rec->n_allele; ++i){
                //            if(i==recd_a1_idx || i==recd_a2_idx || i==tmp_recd_a3_idx || i==tmp_recd_a4_idx){
                //                continue;
                //            }
                //            if(tmp_recd_a5_idx==-1){
                //                tmp_recd_a5_idx=i;
                //                break;
                //            }
                //        }
                //        ASSERT(tmp_recd_a5_idx!=-1);

                //        // 0/0 = combined likelihoods for 0/0, 0/2, 2/2 + 0/3, 2/3, 3/3 + 0/4, 2/4, 3/4, 4/4 
                //        i1lgl1 += log10(pow(10.0,(double)i1_lgls[bcf_alleles2gt(recd_a1_idx, tmp_recd_a5_idx)]) + pow(10.0,(double)i1_lgls[bcf_alleles2gt(tmp_recd_a3_idx, tmp_recd_a5_idx)]) + pow(10.0,(double)i1_lgls[bcf_alleles2gt(tmp_recd_a4_idx, tmp_recd_a5_idx)]) + pow(10.0,(double)i1_lgls[bcf_alleles2gt(tmp_recd_a5_idx, tmp_recd_a5_idx)]));

                //        // 0/1 = combined likelihoods for 0/1, 1/2 + 1/3 + 1/4
                //        i1lgl2 += log10(pow(10.0,(double)i1_lgls[bcf_alleles2gt(recd_a2_idx,tmp_recd_a5_idx)]));

                //    }


                //    if(i1lgl1==i1lgl2 && i1lgl2==i1lgl3){
                //        // giving up
                //        NEVER;
                //    }

                //    // normalize 3 gls
                //    double maxgl = i1lgl1;
                //    if (i1lgl2 > maxgl) maxgl = i1lgl2;
                //    if (i1lgl3 > maxgl) maxgl = i1lgl3;
                //    i1lgl1 -= maxgl;
                //    i1lgl2 -= maxgl;
                //    i1lgl3 -= maxgl;
                //}
            }

            DEVASSERT(site_idx == currcontig_posidx);

            for (size_t i2 = 0;i2 < i1;++i2) {

                if (skipInd[i1] || skipInd[i2]) {
                    ibds->pairs_ibd_scores[pidx][currcontig_posidx] = 0.0;
                    ++pidx;
                    continue;
                }

                //i2_lgls = lgls + (i2 * n_gls);
                //i2lgl1 = (double)i2_lgls[ordered_gt_fetch[0]];
                //i2lgl2 = (double)i2_lgls[ordered_gt_fetch[1]];
                //i2lgl3 = (double)i2_lgls[ordered_gt_fetch[2]];
                i2_lgls = GLDATA_GET_INDPTR_AT_SITE(gldata, i2, site_idx);
                i2lgl1 = (double)i2_lgls[0];
                i2lgl2 = (double)i2_lgls[1];
                i2lgl3 = (double)i2_lgls[2];

                //if(i2lgl1==i2lgl2 && i2lgl2==i2lgl3){
                //    // site is multiallelic 

                //    // only expected with multiallelic sites
                //    ASSERT(vcfd->rec->n_allele>2); 

                //    i2_lgls = lgls + (i2 * n_gls);

                //    int tmp_recd_a3_idx=-1;
                //    if(vcfd->rec->n_allele>=3){

                //        for(int i=0; i<vcfd->rec->n_allele; ++i){
                //            if(i==recd_a1_idx || i==recd_a2_idx){
                //                continue;
                //            }
                //            if(tmp_recd_a3_idx==-1){
                //                tmp_recd_a3_idx=i;
                //                break;
                //            }
                //        }
                //        ASSERT(tmp_recd_a3_idx!=-1);

                //        // 0/0 = combined likelihoods for 0/0, 0/2, 2/2
                //        i2lgl1 = log10(pow(10.0,(double)i2_lgls[bcf_alleles2gt(recd_a1_idx, recd_a1_idx)])+pow(10.0,(double)i2_lgls[bcf_alleles2gt(recd_a1_idx, tmp_recd_a3_idx)]) + pow(10.0,(double)i2_lgls[bcf_alleles2gt(tmp_recd_a3_idx, tmp_recd_a3_idx)]));
                //        // 0/1 = combined likelihoods for 0/1, 1/2
                //        i2lgl2 = log10(pow(10.0,(double)i2_lgls[bcf_alleles2gt(recd_a1_idx, recd_a2_idx)])+pow(10.0,(double)i2_lgls[bcf_alleles2gt(recd_a2_idx, tmp_recd_a3_idx)]));
                //        // 1/1 
                //        i2lgl3 = (double)i2_lgls[bcf_alleles2gt(recd_a2_idx, recd_a2_idx)];
                //    }


                //    int tmp_recd_a4_idx=-1;
                //    if(vcfd->rec->n_allele>=4){

                //        for(int i=0; i<vcfd->rec->n_allele; ++i){
                //            if(i==recd_a1_idx || i==recd_a2_idx || i==tmp_recd_a3_idx){
                //                continue;
                //            }
                //            if(tmp_recd_a4_idx==-1){
                //                tmp_recd_a4_idx=i;
                //                break;
                //            }
                //        }
                //        ASSERT(tmp_recd_a4_idx!=-1);

                //        // 0/0 = combined likelihoods for 0/0, 0/2, 2/2 + 0/3, 2/3, 3/3
                //        i2lgl1 += log10(pow(10.0,(double)i2_lgls[bcf_alleles2gt(recd_a1_idx, tmp_recd_a4_idx)]) + pow(10.0,(double)i2_lgls[bcf_alleles2gt(tmp_recd_a3_idx, tmp_recd_a4_idx)]) + pow(10.0,(double)i2_lgls[bcf_alleles2gt(tmp_recd_a4_idx, tmp_recd_a4_idx)]));
                //        // 0/1 = combined likelihoods for 0/1, 1/2 + 1/3
                //        i2lgl2 += log10(pow(10.0,(double)i2_lgls[bcf_alleles2gt(recd_a2_idx,tmp_recd_a4_idx)]));

                //    }

                //    if(vcfd->rec->n_allele==5){

                //        int tmp_recd_a5_idx=-1;
                //        for(int i=0; i<vcfd->rec->n_allele; ++i){
                //            if(i==recd_a1_idx || i==recd_a2_idx || i==tmp_recd_a3_idx || i==tmp_recd_a4_idx){
                //                continue;
                //            }
                //            if(tmp_recd_a5_idx==-1){
                //                tmp_recd_a5_idx=i;
                //                break;
                //            }
                //        }
                //        ASSERT(tmp_recd_a5_idx!=-1);

                //        // 0/0 = combined likelihoods for 0/0, 0/2, 2/2 + 0/3, 2/3, 3/3 + 0/4, 2/4, 3/4, 4/4 
                //        i2lgl1 += log10(pow(10.0,(double)i2_lgls[bcf_alleles2gt(recd_a1_idx, tmp_recd_a5_idx)]) + pow(10.0,(double)i2_lgls[bcf_alleles2gt(tmp_recd_a3_idx, tmp_recd_a5_idx)]) + pow(10.0,(double)i2_lgls[bcf_alleles2gt(tmp_recd_a4_idx, tmp_recd_a5_idx)]) + pow(10.0,(double)i2_lgls[bcf_alleles2gt(tmp_recd_a5_idx, tmp_recd_a5_idx)]));

                //        // 0/1 = combined likelihoods for 0/1, 1/2 + 1/3 + 1/4
                //        i2lgl2 += log10(pow(10.0,(double)i2_lgls[bcf_alleles2gt(recd_a2_idx,tmp_recd_a5_idx)]));

                //    }

                //    if(i2lgl1==i2lgl2 && i2lgl2==i2lgl3){
                //        // giving up
                //        NEVER;
                //    }

                //    // normalize 3 gls
                //    double maxgl = i2lgl1;
                //    if (i2lgl2 > maxgl) maxgl = i2lgl2;
                //    if (i2lgl3 > maxgl) maxgl = i2lgl3;
                //    i2lgl1 -= maxgl;
                //    i2lgl2 -= maxgl;
                //    i2lgl3 -= maxgl;
                //}


                gls0 = pow(10.0, (i1lgl1 + i2lgl1));
                gls1 = pow(10.0, (i1lgl2 + i2lgl1));
                gls2 = pow(10.0, (i1lgl3 + i2lgl1));
                gls3 = pow(10.0, (i1lgl1 + i2lgl2));
                gls4 = pow(10.0, (i1lgl2 + i2lgl2));
                gls5 = pow(10.0, (i1lgl3 + i2lgl2));
                gls6 = pow(10.0, (i1lgl1 + i2lgl3));
                gls7 = pow(10.0, (i1lgl2 + i2lgl3));
                gls8 = pow(10.0, (i1lgl3 + i2lgl3));

                LIBD1 = 0.0;

                LIBD1 += gls0 * ibdLikes0;
                LIBD1 += gls1 * ibdLikes1;
#if 0
                //LIBD1 += gls2 * ibdLikes2; // == 0.0
#endif
#if 1
                LIBD1 += gls2 * ibdLikes2;
#endif
                LIBD1 += gls3 * ibdLikes3;
                LIBD1 += gls4 * ibdLikes4;
                LIBD1 += gls5 * ibdLikes5;
#if 0
                //LIBD1 += gls6 * ibdLikes6; // == 0.0
#endif
#if 1
                LIBD1 += gls6 * ibdLikes6;
#endif
                LIBD1 += gls7 * ibdLikes7;
                LIBD1 += gls8 * ibdLikes8;

                LIBD0 = 0.0;

                LIBD0 += gls0 * nullLikes0;
                LIBD0 += gls1 * nullLikes1;
                LIBD0 += gls2 * nullLikes2;
                LIBD0 += gls3 * nullLikes3;
                LIBD0 += gls4 * nullLikes4;
                LIBD0 += gls5 * nullLikes5;
                LIBD0 += gls6 * nullLikes6;
                LIBD0 += gls7 * nullLikes7;
                LIBD0 += gls8 * nullLikes8;

                ASSERT(0.0 != LIBD0);

                double r = (LIBD1 / LIBD0);

                ASSERT(r >= 0.0);

                r = (double)log10(r);
                ibds->pairs_ibd_scores[pidx][currcontig_posidx] = r;
                ++pidx;
            }
        }
    }

    goto exit;

exit:

    if (lgls != NULL) {
        FREE(lgls);
    }
    if (data_info_af != NULL) {
        FREE(data_info_af);
    }
    FREE(skipInd);
    return(skip_site);
}


// return 1: skip site for all individuals
static skip_site_reason site_read_GT(jgtmat_t* jgtmat, ibds_t* ibds, vcfData* vcfd, paramStruct* pars, const bool contig_changed, const size_t currcontig_posidx) {

    static const bool use_errorprop = (args->ibd_errorprop != 0.0);

    int32_t* gts = NULL, size_e = 0, n_values = 0;

    //float* data_info_af = NULL;
    //int r_info_af = 0, n_info_af = 0;

    double fa2 = 0.0;

    double lut_ibdScores[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
    //double hbdScores_lut[3] = { 0.0 }; // TODO not using rn, maybe not even calc or start printing hbd segs as well

    int gt1 = -1;
    int gt2 = -1;
    size_t pidx = 0;
    const int nInd = pars->nInd;
    int n_missing_ind = 0;
    int n_missing_alleles = 0;

    int32_t* sample_gts = NULL;

    bool* skipInd = (bool*)malloc(nInd * sizeof(bool));
    ASSERT(skipInd != NULL);
    bool* indIsMissing = (bool*)malloc(nInd * sizeof(bool));
    ASSERT(indIsMissing != NULL);
    int* nder = (int*)malloc(nInd * sizeof(int));
    ASSERT(nder != NULL);

    int a2count = 0;

    // N.B. ibdseq accepts missing gts and sets the score for the pair with at least one individuals with missing gts to 0.0
    static const bool accept_missing_gts = (args->doIbd == ARG_DOIBD_GT_METHOD);

    // N.B. ibdseq accepts genotypes with alleles that are not one of the two alleles to use 
    static const bool accept_nona1a2_gts = (args->doIbd == ARG_DOIBD_GT_METHOD);

    skip_site_reason skip_site = SKIP_SITE_REASON_NOSKIP; // init
    skip_site = read_site_a1a2(vcfd, pars, contig_changed);
    const int recd_a1_idx = pars->a1a2[0]; // filled based on alleles file if exists in read_site_a1a2
    const int recd_a2_idx = pars->a1a2[1]; // filled based on alleles file if exists in read_site_a1a2

    if (SITE_WILL_BE_SKIPPED(skip_site)) {
        goto exit;
    }

    n_values = bcf_get_genotypes(vcfd->hdr, vcfd->rec, &gts, &size_e);
    if (n_values <= 0) {
        ERROR("Could not read GT tag from the VCF file.");
    }
    if (n_values % nInd != 0) {
        NEVER;
    }
    if ((n_values / nInd) != PROGRAM_PLOIDY) {
        ERROR("Ploidy %d not supported.", n_values / nInd);
    }

    n_missing_ind = 0;
    n_missing_alleles = 0;
    for (int indi = 0; indi < nInd; ++indi) {

        nder[indi] = 0; // init
        skipInd[indi] = true; // init
        indIsMissing[indi] = false; // init

        sample_gts = gts + indi * PROGRAM_PLOIDY;

        if (bcf_gt_is_missing(sample_gts[0])) {
            n_missing_alleles++;
            if (bcf_gt_is_missing(sample_gts[1])) {
                n_missing_alleles++;

                indIsMissing[indi] = true;

                // minInd 0
                // only use sites shared across all individuals
                // so skip site when you first encounter nan
                if (0 == args->minInd) {
                    skip_site = SKIP_SITE_REASON_MININD;
                    goto exit;
                }

                n_missing_ind++;

                // all individuals are missing 
                if (nInd == n_missing_ind) {
                    if (accept_missing_gts) {
                        // if we accept missing gts, we do not skip this site
                        //skipInd[indi] = true; // already initted to true
                        DEVASSERT(indi == nInd - 1);
                    } else {
                        skip_site = SKIP_SITE_REASON_MININD;
                        goto exit;
                    }
                }

                // skip site if minInd is defined and #non-missing inds=<nInd
                if (args->minInd != -1) {
                    if ((nInd - n_missing_ind) < args->minInd) {
                        skip_site = SKIP_SITE_REASON_MININD;
                        goto exit;
                    }
                }

                continue;
            } else {
                NEVER;
            }
        }

        // -> indi's gt is not missing
        gt1 = bcf_gt_allele(sample_gts[0]);
        gt2 = bcf_gt_allele(sample_gts[1]);


        if (gt1 == recd_a1_idx) {
            // first base in genotype is ref allele

            if (gt2 == recd_a1_idx) {
                // second base in genotype is ref allele
            } else if (gt2 == recd_a2_idx) {
                nder[indi]++;
            } else {
                if (!accept_nona1a2_gts) {
                    IO::vprint(1, "Skipping individual %d (%s) at site %s:%ld. Reason: Second allele in genotype (%s) is not one of alleles to use (%s, %s).", indi, vcfd->hdr->samples[indi], vcfd->get_contig_name(), vcfd->rec->pos + 1, vcfd->rec->d.allele[gt2], vcfd->rec->d.allele[recd_a1_idx], vcfd->rec->d.allele[recd_a2_idx]);
                    //skipInd[indi] = true; // already initted to true
                    continue;
                }
            }

        } else if (gt1 == recd_a2_idx) {

            nder[indi]++;

            if (gt2 == recd_a1_idx) {
                // second base in genotype 2 is ref allele
            } else if (gt2 == recd_a2_idx) {
                // second base in genotype is alt allele
                nder[indi]++;
            } else {
                if (!accept_nona1a2_gts) {
                    IO::vprint(1, "Skipping individual %d (%s) at site %s:%ld. Reason: Second allele in genotype (%s) is not one of alleles to use (%s, %s).", indi, vcfd->hdr->samples[indi], vcfd->get_contig_name(), vcfd->rec->pos + 1, vcfd->rec->d.allele[gt2], vcfd->rec->d.allele[recd_a1_idx], vcfd->rec->d.allele[recd_a2_idx]);
                    //skipInd[indi] = true; // already initted to true
                    continue;
                }
            }


        } else {
            // gt1 is not one of the two alleles to use

            if (!accept_nona1a2_gts) {
                IO::vprint(1, "Skipping individual %d (%s) at site %s:%ld. Reason: First allele in genotype (%s) is not one of alleles to use (%s, %s).", indi, vcfd->hdr->samples[indi], vcfd->get_contig_name(), vcfd->rec->pos + 1, vcfd->rec->d.allele[gt1], vcfd->rec->d.allele[recd_a1_idx], vcfd->rec->d.allele[recd_a2_idx]);
                //skipInd[indi] = true; // already initted to true
                continue;
            }

            if (gt2 == recd_a1_idx) {
                // second base in genotype 2 is ref allele
            } else if (gt2 == recd_a2_idx) {
                // second base in genotype 2 is alt allele
                nder[indi]++;
            } else {
                if (!accept_nona1a2_gts) {
                    IO::vprint(1, "Skipping individual %d (%s) at site %s:%ld. Reason: Second allele in genotype (%s) is not one of alleles to use (%s, %s).", indi, vcfd->hdr->samples[indi], vcfd->get_contig_name(), vcfd->rec->pos + 1, vcfd->rec->d.allele[gt2], vcfd->rec->d.allele[recd_a1_idx], vcfd->rec->d.allele[recd_a2_idx]);
                    //skipInd[indi] = true; // already initted to true
                    continue;
                }
            }
        }

        // -> indi's gt is usable (all alleles in gt are in pars->a1a2)
        skipInd[indi] = false; // set; only if it has data + passes per ind filters

        a2count += nder[indi];

    }



    // ------------------------------------------------
    // FILTER: --min-a2f
    /// get data: INFO_AF
    /// require: if one of the following conditions is met:
    /// -> --min-a2f is set
    /// -> a2f is needed for performing doIbd
    if ((0.0 != args->min_info_a2f) || (args->doIbd)) {

        if (NULL != args->a2freqs) {

            fa2 = args->a2freqs[pars->totSites];

        } else {

            // calculate a2f like ibdseq does
            // na2 / nNonMissingAlleles
            fa2 = ((double)a2count) / ((double)((nInd * PROGRAM_PLOIDY) - n_missing_alleles));
            DEVASSERT(fa2 > 0.0);

            // -> if we calculate it from AF tag:
            //r_info_af = bcf_get_info_float(vcfd->hdr, vcfd->rec, "AF", &data_info_af, &n_info_af);
            //if (r_info_af <= 0) {
            //    ERROR("AF tag not found.");
            //}

            //if (0 == recd_a2_idx) {
            //    // if recd_a2_idx=0 (i.e. REF), then we need to set a2f to 1-sum(data_info_af) to calculate a2f
            //    fa2 = 1.0;
            //    for (size_t i = 0;i < (size_t)n_info_af;++i) {
            //        fa2 -= data_info_af[i];
            //    }
            //} else {
            //    // if recd_a2_idx=x, then a2f is data_info_af[x-1] (shift indexing by -1 since info_af starts with ALT)
            //    fa2=data_info_af[recd_a2_idx-1];
            //}

        }

        if (fa2 < args->min_info_a2f) {
            skip_site = SKIP_SITE_REASON_MIN_A2F;
            goto exit;
        }
    }


    // ---------------------------------------------------------------------------------
    // DOIBD
    if (ARG_DOIBD_GT_METHOD == args->doIbd) {

        if (fa2 > 0.5) {
            ERROR("Minor allele frequency at site %s:%ld is greater than 0.5.", vcfd->get_contig_name(), vcfd->rec->pos + 1);
        }

        // ------------------------------------------------
        // precalculate the dose pair ibd scores
        //double fa1 = 1.0 - fa2;
        double pfa2, e;
        if (use_errorprop) {
            e = MIN(args->ibd_errorprop * fa2, args->ibd_errormax);
            pfa2 = (fa2 - e) / (1.0 - 2 * e);
        } else {
            e = 0.0;
            pfa2 = fa2;
        }

        //DEVPRINT("At site %s:%ld, fa1: %f fa2: %f e: %f pfa1: %f pfa2: %f", vcfd->get_contig_name(), vcfd->rec->pos + 1, fa1, fa2, e, pfa1, pfa2);

        double score;
        for (int dose1 = 0; dose1 < 3; ++dose1) {
            for (int dose2 = dose1; dose2 < 3; ++dose2) {

                score = log10(get_ibd_likelihood(dose1, dose2, e, pfa2) / get_null_likelihood(dose1, dose2, fa2));
                lut_ibdScores[dose1][dose2] = score;
                lut_ibdScores[dose2][dose1] = score;

                //fprintf(stdout, "dose1: %d dose2: %d fB: %f r: %f\n", dose1, dose2, fa2, score);
                //fprintf(stdout, "dose1: %d dose2: %d fB: %f r: %f ibdlike: %f nulllike:%f \n", dose1, dose2, fa2, score, get_ibd_likelihood(dose1, dose2, e, pfa2), get_null_likelihood(dose1, dose2, fa2));
            }
        }

        // getHbdScores
        //hbdScores_lut[0] = log10((pfa1 + e * e * pfa2) / (fa1 * fa1));
        //hbdScores_lut[1] = log10(e * (1 - e) / (fa1 * fa2));
        //hbdScores_lut[2] = log10((e * e * pfa1 + pfa2) / (fa2 * fa2));

        // ------------------------------------------------

        pidx = 0;

        for (int i1 = 1; i1 < nInd; ++i1) {
            for (int i2 = 0; i2 < i1; ++i2) {

                // N.B. ibdseq accepts missing gts and sets the score for the pair with at least one individuals with missing gts to 0.0
                if (indIsMissing[i1] || indIsMissing[i2]) {
                    ibds->pairs_ibd_scores[pidx][currcontig_posidx] = 0.0;
                    //fprintf(stdout, "ismisVGT1\t%d\t%d\t%d\t%d\t%f\t%f\n",i2,i1,currcontig_posidx, vcfd->rec->pos+1,ibds->pairs_ibd_scores[pidx][currcontig_posidx],ibds->pairs_ibd_scores[pidx][currcontig_posidx]);
                    ++pidx;
                    continue;
                }
                if (skipInd[i1] || skipInd[i2]) {
                    ibds->pairs_ibd_scores[pidx][currcontig_posidx] = 0.0;
                    NEVER;
                    //fprintf(stdout, "isskipVGT1\t%d\t%d\t%d\t%d\t%f\t%f\n",i2,i1,currcontig_posidx, vcfd->rec->pos+1,ibds->pairs_ibd_scores[pidx][currcontig_posidx],ibds->pairs_ibd_scores[pidx][currcontig_posidx]);
                    ++pidx;
                    continue;
                }


                int dose1 = nder[i1];
                int dose2 = nder[i2];

                DEVASSERT(dose1 >= 0 && dose1 <= 2);
                DEVASSERT(dose2 >= 0 && dose2 <= 2);


                // N.B. ibdseq implementation: it seems like from the source code that:
                // if Not correlated (isCor == false): Always calculate the score.
                // if Correlated  (isCor == true) :
                //    - If not IBS  (notIbs == true) : Calculate the score.
                //    - If IBS  (notIbs == false) : Set the score to 0.0.
                // so notIbs is only used if isCor can be true, which is not the case in our current implementation that always assumes isCor is false
                // this implies if the data is filtered for correlated sites beforehand, the ibd analysis may give different results compared to filtering during the ibd analysis
                // this was not apparent from the publication main text or the documentation 
                //bool notIbs = ((dose1 == 0 && dose2 == 2) || (dose1 == 2 && dose2 == 0));

                //ibds->pairs_ibd_scores[pidx][currcontig_posidx] = lut_ibdScores[dose1 + dose2 + (3 * ((dose1 == 0 && dose2 == 2) || (dose1 == 2 && dose2 == 0)))];
                ibds->pairs_ibd_scores[pidx][currcontig_posidx] = lut_ibdScores[dose1][dose2];
                ++pidx;
            }
        }
        DEVASSERT((int)pidx == ((int)nInd * ((int)nInd - 1)) / 2);

    }
    // END DOIBD
    // ---------------------------------------------------------------------------------


    if (jgtmat != NULL) {
        uint64_t** jgtmat_m = jgtmat->m[0]; // original run 
        pidx = 0;
        int k;
        for (int i1 = 1; i1 < nInd; ++i1) {
            if (skipInd[i1]) {
                ++pidx;
                continue;
            }

            k = nder[i1] * 3;

            for (int i2 = 0; i2 < i1; ++i2) {
                if (skipInd[i2]) {
                    ++pidx;
                    continue;
                }

                jgtmat_m[pidx][k + nder[i2]]++;
                ++pidx;
            }

        }
    }

    goto exit;


exit:

    if (gts != NULL) {
        FREE(gts);
    }
    FREE(skipInd);
    FREE(nder);
    FREE(indIsMissing);
    return(skip_site);
}

vcfData* vcfData_init(paramStruct* pars, metadata_t* metadata) {

    vcfData* vcfd = new vcfData;

    vcfd->DO_BCF_UNPACK = require_unpack();

    vcfd->in_fp = bcf_open(args->in_vcf_fn, "r");
    if (vcfd->in_fp == NULL) {
        ERROR("Could not open bcf file: %s\n", args->in_vcf_fn);
    }

    vcfd->hdr = bcf_hdr_read(vcfd->in_fp);
    vcfd->rec = bcf_init();

    if (require_index() & IDX_CSI) {
        vcfd->idx = IO::load_bcf_csi_idx(args->in_vcf_fn);
    }

    if (require_index() & IDX_TBI) {
        vcfd->tbx = IO::load_vcf_tabix_idx(args->in_vcf_fn);
        vcfd->itr = tbx_itr_querys(vcfd->tbx, args->in_region);
        if (NULL == vcfd->itr) {
            ERROR("Could not parse region: %s. Please make sure region exists and defined in the correct format.", args->in_region);
        }
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


    const size_t in_vcf_nind = bcf_hdr_nsamples(vcfd->hdr);

    //TODO match these names mtd and vcf and etc (check doc/md files for detaisl)
    if (metadata == NULL) {
        pars->names = strArray_alloc(in_vcf_nind);
        for (size_t i = 0;i < in_vcf_nind;++i) {
            pars->names->add(vcfd->hdr->samples[i]);
        }
        DEVASSERT(pars->names->len == in_vcf_nind);
        LOG("Found %ld individuals in the VCF file", pars->names->len);
        ASSERT(pars->names->len > 0);
    } else {
        // -> set pars with vcfdata
        const size_t mtd_nind = metadata->indNames->len;

        size_t midx;

        bool badorder = false;
        uint64_t* samples_metad2vcf = (uint64_t*)malloc(mtd_nind * sizeof(uint64_t));
        ASSERT(samples_metad2vcf != NULL);

        size_t newvidx = 0;
        for (size_t vidx = 0;vidx < in_vcf_nind;++vidx) {
            if (metadata->indNames->find(vcfd->hdr->samples[vidx], &midx)) {
                samples_metad2vcf[midx] = vidx;
                if (newvidx != midx) {
                    badorder = true;
                    break;
                }
                ++newvidx;
            } else {
                LOG("Skipping individual %s in the VCF file. Reason: Individual not present in metadata file.", vcfd->hdr->samples[vidx]);
                // TODO: does this actually work? are we reading data for these inds?  addtest
            }
        }
        if (newvidx != mtd_nind) {
            ERROR("All individuals listed in metadata file must be present in the VCF file. Program could only find %ld out of %ld individuals.", newvidx, mtd_nind);
        }

        if (badorder) {
            ERROR("The order of individuals in the VCF file does not match the order of individuals in the metadata file. Please sort the individuals in the metadata file to match the order of individuals in the VCF file.");
        }

        pars->names = metadata->indNames;

        kstring_t tmp = KS_INIT;
        for (size_t i = 0;i < mtd_nind;++i) {
            ksprintf(&tmp, "%s", vcfd->hdr->samples[samples_metad2vcf[i]]);
            if (i != mtd_nind - 1) {
                kputc(',', &tmp);
            }
        }
        if (bcf_hdr_set_samples(vcfd->hdr, tmp.s, 0)) {
            ERROR("Could not set samples in bcf header.");
        }
        ks_free(&tmp);
        int nsamples = bcf_hdr_nsamples(vcfd->hdr); // sanity check
        ASSERT(nsamples == (int)metadata->indNames->len);

        if ((size_t)nsamples != in_vcf_nind) {
            LOG("Will use %d individuals (out of %ld) found in the VCF file.", nsamples, in_vcf_nind);
        } else {
            LOG("Found %d individuals in the VCF file", nsamples);
        }

        FREE(samples_metad2vcf);
    }

    pars->nInd = pars->names->len;
    pars->nIndPairs = (pars->nInd * (pars->nInd - 1)) / 2; // LTED or UTED

    vcfd->nContigs = vcfd->hdr->n[BCF_DT_CTG];
    ASSERT(vcfd->nContigs > 0);

    if ((int)pars->nInd < args->minInd) {
        ERROR("Number of individuals in the VCF file is less than the minimum number of individuals provided by -minInd (%d).", args->minInd);
    }

    if (args->minInd == (int)pars->nInd) {
        DEVPRINT("-minInd %d is equal to the number of individuals found in file: %ld. Setting -minInd to 0 (all).\n", args->minInd, pars->nInd);
        args->minInd = 0;
    }

    if (pars->nInd == 1) {
        ERROR("Only one individual found in the VCF file. At least two individuals are required to perform the analysis.");
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

    //TODO add checks for other info tags for early exit
    // if (args->rmInvarSites && PROGRAM_WILL_USE_BCF_FMT_GL) {
    //     int tag_id = bcf_hdr_id2int(vcfd->hdr, BCF_DT_ID, "AD");
    //     if (!bcf_hdr_idinfo_exists(vcfd->hdr, BCF_HL_FMT, tag_id)) {
    //         ERROR("AD tag is required to remove invariant sites when using GL values (--rm-invar-sites).");
    //     }
    // }


    // calculate the total number of sites from vcf header 
    size_t maxnSites = 0;
    size_t max_perContig_nSites = 0;
    for (int i = 0; i < vcfd->nContigs; ++i) {
        size_t perContig_nSites = vcfd->hdr->id[BCF_DT_CTG][i].val->info[0];
        if (perContig_nSites > max_perContig_nSites) {
            max_perContig_nSites = perContig_nSites;
        }
        maxnSites += perContig_nSites;
    }
    ASSERT(maxnSites > 0);
    ASSERT(max_perContig_nSites > 0);
    vcfd->max_nsites = maxnSites;
    vcfd->max_percontig_nsites = max_perContig_nSites;

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




void readSites(vcfData* vcfd, paramStruct* pars, jgtmat_t* jgtmat, bblocks_t* bblocks, ibds_t* ibds, gldata_t* gldata) {
    BEGIN_LOGSECTION_MSG("(--in-vcf) Read sites from VCF file");

    // while we read the data into memory, all we do for bootstrap is to keep track of the number of sites in each block
    static const bool doBlocks = (bblocks != NULL);
    static const size_t nBlocks = (doBlocks) ? bblocks->n_blocks : 1;
    uint64_t* nsites_per_block = NULL; // number of unskipped sites in each block
    if (doBlocks) {
        nsites_per_block = bblocks->nsites_per_block;
    }

    const char* this_contigName = NULL;

    int contig_i = -1;
    size_t block_i = 0;

    bool isFirstRecord = true; // init
    bool newBlock;

    bool contig_changed = false;

    size_t block_start_pos;
    size_t block_end_pos;
    size_t this_pos;

    size_t currcontig_posidx = 0;

    kstring_t ks_prev_contigName = KS_INIT;

    skip_site_reason skip_site = SKIP_SITE_REASON_NOSKIP; // init

    while (1 == (vcfd->records_next())) {

        skip_site = SKIP_SITE_REASON_NOSKIP; // init


        if (vcfd->rec->pos >= 0) {
            this_pos = vcfd->rec->pos;
        } else {
            NEVER;
        }

        if (gldata != NULL) {
            //TODO instead check if nSites%%gldata->step==0?
            if (pars->nSites == gldata->n_sites) {
                gldata_expand(gldata);
            } else if (pars->nSites > gldata->n_sites) {
                NEVER;
            }
        }


        if (ibds != NULL) {
            if (currcontig_posidx == ibds->size) {
                ibds_realloc(ibds);
            } else if (currcontig_posidx > ibds->size) {
                NEVER;
            }
        }


        newBlock = false;
        contig_changed = false; // -> reset

        ASSERT(NULL != (this_contigName = bcf_hdr_id2name(vcfd->hdr, vcfd->rec ? vcfd->rec->rid : -1)));


        // TODO check: what happens is the first site is skipped? does isFirstRecord and its expectations etc play nicely ?

        if (!isFirstRecord) {

            if (0 != strcmp(ks_prev_contigName.s, this_contigName)) {

                // ** contig change **
                contig_changed = true; // set

                if (args->doIbd) {
                    // snprintf(ibds->contig_name,sizeof(ibds->contig_name),"%s",ks_prev_contigName);
                    ks_clear(&ibds->ks_contig_name);
                    ksprintf(&ibds->ks_contig_name, "%s", ks_prev_contigName.s);
                }

                DEVASSERT(contig_i >= 0);

                ++contig_i;
                currcontig_posidx = 0;

                ks_clear(&ks_prev_contigName);
                ksprintf(&ks_prev_contigName, "%s", this_contigName);


                if (doBlocks) {
                    bblocks_match_contig(bblocks, this_contigName, contig_i);
                    newBlock = true;
                }
            }

        } else {

            // if at the very beginning; first ever pos
            // contig_changed = false; // already set before site reading loop at declaration
            contig_i = 0;
            newBlock = true;

            ksprintf(&ks_prev_contigName, "%s", this_contigName);

            isFirstRecord = false; // change once

            if (doBlocks) {
                bblocks_match_contig(bblocks, this_contigName, contig_i);
            }

        }


        // ------------------------------------------------- //
        if (contig_changed) {
            // -> perform per-contig data flush/analyses before reading in the next contig
            launch_percontig_jobs(ibds, pars);
        }

        // ------------------------------------------------- //
        // -> read data into memory



        // -> first, if doBlocks, find which block we are in
        if (doBlocks) {

            while (1) {
                block_start_pos = bblocks->block_start_pos[block_i];
                block_end_pos = bblocks->block_end_pos[block_i];
                ASSERT(block_start_pos <= block_end_pos);

                if (block_i < nBlocks) {
                    // -> compare this_contigName with the contig name of the block
                    // if they are different, then we need to move to the next block
                    if (0 != strcmp(this_contigName, bcf_hdr_id2name(vcfd->hdr, bblocks->block_contig[block_i]))) {
                        ++block_i;
                        newBlock=true;
                        continue;
                    } else {
                        // if they are the same, then we need to find the block that contains the position

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
                            ERROR("Found unexpected position (%ld). Please make sure your input file is sorted, and contact the developers if the issue persists.", this_pos + 1);
                        }

                    }
                } else if (block_i == nBlocks) {
                    // out of blocks but still not found
                    ERROR("Could not find the block for the current position %ld in contig %s", this_pos + 1, this_contigName);
                } else {
                    NEVER;
                }

            }

            ASSERT(bblocks->block_contig[block_i] == (size_t)contig_i);

            if (newBlock) {
                bblocks->block_start_siteidx[block_i] = pars->nSites;
                ASSERT(bblocks->block_contig[block_i] == (size_t)contig_i);
            }

        }
        // <- block bootstrapping prep


        if (doBlocks && !newBlock && 0 == nsites_per_block[block_i]) {
            //TODO add testcase forthis
            // a previous site that belongs to this block was found before but was skipped for all indvs
            // so treat it as a new block
            newBlock = true;
        }

        do {

            if (SITE_WILL_BE_SKIPPED(skip_site)) {
                break;
            }

            // -> then, read the needed bcf format/data fields 
            if (PROGRAM_WILL_USE_BCF_FMT_GL) {
                //skip_site = site_read_3GL(vcfd, pars, block_i, contig_changed, currcontig_posidx, ibds, gldata);
                skip_site = site_read_3GL(vcfd, pars, contig_changed, currcontig_posidx, ibds, gldata);
            } else if (PROGRAM_WILL_USE_BCF_FMT_GT) {
                skip_site = site_read_GT(jgtmat, ibds, vcfd, pars, contig_changed, currcontig_posidx);
            } else {
                NEVER;
            }

        } while (0);


        if (SITE_WILL_BE_SKIPPED(skip_site)) {
            IO::vprint(1, "Skipping site at %s:%ld. Reason: %s\n", vcfd->get_contig_name(), this_pos + 1, get_skip_site_reason_str(skip_site));
            pars->totSites++;
            continue;
        }

        // ------------------------------------------------- //
        // -> if site is NOT skipped

        if (doBlocks) {
            nsites_per_block[block_i]++;
        }

        if (args->doIbd) {
            ibds->pos0[currcontig_posidx] = this_pos;
        }


        currcontig_posidx++;
        pars->nSites++;
        pars->totSites++;

        // ------------------------------------------------- //

    } // end sites loop


    // -> perform per-contig data flush/analyses for the last contig
    if (args->doIbd) {
        ks_clear(&ibds->ks_contig_name);
        ksprintf(&ibds->ks_contig_name, "%s", ks_prev_contigName.s);
    }
    launch_percontig_jobs(ibds, pars);


    if (doBlocks) {
        for (size_t i = 0;i < nBlocks;++i) {
            //if (nsites_per_block[i] == 0) {
                // will include the block in bootsrap sampling but will exclude it from the bootstrap downstream analyses by checking its nsites_per_block
            //}
            if (PROGRAM_VERBOSITY_LEVEL >= 1) {
                LOG("Found %ld sites in block %ld (%s:%ld-%ld)", nsites_per_block[i], i, bblocks->contig_names->get(bblocks->block_contig[i]), bblocks->block_start_pos[i] + 1, bblocks->block_end_pos[i]);

            }
        }
    }


    if (args->bcfSrc & ARG_INTPLUS_BCFSRC_FMT_GT) {

        if (jgtmat != NULL) {

            uint64_t snSites;

            size_t pidx = 0;
            // LTED loop
            for (size_t i1 = 1;i1 < pars->nInd;++i1) {
                for (size_t i2 = 0;i2 < i1;++i2) {
                    snSites = 0;
                    for (size_t i = 0;i < jgtmat->n_gc;++i) {
                        snSites += jgtmat->m[0][pidx][i];
                    }

                    if ((int)snSites < args->pair_min_n_sites) {
                        if (args->allow_mispairs) {
                            IO::vprint(1, "Dropping pair %ld (%s, %s) due to insufficient number of shared sites (%ld)", pidx, pars->names->d[i1], pars->names->d[i2], snSites);
                            jgtmat->drop[0][pidx] = true;
                        } else {
                            ERROR("Individuals %s and %s have %ld shared sites, which is less than the minimum number of shared sites required (%d).", pars->names->d[i1], pars->names->d[i2], snSites, args->pair_min_n_sites);
                        }
                    }
                    jgtmat->snsites[0][pidx] = snSites;
                    ++pidx;
                }
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

    ks_free(&ks_prev_contigName);
    return;
}

