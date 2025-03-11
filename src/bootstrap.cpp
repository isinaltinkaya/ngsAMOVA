#include "bootstrap.h"
#include <algorithm>

#include "io.h"

bblocks_t* bblocks_init(void) {
    bblocks_t* bblocks = (bblocks_t*)malloc(sizeof(bblocks_t));
    ASSERT(bblocks != NULL);
    bblocks->n_blocks = 0;
    bblocks->n_contigs = 0;
    bblocks->n_ind = 0;
    bblocks->nsites_per_block = NULL;
    bblocks->block_start_pos = NULL;
    bblocks->block_end_pos = NULL;
    bblocks->block_contig = NULL;
    bblocks->block_start_siteidx = NULL;
    bblocks->contig_names = NULL;
    return(bblocks);
}

void bblocks_destroy(bblocks_t* bblocks) {
    for (size_t i = 0;i < (size_t)args->nBootstraps;++i) {
        FREE(bblocks->rblocks[i]);
    }
    FREE(bblocks->rblocks);
    FREE(bblocks->nsites_per_block);
    FREE(bblocks->block_start_pos);
    FREE(bblocks->block_end_pos);
    FREE(bblocks->block_contig);
    FREE(bblocks->block_start_siteidx);
    strArray_destroy(bblocks->contig_names);
    FREE(bblocks);
    return;
}


void bblocks_sample_with_replacement(bblocks_t* bblocks) {
    static const size_t nReps = args->nBootstraps;
    ASSERT(bblocks!=NULL);
    static const size_t nBlocks = bblocks->n_blocks;
    bblocks->rblocks = (size_t**)malloc(nReps * sizeof(size_t*));
    ASSERT(bblocks->rblocks != NULL);
    for (size_t rep = 0;rep < (size_t)nReps;++rep) {
        bblocks->rblocks[rep] = (size_t*)malloc(nBlocks * sizeof(size_t));
        ASSERT(bblocks->rblocks[rep] != NULL);
        for (size_t block = 0;block < (size_t)nBlocks;++block) {
            bblocks->rblocks[rep][block] = (size_t) (nBlocks * drand48());
        }
    }
    return;
}

void bblocks_print_bootstrap_samples(bblocks_t* bblocks, outfile_t* outfile) {

    kstring_t* kbuf = &outfile->kbuf;
    LOG("(%s) Writing bootstrap samples to file: %s\n", "--print-bootstrap", outfile->fn);

    // N.B. printing in 1-based [start, end] format
    ksprintf(kbuf, "Rep\tPos\tBlockID\tBlockContig\tBlockStart\tBlockEnd\tNSites\n");
    for (size_t rep = 0;rep < (size_t)args->nBootstraps;++rep) {
        for (size_t block = 0;block < (size_t)bblocks->n_blocks;++block) {
            ksprintf(kbuf, "%ld\t%ld\t%ld\t%s\t%ld\t%ld\t%ld\n", rep, block, bblocks->rblocks[rep][block], bblocks->contig_names->d[bblocks->block_contig[bblocks->rblocks[rep][block]]], bblocks->block_start_pos[bblocks->rblocks[rep][block]] + 1, bblocks->block_end_pos[bblocks->rblocks[rep][block]], bblocks->nsites_per_block[bblocks->rblocks[rep][block]]);
        }
    }
    return;
}

// blocks tsv file = 1-based, [start, end]
// internal representation = 0-based, [start, end)
// conversion from internal to blocks tab: start+1,end
void bblocks_print_blocks_tab(bblocks_t* bblocks, outfile_t* outfile) {
    kstring_t* kbuf = &outfile->kbuf;

    LOG("(%s) Writing the bootstrapping blocks to tsv file: %s", "--print-blocks", outfile->fn);

    size_t ci;
    ksprintf(kbuf, "Contig\tStart\tEnd\tNSites\n");
    for (size_t bi = 0;bi < bblocks->n_blocks;++bi) {
        ci = bblocks->block_contig[bi];
        ksprintf(kbuf, "%s\t%ld\t%ld\t%ld\n", bblocks->contig_names->d[ci], bblocks->block_start_pos[bi] + 1, bblocks->block_end_pos[bi], bblocks->nsites_per_block[bi]);
    }
    return;
}

static void bblocks_generate_blocks_with_size(bblocks_t* bblocks, vcfData* vcfd, paramStruct* pars, const uint64_t blockSize) {

    BEGIN_LOGSECTION_MSG("(--block-size) Generating blocks with given block size");

    if (vcfd->idx == NULL) {
        ERROR("-doBlockBootstrap with --block-size requires a VCF file with an index. Please run 'bcftools index' on your VCF file and try again.");
    }

    LOG("(--block-size) Targeted block size: %ld", blockSize);

    const size_t nInd = pars->nInd;
    const size_t nContigs = vcfd->nContigs;

    const int targeted_blockSize = args->blockSize;
    int current_blockSize = targeted_blockSize;

    bblocks->n_contigs = nContigs;
    bblocks->contig_names = strArray_alloc(nContigs);

    uint64_t* contig_mapped = (uint64_t*)calloc(nContigs, sizeof(uint64_t));
    uint64_t* contig_unmapped = (uint64_t*)calloc(nContigs, sizeof(uint64_t));
    ASSERT(contig_mapped != NULL && contig_unmapped != NULL);

    size_t totnBlocks = 0;

    htsFile* vcfd_in_fp = bcf_open(args->in_vcf_fn, "r");
    ASSERT(vcfd_in_fp != NULL);
    bcf_hdr_t* vcfd_hdr = bcf_hdr_dup(vcfd->hdr);
    ASSERT(vcfd_hdr != NULL);
    hts_idx_t* vcfd_idx = IO::load_bcf_csi_idx(args->in_vcf_fn);
    ASSERT(vcfd_idx != NULL);
    bcf1_t* vcfd_rec=bcf_init();
    ASSERT(vcfd_rec != NULL);

    size_t tmp_nBlocks = 128;
    bblocks->block_start_pos = (size_t*)malloc(tmp_nBlocks * sizeof(size_t));
    ASSERT(bblocks->block_start_pos != NULL);
    bblocks->block_end_pos = (size_t*)malloc(tmp_nBlocks * sizeof(size_t));
    ASSERT(bblocks->block_end_pos != NULL);
    bblocks->block_contig = (size_t*)malloc(tmp_nBlocks * sizeof(size_t));
    ASSERT(bblocks->block_contig != NULL);
    for (size_t i = 0;i < tmp_nBlocks;++i) {
        bblocks->block_start_pos[i] = 0;
        bblocks->block_end_pos[i] = 0;
        bblocks->block_contig[i] = 0;
    }


    size_t* tmp = NULL;
    for (size_t ci = 0; ci < nContigs; ++ci) {
        const char* contigName = vcfd_hdr->id[BCF_DT_CTG][ci].key;

        if (hts_idx_get_stat(vcfd_idx, ci, &contig_mapped[ci], &contig_unmapped[ci]) == -1) {
            ERROR("Failed to get statistics for contig %s. Please make sure your VCF file contains sites for this contig.", contigName);
        }
        if (contig_mapped[ci] == 0 && contig_unmapped[ci] == 0) {
            NEVER;
        }
        DEVPRINT("Contig %s has %ld sites in the VCF file.", contigName, contig_mapped[ci]);

        bblocks->contig_names->add(contigName);
        const uint64_t contigSize = vcfd_hdr->id[BCF_DT_CTG][ci].val->info[0];
        LOG("Found contig %s with size %ld", contigName, contigSize);

        if ((size_t)targeted_blockSize > contigSize) {
            WARN("Size of the contig %s (%ld) is smaller than the targeted block size (%d). Setting the block size to the size of the contig.", contigName, contigSize, targeted_blockSize);
            current_blockSize = contigSize;
        } else {
            current_blockSize = targeted_blockSize;
        }

        // +1 to account for the remainder (last block)
        const int init_nBlocks_perContig = (contigSize % current_blockSize == 0) ? (contigSize / current_blockSize) : ((contigSize / current_blockSize) + 1);
        DEVASSERT(init_nBlocks_perContig > 0);

        bool has_valid_blocks = false;
        size_t valid_nBlocks_perContig = 0;
        for (size_t bi = 0; bi < (size_t)init_nBlocks_perContig; ++bi) {
            size_t start = bi * current_blockSize;
            size_t end = start + current_blockSize;
            end = end>contigSize ? contigSize : end;

            // -> check for the presence of sites in the block
            hts_itr_t* itr = bcf_itr_queryi(vcfd_idx, ci, start, end);
            if(itr == NULL) {
                ERROR("Failed to get iterator for contig %s. Please make sure your VCF file contains sites for this contig.", contigName);
            }

            // if no sites are found, skip the block
            if(bcf_itr_next(vcfd_in_fp, itr, vcfd_rec) == -1) {
                IO::vprint(1, "Skipping block %s:%ld-%ld. Reason: No sites found in the block.", contigName, start, end);
                continue;
            }
            ////print number of sites in the block
            ////to do this, disable bcf_itr_next() above since it will advance the iterator
            //int n_sites_in_block = 0;
            //while (bcf_itr_next(vcfd->in_fp, itr, rec) >= 0) {
            //    n_sites_in_block++;
            //}
            //LOG("Found %d sites in the block %s:%ld-%ld", n_sites_in_block, contigName, start, end);

            // if sites are found, set has_valid_blocks to true
            has_valid_blocks = true;
            hts_itr_destroy(itr);
            ++valid_nBlocks_perContig;

            if(valid_nBlocks_perContig > tmp_nBlocks) {
                tmp_nBlocks += 128;

                tmp=NULL;
                tmp = (size_t*)realloc(bblocks->block_start_pos, tmp_nBlocks * sizeof(size_t));
                ASSERT(tmp != NULL);
                bblocks->block_start_pos = tmp;

                tmp = NULL;
                tmp = (size_t*)realloc(bblocks->block_end_pos, tmp_nBlocks * sizeof(size_t));
                ASSERT(tmp != NULL);
                bblocks->block_end_pos = tmp;

                tmp = NULL;
                tmp = (size_t*)realloc(bblocks->block_contig, tmp_nBlocks * sizeof(size_t));
                ASSERT(tmp != NULL);
                bblocks->block_contig = tmp;
            }
            ++totnBlocks;
            size_t block_idx=totnBlocks-1;
            bblocks->block_start_pos[block_idx] = start;
            bblocks->block_end_pos[block_idx] = end;
            bblocks->block_contig[block_idx] = ci;
        }

        // if no sites are found for any block in this contig, exit with error
        if(!has_valid_blocks) {
            ERROR("No valid blocks found for contig %s. Please make sure your VCF file contains sites for this contig.", contigName);
        }

    }

    // if there is only one block, exit with error
    if (totnBlocks == 1) {
        ERROR("Cannot perform block bootstrapping with only one block. Please decrease the block size and try again.");
    }

    // -> shrink the blocks array to the number of valid blocks
    tmp = NULL;
    tmp = (size_t*)realloc(bblocks->block_start_pos, totnBlocks * sizeof(size_t));
    ASSERT(tmp != NULL);
    bblocks->block_start_pos = tmp;

    tmp = NULL;
    tmp = (size_t*)realloc(bblocks->block_end_pos, totnBlocks * sizeof(size_t));
    ASSERT(tmp != NULL);
    bblocks->block_end_pos = tmp;

    tmp = NULL;
    tmp = (size_t*)realloc(bblocks->block_contig, totnBlocks * sizeof(size_t));
    ASSERT(tmp != NULL);
    bblocks->block_contig = tmp;
    // <- shrink

    bblocks->n_blocks = totnBlocks;
    bblocks->n_ind = nInd;

    bblocks->nsites_per_block = (uint64_t*)malloc(totnBlocks * sizeof(uint64_t));
    ASSERT(bblocks->nsites_per_block != NULL);
    bblocks->block_start_siteidx = (size_t*)malloc(totnBlocks * sizeof(size_t));
    ASSERT(bblocks->block_start_siteidx != NULL);
    for (size_t i = 0;i < totnBlocks;++i) {
        bblocks->nsites_per_block[i] = 0;
        bblocks->block_start_siteidx[i] = 0;
    }

#if DEV==1
    printf("Generated %ld blocks", totnBlocks);
    //print each block
    for (size_t i = 0;i < totnBlocks;++i) {
        printf("Block %ld: %s:%ld-%ld", i, bblocks->contig_names->d[bblocks->block_contig[i]], bblocks->block_start_pos[i], bblocks->block_end_pos[i]);
    }
#endif

    // reset the iteration for vcf record so that it can be used again
    bcf_destroy(vcfd_rec);
    bcf_hdr_destroy(vcfd_hdr);
    hts_idx_destroy(vcfd_idx);
    hts_close(vcfd_in_fp);


    FREE(contig_mapped);
    FREE(contig_unmapped);

    END_LOGSECTION_MSG("(--block-size) Generating blocks with given block size");
    return;
}


// blocks tab file = 1-based, [start, end]
// internal representation = 0-based, [start, end)
static void bblocks_read_tab(bblocks_t* bblocks, const char* fn, vcfData* vcfd, paramStruct* pars) {

    LOG("(in-blocks-tab) Reading blocks tab file: %s", fn);


    FILE* fp = IO::getFile(fn, "r");
    char* firstLine = IO::readFile::getFirstLine(fp);
    int nCols = IO::inspectFile::count_nCols(firstLine, "\t");
    if (nCols != 3) {
        ERROR("Blocks tab file must have 3 columns. Found %d columns.", nCols);
    }
    FREE(firstLine);

    ASSERT(fseek(fp, 0, SEEK_SET) == 0);

    char* tok = NULL;
    char chr[256];
    int64_t start;
    int64_t end;


    size_t nContigsFound = 0;
    char lastChr[256];
    bool contigChanged = false;

    int nContigsVcf;
    const char** vcfContigs = bcf_hdr_seqnames(vcfd->hdr, &nContigsVcf);
    DEVASSERT(nContigsVcf == vcfd->nContigs);

    size_t tmp_nBlocks = 128;

    bblocks->block_start_pos = (size_t*)malloc(tmp_nBlocks * sizeof(size_t));
    ASSERT(bblocks->block_start_pos != NULL);
    bblocks->block_end_pos = (size_t*)malloc(tmp_nBlocks * sizeof(size_t));
    ASSERT(bblocks->block_end_pos != NULL);
    for (size_t i = 0;i < tmp_nBlocks;++i) {
        bblocks->block_start_pos[i] = 0;
        bblocks->block_end_pos[i] = 0;
    }

    bblocks->contig_names = strArray_init();

    size_t nBlocks = 0;
    while (EOF != fscanf(fp, "%s\t%ld\t%ld", chr, &start, &end)) {

        // tab file positions are 1-based, but bblocks_t positions are 0-based
        // so use -1 to make it 0-based
        // both tab file start and bblocks_t start are inclusive
        if (start <= 0) {
            ERROR("Start position must be greater than 0. Note: Tab files have 1-based indexing with [start:inclusive, end:inclusive].");
        }
        bblocks->block_start_pos[nBlocks] = start - 1;

        // tab file end is inclusive, but bblocks_t end is exclusive
        // +1 to make it exclusive
        // -1+1 = 0 // so no change
        if (end <= 0) {
            ERROR("End position must be greater than 0. Note: Tab files have 1-based indexing with [start:inclusive, end:inclusive].");
        }
        bblocks->block_end_pos[nBlocks] = end;

        if (end <= start) {
            ERROR("End position must be greater than start position. Found start: %ld, end: %ld at line %ld", start, end, nBlocks);
        }

        // DEVPRINT("Found block %ld with size %ld at %s:%ld-%ld", nBlocks, end - start+1, chr, start, end);

        if (nContigsFound == 0) {
            // first loop
            strcpy(lastChr, chr);
            ++nContigsFound;
            contigChanged = true;


        } else {
            if (strcmp(lastChr, chr) == 0) {
                contigChanged = false;
            } else {
                contigChanged = true;
                strcpy(lastChr, chr);
                ++nContigsFound;
            }
        }

        // if found new contig:
        // -> check if the contig in input blocks file exists in the VCF file
        if (contigChanged) {
            size_t contig_i;
            for (contig_i = 0; contig_i < (size_t)nContigsVcf; ++contig_i) {
                if (strcmp(chr, vcfContigs[contig_i]) == 0) {
                    break;
                }
                if ((int)contig_i == nContigsVcf - 1) {
                    ERROR("Contig %s in input blocks file does not exist in the VCF file.", chr);
                }
            }
            bblocks->contig_names->add(chr);
            size_t* tmp = (size_t*)realloc(bblocks->block_contig, (nBlocks + 1) * sizeof(size_t));
            ASSERT(tmp != NULL);
            bblocks->block_contig = tmp;
            bblocks->block_contig[nBlocks] = contig_i;

        }

        ++nBlocks;
        if (tmp_nBlocks == nBlocks) {
            tmp_nBlocks += 256;
            bblocks->block_start_pos = (size_t*)realloc(bblocks->block_start_pos, tmp_nBlocks * sizeof(size_t));
            ASSERT(bblocks->block_start_pos != NULL);
            bblocks->block_end_pos = (size_t*)realloc(bblocks->block_end_pos, tmp_nBlocks * sizeof(size_t));
            ASSERT(bblocks->block_end_pos != NULL);
            for (size_t i = nBlocks;i < tmp_nBlocks;++i) {
                bblocks->block_start_pos[i] = 0;
                bblocks->block_end_pos[i] = 0;
            }
        }
    }

    bblocks->n_blocks = nBlocks;
    bblocks->n_contigs = nContigsFound;
    const size_t nInd = pars->nInd;
    bblocks->n_ind = nInd;
    bblocks->nsites_per_block = (uint64_t*)malloc(nBlocks * sizeof(uint64_t));
    ASSERT(bblocks->nsites_per_block != NULL);
    for (size_t i = 0;i < nBlocks;++i) {
        bblocks->nsites_per_block[i] = 0;
    }

    FREE(tok);
    FCLOSE(fp);
    return;
}


static void bblocks_read_bed(bblocks_t* bblocks, const char* fn, vcfData* vcfd, paramStruct* pars) {

    LOG("(in-blocks-bed) Reading blocks bed file: %s", fn);

    FILE* fp = IO::getFile(fn, "r");
    char* firstLine = IO::readFile::getFirstLine(fp);
    int nCols = IO::inspectFile::count_nCols(firstLine, "\t");
    if (nCols != 3) {
        ERROR("Blocks tab file must have 3 columns. Found %d columns.", nCols);
    }
    FREE(firstLine);

    ASSERT(fseek(fp, 0, SEEK_SET) == 0);

    char* tok = NULL;
    char chr[256];
    int64_t start;
    int64_t end;


    size_t nContigsFound = 0;
    char lastChr[256];
    bool contigChanged = false;

    int nContigsVcf;
    const char** vcfContigs = bcf_hdr_seqnames(vcfd->hdr, &nContigsVcf);
    DEVASSERT(nContigsVcf == vcfd->nContigs);

    size_t tmp_nBlocks = 256;

    bblocks->block_start_pos = (size_t*)malloc(tmp_nBlocks * sizeof(size_t));
    ASSERT(bblocks->block_start_pos != NULL);
    bblocks->block_end_pos = (size_t*)malloc(tmp_nBlocks * sizeof(size_t));
    ASSERT(bblocks->block_end_pos != NULL);
    for (size_t i = 0;i < tmp_nBlocks;++i) {
        bblocks->block_start_pos[i] = 0;
        bblocks->block_end_pos[i] = 0;
    }

    bblocks->block_contig = (size_t*)malloc(sizeof(size_t));
    ASSERT(bblocks->block_contig != NULL);

    bblocks->contig_names = strArray_init();

    size_t nBlocks = 0;
    while (EOF != fscanf(fp, "%s\t%ld\t%ld", chr, &start, &end)) {

        // both bed file and bblocks_t positions are 0-based
        // both bed file start and bblocks_t start are inclusive
        // both bed file end and bblocks_t end are exclusive
        // so no change

        bblocks->block_start_pos[nBlocks] = start;

        if (end <= 0) {
            ERROR("End position must be greater than 0. Note: Tab files have 1-based indexing with [start:inclusive, end:inclusive].");
        }
        bblocks->block_end_pos[nBlocks] = end;

        if (end <= start) {
            ERROR("End position must be greater than start position. Found start: %ld, end: %ld at line %ld", start, end, nBlocks);
        }

        // DEVPRINT("Found block %ld with size %ld at %s:%ld-%ld", nBlocks, end - start, chr, start, end);

        if (nContigsFound == 0) {
            // first loop
            strcpy(lastChr, chr);
            ++nContigsFound;
            contigChanged = true;

        } else {
            if (strcmp(lastChr, chr) == 0) {
                contigChanged = false;
            } else {
                contigChanged = true;
                strcpy(lastChr, chr);
                ++nContigsFound;
            }
        }

        // if found new contig:
        // -> check if the contig in input blocks file exists in the VCF file
        if (contigChanged) {
            size_t contig_i;
            for (contig_i = 0; contig_i < (size_t)nContigsVcf; ++contig_i) {
                if (strcmp(chr, vcfContigs[contig_i]) == 0) {
                    break;
                }
                if ((int)contig_i == nContigsVcf - 1) {
                    ERROR("Contig %s in input blocks file does not exist in the VCF file.", chr);
                }
            }
            bblocks->contig_names->add(chr);
            size_t* tmp = (size_t*)realloc(bblocks->block_contig, (nBlocks + 1) * sizeof(size_t));
            ASSERT(tmp != NULL);
            bblocks->block_contig = tmp;
            bblocks->block_contig[nBlocks] = contig_i;

        }

        ++nBlocks;
        if (tmp_nBlocks == nBlocks) {
            tmp_nBlocks *= 2;
            bblocks->block_start_pos = (size_t*)realloc(bblocks->block_start_pos, tmp_nBlocks * sizeof(size_t));
            ASSERT(bblocks->block_start_pos != NULL);
            bblocks->block_end_pos = (size_t*)realloc(bblocks->block_end_pos, tmp_nBlocks * sizeof(size_t));
            ASSERT(bblocks->block_end_pos != NULL);
            for (size_t i = nBlocks;i < tmp_nBlocks;++i) {
                bblocks->block_start_pos[i] = 0;
                bblocks->block_end_pos[i] = 0;
            }
        }
    }

    bblocks->n_blocks = nBlocks;
    bblocks->n_contigs = nContigsFound;
    const size_t nInd = pars->nInd;
    bblocks->n_ind = nInd;
    bblocks->nsites_per_block = (uint64_t*)malloc(nBlocks * sizeof(uint64_t));
    ASSERT(bblocks->nsites_per_block != NULL);
    for (size_t i = 0;i < nBlocks;++i) {
        bblocks->nsites_per_block[i] = 0;
    }

    FREE(tok);
    FCLOSE(fp);
    return;
}


void bblocks_get(bblocks_t* bblocks, vcfData* vcfd, paramStruct* pars) {

    if (PROGRAM_HAS_INPUT_BLOCKS) {
        if (args->in_blocks_bed_fn != NULL) {
            bblocks_read_bed(bblocks, args->in_blocks_bed_fn, vcfd, pars);
        } else if (args->in_blocks_tab_fn != NULL) {
            bblocks_read_tab(bblocks, args->in_blocks_tab_fn, vcfd, pars);
        } else {
            NEVER;
        }
    } else if (args->blockSize > 0) {
        bblocks_generate_blocks_with_size(bblocks, vcfd, pars, args->blockSize);
    } else {
        NEVER;
    }

    return;
}