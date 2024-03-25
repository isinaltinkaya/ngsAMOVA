#include "bootstrap.h"

#include <algorithm>

// TODO add test case

void bblocks_destroy(bblocks_t* bblocks) {
    for (size_t i = 0;i < args->nBootstraps;++i) {
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
    DEVASSERT(bblocks != NULL);

    const size_t nReps = args->nBootstraps;
    const int nBlocks = bblocks->n_blocks;


    for (size_t rep = 0;rep < nReps;++rep) {
        for (size_t block = 0;block < nBlocks;++block) {
            bblocks->rblocks[rep][block] = (int)(nBlocks * drand48());
        }
    }

}

void bblocks_print_bootstrap_samples(bblocks_t* bblocks) {

    kstring_t* kbuf= kbuf_init();

    ksprintf(kbuf, "Replicate,BlockNo,BlockID\n");

    for (size_t rep = 0;rep < args->nBootstraps;++rep) {
        for (size_t block = 0;block < bblocks->n_blocks;++block) {
            ksprintf(kbuf, "%ld,%ld,%ld\n", rep, block, bblocks->rblocks[rep][block]);
        }
    }

    fprintf(stdout,"%s", kbuf->s);

    kbuf_destroy(kbuf);


    return;
}



// blocks tab file = 1-based, [start, end]
// internal representation = 0-based, [start, end)
// conversion from internal to blocks tab: start+1,end
void bblocks_print_blocks_tab(bblocks_t* bblocks) {
    LOGADD("(%s) Writing blocks tab file: %s", "-printBlocksTab", outFiles->out_blockstab_fs->fn);

    outFiles->out_blockstab_fs->kbuf = kbuf_init();

    size_t ci;
    for (size_t bi = 0;bi < bblocks->n_blocks;++bi) {

        ci = bblocks->block_contig[bi];

        ksprintf(outFiles->out_blockstab_fs->kbuf, "%s\t%ld\t%ld\n", bblocks->contig_names->d[ci], bblocks->block_start_pos[bi] + 1, bblocks->block_end_pos[bi]);
    }
    outFiles->out_blockstab_fs->kbuf_write();
    return;
}

void bblocks_generate_blocks_with_size(bblocks_t* bblocks, vcfData* vcfd, paramStruct* pars, const uint64_t blockSize) {

    BEGIN_LOGSECTION_MSG("Generating blocks for block bootstrapping");

    LOGADD("(%s) Targeted block size: %d", "--block-size", blockSize);

    const size_t nInd = pars->names->len;
    const size_t nContigs = vcfd->nContigs;

    const int targeted_blockSize = args->blockSize;
    int current_blockSize = targeted_blockSize;


    bblocks->n_contigs = nContigs;
    bblocks->contig_names = strArray_alloc(nContigs);

    uint64_t prev_nBlocks;
    uint64_t totnBlocks = 0;
    for (size_t ci = 0;ci < nContigs;++ci) {

        const uint64_t contigSize = vcfd->hdr->id[BCF_DT_CTG][ci].val->info[0];
        const char* contigName = vcfd->hdr->id[BCF_DT_CTG][ci].key;

        bblocks->contig_names->add(contigName);
        LOGADD("Found contig %s with size %d", contigName, contigSize);

        if (targeted_blockSize > contigSize) {
            LOGADD("Size of the contig %s (%d) is smaller than the targeted block size (%d). Setting the block size to the size of the contig.", contigName, contigSize, targeted_blockSize);
            current_blockSize = contigSize;
        } else {
            current_blockSize = targeted_blockSize;
        }

        // +1 to account for the remainder (last block)
        const int nBlocks_perContig = (contigSize % current_blockSize == 0) ? (contigSize / current_blockSize) : ((contigSize / current_blockSize) + 1);
        DEVASSERT(nBlocks_perContig > 0);

        LOGADD("Generating %d block%c of size %d for contig %s", nBlocks_perContig, (nBlocks_perContig > 1) ? 's' : '\0', current_blockSize, contigName);

        prev_nBlocks = totnBlocks;
        totnBlocks += nBlocks_perContig;

        if (ci == 0) {
            bblocks->block_start_pos = (size_t*)malloc(totnBlocks * sizeof(size_t));
            ASSERT(bblocks->block_start_pos != NULL);
            bblocks->block_end_pos = (size_t*)malloc(totnBlocks * sizeof(size_t));
            ASSERT(bblocks->block_end_pos != NULL);
            bblocks->block_contig = (size_t*)malloc(totnBlocks * sizeof(size_t));
            ASSERT(bblocks->block_contig != NULL);
            bblocks->block_start_siteidx = (size_t*)malloc(totnBlocks * sizeof(size_t));
            ASSERT(bblocks->block_start_siteidx != NULL);
            for(size_t i = 0;i < totnBlocks;++i) {
                bblocks->block_start_pos[i] = 0;
                bblocks->block_end_pos[i] = 0;
                bblocks->block_contig[i] = 0;
                bblocks->block_start_siteidx[i] = 0;
            }

        } else {
            size_t* tmp = NULL;
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
            bblocks->block_contig =tmp;
            tmp = NULL;

            tmp= (size_t*)realloc(bblocks->block_start_siteidx, totnBlocks * sizeof(size_t));
            ASSERT(tmp != NULL);
            bblocks->block_start_siteidx = tmp;
            tmp = NULL;

            for(size_t i = prev_nBlocks;i < totnBlocks;++i) {
                bblocks->block_start_pos[i] = 0;
                bblocks->block_end_pos[i] = 0;
                bblocks->block_contig[i] = 0;
                bblocks->block_start_siteidx[i] = 0;
            }

        }

        // -> set
        for (size_t bi = 0;bi < nBlocks_perContig;++bi) {
            bblocks->block_start_pos[bi + prev_nBlocks] = bi * current_blockSize;
            bblocks->block_end_pos[bi + prev_nBlocks] = ((bi + 1) * current_blockSize);
            bblocks->block_contig[bi + prev_nBlocks] = ci;
        }

        // <- set

    }

    if (totnBlocks == 1) {
        ERROR("Cannot perform block bootstrapping with only one block. Please decrease the block size and try again.");
    }

    bblocks->n_blocks = totnBlocks;
    bblocks->n_ind = nInd;

    bblocks->rblocks = (size_t**)malloc(args->nBootstraps * sizeof(size_t*));
    ASSERT(bblocks->rblocks != NULL);
    for (size_t i = 0;i < args->nBootstraps;++i) {
        bblocks->rblocks[i] = (size_t*)malloc(totnBlocks * sizeof(size_t));
        ASSERT(bblocks->rblocks[i] != NULL);
        for (size_t j = 0;j < totnBlocks;++j) {
            bblocks->rblocks[i][j] = 0;
        }
    }

    bblocks->nsites_per_block = (uint64_t*)malloc(totnBlocks * sizeof(uint64_t));
    ASSERT(bblocks->nsites_per_block != NULL);
    for (size_t i = 0;i < totnBlocks;++i) {
        bblocks->nsites_per_block[i] = 0;
    }

    return;
}


// blocks tab file = 1-based, [start, end]
// internal representation = 0-based, [start, end)
void bblocks_read_tab(bblocks_t* bblocks, const char* fn, vcfData* vcfd, paramStruct* pars) {

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

    size_t tmp_nBlocks = 512;

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
            for (contig_i = 0; contig_i < nContigsVcf; ++contig_i) {
                if (strcmp(chr, vcfContigs[contig_i]) == 0) {
                    break;
                }
                if (contig_i == nContigsVcf - 1) {
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
    const size_t nInd = pars->names->len;
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


void bblocks_read_bed(bblocks_t* bblocks, const char* fn, vcfData* vcfd, paramStruct* pars) {

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
            for (contig_i = 0; contig_i < nContigsVcf; ++contig_i) {
                if (strcmp(chr, vcfContigs[contig_i]) == 0) {
                    break;
                }
                if (contig_i == nContigsVcf - 1) {
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
    const size_t nInd = pars->names->len;
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
    DEVASSERT(bblocks != NULL);
    DEVASSERT(vcfd != NULL);


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
