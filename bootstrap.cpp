#include "bootstrap.h"

#include <algorithm>

blobStruct *blobStruct_get(vcfData *vcf, paramStruct *pars) {
    fprintf(stderr, "\n\t-> --nBootstraps %d is set, will perform %d bootstraps for AMOVA significance testing.\n", args->nBootstraps, args->nBootstraps);

    if (args->nBootstraps < 0) {
        ERROR("nBootstraps should be a positive integer or 0, but found %d.", args->nBootstraps);
    }

    blobStruct *blob = NULL;

    if (args->blockSize != 0) {
        fprintf(stderr, "\n\t-> blockSize is set, will perform block bootstrapping with blocks of size %d.\n", args->blockSize);
        blob = blobStruct_populate_blocks_withSize(vcf);
    } else if (args->in_blocks_tab_fn != NULL) {
        blob = blobStruct_read_tab(args->in_blocks_tab_fn);
    } else if (args->in_blocks_bed_fn != NULL) {
        blob = blobStruct_read_bed(args->in_blocks_bed_fn);
    } else {
        ERROR("--nBootstraps %d requires --block-size, --in_blocks_tab_fn or --in_blocks_bed_fn to be set.", args->nBootstraps);
    }

    blob->bootstraps = bootstrapDataset_get(vcf, pars, blob);

    blob->print();

    blob->bootstraps->print();

    return blob;
}

bootstrapDataset *bootstrapDataset_get(vcfData *vcfd, paramStruct *pars, blobStruct *blobSt) {
    bootstrapDataset *bootstraps = new bootstrapDataset(pars, args->nBootstraps, blobSt->nBlocks);

    int rblock = -1;

    for (int rep = 0; rep < args->nBootstraps; ++rep) {
        for (int block = 0; block < blobSt->nBlocks; ++block) {
            rblock = sample_with_replacement(blobSt->nBlocks);
            bootstraps->replicates[rep]->rBlocks[block] = rblock;
        }
    }

    return bootstraps;
}

// / @brief sample an index from a set of size n
/// @param n
/// @return
int sample_with_replacement(int n) {
    ASSERT(n > 0);
    return (int)(n * drand48());
}

void blobStruct::_print() {
    for (int i = 0; i < nBlocks; ++i) {
        fprintf(stderr, "%s\t%d\t%d\t%d\n", blocks[i]->chr, blocks[i]->start, blocks[i]->end, blocks[i]->len);
    }
}

blobStruct::~blobStruct() {
    DEL(bootstraps);
    for (int i = 0; i < nBlocks; ++i) {
        FREE(blocks[i]->chr);
        FREE(blocks[i]);
        FREE(blockPtrs[i]);
    }
    FREE(blocks);
    FREE(blockPtrs);
}

void blobStruct::addBlock() {
    ++nBlocks;
    if (nBlocks > 1) {
        blocks = (blockStruct **)realloc(blocks, nBlocks * sizeof(blockStruct *));
        ASSERT(blocks != NULL);
        blockPtrs = (int **)realloc(blockPtrs, nBlocks * sizeof(int *));
        ASSERT(blockPtrs != NULL);
    } else {
        blocks = (blockStruct **)malloc(nBlocks * sizeof(blockStruct *));
        blockPtrs = (int **)malloc(nBlocks * sizeof(int *));
    }
    blocks[nBlocks - 1] = (blockStruct *)malloc(sizeof(blockStruct));
    blockPtrs[nBlocks - 1] = (int *)malloc(sizeof(int));
}

//     - 1-based
//     - [start:included, end:included]
blobStruct *blobStruct_read_tab(const char *fn) {
    FILE *fp = IO::getFile(fn, "r");
    char *firstLine = IO::readFile::getFirstLine(fp);
    int nCols = IO::inspectFile::count_nCols(firstLine, "\t");
    if (nCols != 3) {
        ERROR("Blocks tab file must have 3 columns. Found %d columns.", nCols);
    }

    ASSERT(fseek(fp, 0, SEEK_SET) == 0);
    int nBlocks = 0;

    char *tok = NULL;
    char chr[100];
    char start[100];
    char end[100];
    blobStruct *blob = new blobStruct();

    while (EOF != fscanf(fp, "%s\t%s\t%s", chr, start, end)) {
        blob->addBlock();

        ASSERTM(strIsNumeric(start), "Start position must be numeric.");
        int start_int = atoi(start);

        ASSERTM(strIsNumeric(end), "End position must be numeric.");
        int end_int = atoi(end);

        IO::validateString(chr);
        blob->blocks[nBlocks]->chr = strdup(chr);

        // tab file positions are 1-based, but blockStruct positions are 0-based
        // so -1 to make it 0-based

        // both tab file start and blockStruct start are inclusive
        blob->blocks[nBlocks]->start = start_int - 1;
        ASSERTM(start_int > 0, "Start position must be greater than 0. Note: Tab files have 1-based indexing with [start:inclusive, end:inclusive].");

        // tab file end is inclusive, but blockStruct end is exclusive
        // +1 to make it exclusive
        // -1+1 = 0
        blob->blocks[nBlocks]->end = end_int;

        ASSERTM(end_int > 0, "End position must be greater than 0. Note: Tab files have 1-based indexing with [start:inclusive, end:inclusive].");

        blob->blocks[nBlocks]->len = blob->blocks[nBlocks]->end - blob->blocks[nBlocks]->start;
        ASSERTM(blob->blocks[nBlocks]->len > 0, "Block length must be greater than 0.");

        ++nBlocks;
        // fprintf(stderr, "%s\t%d\t%d\n", chr, start_int, end_int);
    }

    ASSERT(blob->nBlocks == nBlocks);

    FREE(firstLine);
    FREE(tok);
    FCLOSE(fp);

    return blob;
}

//     - 0-based
//     - [start:included, end:excluded)
blobStruct *blobStruct_read_bed(const char *fn) {
    FILE *fp = IO::getFile(fn, "r");
    char *firstLine = IO::readFile::getFirstLine(fp);
    int nCols = IO::inspectFile::count_nCols(firstLine, "\t");
    if (nCols != 3) {
        ERROR("Blocks bed file must have 3 columns. Found %d columns.", nCols);
    }

    ASSERT(fseek(fp, 0, SEEK_SET) == 0);
    int nBlocks = 0;

    char *tok = NULL;
    char chr[100];
    char start[100];
    char end[100];
    blobStruct *blob = new blobStruct();

    while (EOF != fscanf(fp, "%s\t%s\t%s", chr, start, end)) {
        blob->addBlock();

        ASSERTM(strIsNumeric(start), "Start position must be numeric.");
        int start_int = atoi(start);

        ASSERTM(strIsNumeric(end), "End position must be numeric.");
        int end_int = atoi(end);

        IO::validateString(chr);
        blob->blocks[nBlocks]->chr = strdup(chr);

        // both bed file and blockStruct positions are 0-based

        // both tab file start and blockStruct start are inclusive
        blob->blocks[nBlocks]->start = start_int;
        ASSERTM(start_int > 0, "Start position must be greater than 0. Note: Bed files have 0-based indexing with [start:inclusive, end:exclusive).");

        // both tab file end and blockStruct end are exclusive
        blob->blocks[nBlocks]->end = end_int;
        ASSERTM(end_int > 0, "End position must be greater than 0. Note: Bed files have 0-based indexing with [start:inclusive, end:exclusive).");

        blob->blocks[nBlocks]->len = blob->blocks[nBlocks]->end - blob->blocks[nBlocks]->start;
        ASSERTM(blob->blocks[nBlocks]->len > 0, "Block length must be greater than 0.");

        ++nBlocks;
        // fprintf(stderr, "%s\t%d\t%d\n", chr, start_int, end_int);
    }

    ASSERT(blob->nBlocks == nBlocks);

    FREE(firstLine);
    FREE(tok);
    FCLOSE(fp);

    return blob;
}

void blobStruct::print() {
    if (0 == args->printBlocksTab)
        return;

    fprintf(stderr, "\n[INFO]\t-> Writing blocks tab file: %s\n", outFiles->out_blockstab_fs->fn);
    outFiles->out_blockstab_fs->kbuf = kbuf_init();

    // blocks tab file = 1-based, inclusive start, inclusive end
    // internal representation = 0-based, inclusive start, exclusive end
    // convert internal representation to blocks tab file representation
    for (int bi = 0; bi < nBlocks; ++bi) {
        ksprintf(outFiles->out_blockstab_fs->kbuf, "%s\t%d\t%d\n", blocks[bi]->chr, blocks[bi]->start + 1, blocks[bi]->end);
    }
    outFiles->out_blockstab_fs->kbuf_write();
}

blobStruct *blobStruct_populate_blocks_withSize(vcfData *vcf) {
    // TODO make it work with region and regions, we need to get ncontigs from the region filtered vcf
    // maybe just add a lazy check afterwards to skip empty blocks due to site filtering?
    if (args->in_regions_tab_fn != NULL || args->in_regions_bed_fn != NULL || args->in_region != NULL) {
        ERROR("Block definitions cannot be used with region definitions, yet.")
    }

    const int blockSize = args->blockSize;

    blobStruct *blob = new blobStruct();

    int nBlocks = 0;
    for (int ci = 0; ci < vcf->nContigs; ci++) {
        int nBlocks_perContig = 0;
        const int contigSize = vcf->hdr->id[BCF_DT_CTG][ci].val->info[0];

        if (blockSize < contigSize) {
            if (contigSize % blockSize == 0) {
                nBlocks_perContig = contigSize / blockSize;

            } else {
                // +1 to account for the remainder (last block)
                nBlocks_perContig = (contigSize / blockSize) + 1;
            }

        } else {
            ERROR("Contig \'%s\' is smaller than the given block size (%d). Please use a smaller block size or exclude this contig from the analysis.", vcf->hdr->id[BCF_DT_CTG][ci].key, blockSize);
        }

        for (int bi = 0; bi < nBlocks_perContig; bi++) {
            blob->addBlock();
            int blockStart = bi * blockSize;

            blob->blocks[bi]->start = blockStart;

            if (bi == nBlocks_perContig - 1 && ci == vcf->nContigs - 1) {
                // last block of the last contig
                // make sure it ends at the end of the contig
                // since the size of the last block might be smaller than the block size
                blob->blocks[bi]->end = contigSize;
            } else {
                blob->blocks[bi]->end = blockStart + blockSize;
            }

            blob->blocks[bi]->chr = strdup(vcf->hdr->id[BCF_DT_CTG][ci].key);
        }
        nBlocks += nBlocks_perContig;
    }
    ASSERT(nBlocks == blob->nBlocks);

    return blob;
}

void bootstrapDataset::print() {
    if (1 != DEV)
        return;

    fprintf(stderr, "\n[INFO]\t-> Writing bootstrap info file: %s\n", outFiles->out_v_bootstrapRep_fs->fn);
    outFiles->out_v_bootstrapRep_fs->kbuf = kbuf_init();

    ksprintf(outFiles->out_v_bootstrapRep_fs->kbuf, "Replicate,ReplicateBlockIndex,BlockName\n");

    for (int b = 0; b < nReplicates; ++b) {
        for (int i = 0; i < nBlocks; ++i) {
            ksprintf(outFiles->out_v_bootstrapRep_fs->kbuf, "%d,%d,%d\n", b, i, replicates[b]->rBlocks[i]);
        }
    }
    outFiles->out_v_bootstrapRep_fs->kbuf_write();
}

bootstrapDataset::bootstrapDataset(paramStruct *pars, int nBootstraps_, int nBlocks_) {
    nReplicates = nBootstraps_;
    nBlocks = nBlocks_;
    replicates = new bootstrapReplicate *[nReplicates];
    for (int b = 0; b < nReplicates; ++b) {
        replicates[b] = new bootstrapReplicate(nBlocks);
    }
}

bootstrapDataset::~bootstrapDataset() {
    DEL2D(replicates, nReplicates);
    FREE2D(phiValues, nPhiValues);
}

void bootstrapDataset::print_confidenceInterval(FILE *fp) {
    double mean = 0.0;
    double sd = 0.0;

    double ci = 0.95;
    int n_lower = floor((double)((double)(1 - ci) / 2) * nReplicates);
    int n_upper = nReplicates - n_lower - 1;

    fprintf(fp, "PhiStatistic,nReplicates,Mean,SD,CI_lower,CI_upper\n");
    for (int i = 0; i < nPhiValues; i++) {
        std::sort(phiValues[i], phiValues[i] + nReplicates);

        mean = MATH::MEAN(phiValues[i], nReplicates);
        sd = MATH::SD(phiValues[i], nReplicates);

        fprintf(fp, "%d,%d,%f,%f,%f,%f\n", i, nReplicates, mean, sd, phiValues[i][n_lower], phiValues[i][n_upper]);
    }
}

bootstrapReplicate::bootstrapReplicate(int nBlocks) {
    rBlocks = (int *)malloc(nBlocks * sizeof(int));
    for (int i = 0; i < nBlocks; ++i) {
        rBlocks[i] = -1;
    }
}

bootstrapReplicate::~bootstrapReplicate() {
    FREE(rBlocks);
    DEL1D(amova);
    DEL1D(distanceMatrix);
}
