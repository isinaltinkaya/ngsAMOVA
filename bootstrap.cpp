#include "bootstrap.h"

void blobStruct::_print() {
    for (int i = 0; i < nBlocks; ++i) {
        fprintf(stderr, "%s\t%d\t%d\t%d\n", blocks[i]->chr, blocks[i]->start, blocks[i]->end, blocks[i]->len);
    }
}

blobStruct::~blobStruct() {
    for (int i = 0; i < nBlocks; ++i) {
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
        blockPtrs = (int **)realloc(blockPtrs, nBlocks * sizeof(int *));
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
        fprintf(stderr, "%s\t%d\t%d\n", chr, start_int, end_int);
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
        fprintf(stderr, "%s\t%d\t%d\n", chr, start_int, end_int);
    }

    ASSERT(blob->nBlocks == nBlocks);

    FREE(firstLine);
    FREE(tok);
    FCLOSE(fp);

    return blob;
}

void blobStruct::print(IO::outputStruct *out_blockstab_fs) {
    fprintf(stderr, "\n[INFO]\t-> Writing blocks tab file: %s\n", out_blockstab_fs->fn);
    kstring_t *kbuf = kbuf_init();

    // blocks tab file = 1-based, inclusive start, inclusive end
    // internal representation = 0-based, inclusive start, exclusive end
    // convert internal representation to blocks tab file representation
    for (int bi = 0; bi < nBlocks; ++bi) {
        ksprintf(kbuf, "%s\t%d\t%d\n", blocks[bi]->chr, blocks[bi]->start + 1, blocks[bi]->end);
    }
    out_blockstab_fs->write(kbuf);
    kbuf_destroy(kbuf);
}

blobStruct *blobStruct_populate_blocks_withSize(vcfData *vcf, argStruct *args) {
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

blobStruct *blobStruct_get(vcfData *vcf, paramStruct *pars, argStruct *args, distanceMatrixStruct *dMS, metadataStruct *mS, formulaStruct *formulaSt) {
    fprintf(stderr, "\n\t-> --nBootstraps %d is set, will perform %d bootstraps for AMOVA significance testing.\n", args->nBootstraps, args->nBootstraps);

    if (args->nBootstraps < 0) {
        ERROR("nBootstraps should be a positive integer or 0, but found %d.", args->nBootstraps);
    }

    blobStruct *blob = NULL;

    if (args->blockSize != 0) {
        fprintf(stderr, "\n\t-> blockSize is set, will perform block bootstrapping with blocks of size %d.\n", args->blockSize);
        blob = blobStruct_populate_blocks_withSize(vcf, args);
    } else if (args->in_blocks_tab_fn != NULL) {
        blob = blobStruct_read_tab(args->in_blocks_tab_fn);
    } else if (args->in_blocks_bed_fn != NULL) {
        blob = blobStruct_read_bed(args->in_blocks_bed_fn);
    } else {
        ERROR("--nBootstraps %d requires --block-size, --in_blocks_tab_fn or --in_blocks_bed_fn to be set.", args->nBootstraps);
    }

    bootstrapDataset_get(vcf, pars, args, dMS, mS, blob);
    ASSERT(blob != NULL);

    return blob;
}

bootstrapData *prepare_bootstrap_block_1level(vcfData *vcfd, metadataStruct *mS, blobStruct *blobSt) {
    bootstrapData *boot = (bootstrapData *)malloc(sizeof(bootstrapData));
    boot->vblocks = (int **)malloc(mS->nInd * sizeof(int *));
    for (int i = 0; i < mS->nInd; i++) {
        boot->vblocks[i] = (int *)malloc(blobSt->nBlocks * sizeof(int));
    }

    // shuffle genomic blocks among all individuals in given level
    // vblock= block from an individual at this_block
    // range [0, mS->nIndMetadata-1]

    // nLevels = 1
    //      Shuffle all blocks among all groups
    //      Sample set = all samples
    //
    //      to create an artificial individual set of size nInd

    // i : index of the artificial individual created by shuffling blocks
    for (int i = 0; i < mS->nInd; i++) {
        int vb = -1;

        for (int block = 0; block < blobSt->nBlocks; ++block) {
            // randomly sample an individual to use the block from
            // repeat sampling for each block
            // store sampled individual indices in an array of size nBlocks

            // from a sample set that is the same as the current group assignment
            vb = sample_block_variant(mS, 0, 0);
            boot->vblocks[i][block] = vb;
        }
    }
    return boot;
}

bootstrapData *prepare_bootstrap_blocks_multilevel(vcfData *vcfd, paramStruct *pars, argStruct *args, distanceMatrixStruct *dMS, metadataStruct *mS, formulaStruct *formulaSt, blobStruct *blobSt) {
    bootstrapData *boot = (bootstrapData *)malloc(sizeof(bootstrapData));

    // multiple levels; shuffle blocks among groups of samples based on current group assignment

    // boot->vblocks = (int**)malloc(mS->nInd * sizeof(int*));

    // for (int i = 0; i < mS->nInd; i++)
    // {
    // 	boot->vblocks[i] = (int*)malloc(blobSt->nBlocks * sizeof(int));
    // }

    // // more than one level; shuffle blocks within each level
    // for (int level = 0; level < mS->nLevels; level++)
    // {
    return boot;
}

// / @brief sample an index from a set of size n
/// @param n
/// @return
int sample_from_set(int n) {
    ASSERT(n > 0);
    return (int)(n * drand48());
}

/// BLOCK BOOTSTRAPPING for AMOVA
///
/// Bootstrap data blocks among samples, while keeping the order of blocks the same
/// 1. Collect pointers to data blocks
/// 2. Shuffle pointers among groups of samples based on current group assignment

//
// nLevels = 1
//      Shuffle all blocks among all groups
//      Sample set = all samples
// nLevels = 2
//      Shuffle blocks among groups of samples based on current group assignment
//      Sample set = all samples in a given group
int sample_block_variant(metadataStruct *mtd, const int lvl, const int local_group_idx) {
    int n = 0;
    if (lvl == 0 && mtd->nLevels == 1) {
        n = mtd->nInd;

    } else {
        n = mtd->nIndPerStrata[lvl][local_group_idx];
    }
    return sample_from_set(n);
}

// Ind1  [block0][block1][block2][block3]
// Ind2  [block0][block1][block2][block3]
// Ind3  [block0][block1][block2][block3]
// Ind4  [block0][block1][block2][block3]

// vblock == the index of the 'block variant' among individuals for a given block
// e.g. synthetic individual 1
// s_ind1 [0,2] [1,1] [2,3] [3,0]
//     consists of the following block variants:
//       block0: [0,2]
//           data from individual with idx 2 (Ind3) for block 0
//       block1: [1,1]
//           data from individual with idx 1 (Ind2) for block 1
//       block2: [2,3]
//           data from individual with idx 3 (Ind4) for block 2
//       block3: [3,0]
//           data from individual with idx 0 (Ind1) for block 3
bootstrapDataset::bootstrapDataset(int nBootstraps_, int nInd_, int nBlocks_) {
    nBootstraps = nBootstraps_;
    nInd = nInd_;
    nBlocks = nBlocks_;
    bdata = (bootstrapData **)malloc(nBootstraps * sizeof(bootstrapData *));
}

bootstrapDataset::~bootstrapDataset() {
    for (int b = 0; b < nBootstraps; ++b) {
        for (int i = 0; i < nInd; ++i) {
            FREE(bdata[b]->vblocks[i]);
        }
        FREE(bdata[b]->vblocks);
        FREE(bdata[b]);
    }
    FREE(bdata);
}

bootstrapDataset *bootstrapDataset_get(vcfData *vcfd, paramStruct *pars, argStruct *args, distanceMatrixStruct *dMS, metadataStruct *mS, blobStruct *blobSt) {
    bootstrapDataset *bootstraps = new bootstrapDataset(args->nBootstraps, mS->nInd, blobSt->nBlocks);

    if (mS->nLevels == 1) {
        // only one level; shuffle all individual blocks
        for (int boot = 0; boot < args->nBootstraps; ++boot) {
            bootstraps->bdata[boot] = prepare_bootstrap_block_1level(vcfd, mS, blobSt);
        }
    } else {
        // multiple levels; shuffle blocks among groups of samples based on current group assignment
        for (int boot = 0; boot < args->nBootstraps; ++boot) {
            // bootstraps->bdata[boot] = prepare_bootstrap_blocks_multilevel(vcfd, pars, args, dMS, mS, blobSt);
        }
    }

    int from_ind = -1;

    fprintf(stderr, "\n\n~~~~\n\n\n");
    for (int i = 0; i < vcfd->nInd; i++) {
        for (int b = 0; b < bootstraps->nBlocks; b++) {
            // fprintf(stderr,"ind%d",i);
            // fprintf(stderr,",");
            // fprintf(stderr,"boot%d",b);
            // fprintf(stderr,",");
            // //site 0
            // from_ind=bootstraps->bdata[0]->vblocks[i][b];
            // fprintf(stderr,"%d",from_ind);
            fprintf(stdout, "%s,%d,%d\n", vcfd->indNames[i], b, bootstraps->bdata[0]->vblocks[i][b]);
        }
    }
    fprintf(stderr, "\n\n~~~~\n\n\n");

    return bootstraps;
}
