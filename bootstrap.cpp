#include "bootstrap.h"

bootstrapData *prepare_bootstrap_block_1level(vcfData *vcfd, paramStruct *pars, argStruct *args, distanceMatrixStruct *dMS, metadataStruct *mS, formulaStruct *formulaSt, blobStruct *blobSt)
{

	bootstrapData *boot = (bootstrapData *)malloc(sizeof(bootstrapData));
	boot->vblocks = (int **)malloc(mS->nInd * sizeof(int *));
	for (int i = 0; i < mS->nInd; i++)
	{
		boot->vblocks[i] = (int *)malloc(blobSt->nBlocks * sizeof(int));
	}

	// nLevels = 1
	//      Shuffle all blocks among all groups
	//      Sample set = all samples
	//
	//      to create an artificial individual set of size nInd

	// i : index of the artificial individual created by shuffling blocks
	for (int i = 0; i < mS->nInd; i++)
	{

		for (int block = 0; block < blobSt->nBlocks; ++block)
		{

			// randomly sample an individual to use the block from
			// repeat sampling for each block
			// store sampled individual indices in an array of size nBlocks
			int vb = sample_block_variant(mS, 0);
			boot->vblocks[i][block] = vb;
		}

		// int blockStart = blobSt->contigBlockStarts[contig][block];

		// shuffle genomic blocks among all individuals in given level

		// choose a random individual among the individual set in the shuffling level
		// then set the block to the chosen individual's block at the same position
		// iblock= block from an individual at this_block
		// range [0, mS->nIndMetadata-1]

		// randomly sample an individual to use the block from
		// repeat sampling for each block
		// store sampled individual indices in an array of size nBlocks

		// 	for (int i = 0; i < mS->nInd; i++)
		// 	{
		// 		// int chosen_iblock= rand() % mS->hierArr[level]->nIndPerStrata[i];
		// 	}
		// }

		// // from a sample set that is the same as the current group assignment

		// // fprintf(stderr, "\n\n\n\n######## chosen_block: %d", chosen_iblock);

		// // distanceMatrixStruct *distanceMatrixStruct_read_csv(FILE *in_dm_fp, paramStruct *pars, argStruct *args, metadataStruct *metadataSt);

		// // to access the data
		// // if this block is not the last block in the contig
		// // 		loop: from this_blockStart to next_blockStart
		// // if this block is the last block in the contig
		// // 		loop: from this_blockStart to the end of the contig
		// for (int gti = 0; gti < 3; gti++)
		// {
		// 	double x = vcfd->lngl[0][3 * chosen_iblock + gti];
		// 	fprintf(stderr, "\n\n\n\n######## x: %f", x);
		// }
	}
	return boot;
}

bootstrapData *prepare_bootstrap_blocks_multilevel(vcfData *vcfd, paramStruct *pars, argStruct *args, distanceMatrixStruct *dMS, metadataStruct *mS, formulaStruct *formulaSt, blobStruct *blobSt)
{
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
int sample_block_variant(metadataStruct *mtd, const int r_lvl_idx)
{
	if (r_lvl_idx == 0 && mtd->nLevels == 1)
	{
		return (int)(mtd->nInd * (rand() / (RAND_MAX + 1.0)));
	}
	else
	{
		// mS->hierArr[lvl]->nIndPerStrata[sti];
		// return (int) (mS->hierArr[r_lvl_idx]->nInd * (rand() / (RAND_MAX + 1.0)));
		return -1;
	}
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
bootstrapDataset::bootstrapDataset(int nBootstraps_, int nInd_, int nBlocks_)
{
	nBootstraps = nBootstraps_;
	nInd = nInd_;
	nBlocks = nBlocks_;
	bdata = (bootstrapData **)malloc(nBootstraps * sizeof(bootstrapData *));
}

bootstrapDataset::~bootstrapDataset()
{
	for (int b = 0; b < nBootstraps; ++b)
	{
		for (int i = 0; i < nInd; ++i)
		{
			FREE(bdata[b]->vblocks[i]);
		}
		FREE(bdata[b]->vblocks);
		FREE(bdata[b]);
	}
	FREE(bdata);
}

bootstrapDataset *bootstrapDataset_get(vcfData *vcfd, paramStruct *pars, argStruct *args, distanceMatrixStruct *dMS, metadataStruct *mS, formulaStruct *formulaSt, blobStruct *blobSt)
{

	bootstrapDataset *bootstraps = new bootstrapDataset(args->nBootstraps, mS->nInd, blobSt->nBlocks);

	if (mS->nLevels == 1)
	{
		// only one level; shuffle all individual blocks
		for (int boot = 0; boot < args->nBootstraps; ++boot)
		{
			bootstraps->bdata[boot] = prepare_bootstrap_block_1level(vcfd, pars, args, dMS, mS, formulaSt, blobSt);
		}
	}
	else
	{
		// multiple levels; shuffle blocks among groups of samples based on current group assignment
		for (int boot = 0; boot < args->nBootstraps; ++boot)
		{
			bootstraps->bdata[boot] = prepare_bootstrap_blocks_multilevel(vcfd, pars, args, dMS, mS, formulaSt, blobSt);
		}
	}

	return bootstraps;
}
