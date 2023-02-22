#include "bootstrap.h"


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


int r_sample_ind(metadataStruct* mtd, const int r_lvl_idx)
{
    if(r_lvl_idx==0 && mtd->nLevels==1)
    {
	    return (int) (mtd->nIndMetadata * (rand() / (RAND_MAX + 1.0)));
    }
    else
    {
		// mS->hierArr[lvl]->nIndPerStrata[sti];
        // return (int) (mS->hierArr[r_lvl_idx]->nInd * (rand() / (RAND_MAX + 1.0)));
        return -1;
    }
}

/// @brief bootstrap data blocks
/// @param argStruct pointer to input arguments
/// @param pairStruct pointer to pair of samples
/// @param outFilesStruct pointer to output files
// void bootstrap()
void prepare_bootstrap_blocks(vcfData *vcfd, paramStruct *pars, argStruct *args, distanceMatrixStruct *dMS, sampleStruct *sampleSt, metadataStruct *mS, formulaStruct *formulaSt, IO::outFilesStruct *outSt, blobStruct *blobSt)
{

		// for (int sti = 0; sti < mS->hierArr[lvl]->nStrata; sti++)
		// {
		// 	n = mS->hierArr[lvl]->nIndPerStrata[sti];
		// 	for (int i1 = 0; i1 < dMS->nInd - 1; i1++)
		// 	{
		// 		for (int i2 = i1 + 1; i2 < dMS->nInd; i2++)
		// 		{
		// 			// include only pairs that are in the same strata at this hierarchical level
		// 			if (mS->pairInStrataAtLevel(i1, i2, lvl, sti) == 1)
		// 			{
		// 				// the distance is already stored in squared form, so no need to square it again
		// 				sum += dMS->M[lut_indsToIdx[i1][i2]] / n;
		// 			}
		// 		}
		// 	}
		// }

    // nLevels = 1
    //      Shuffle all blocks among all groups
    //      Sample set = all samples
    //
    //      to create an artificial individual set of size nInd
    


	// only one level; shuffle all individual blocks
	if (mS->nLevels == 1)
	{
        
        // r_iblock_dataset[nInd][nContigs][nBlocks]
        int*** r_iblock_dataset = (int***)malloc(mS->nIndMetadata * sizeof(int**));

        // i : artificial individual created by shuffling blocks
        for (int i = 0; i < mS->nIndMetadata; i++)
        {
            r_iblock_dataset[i] = (int**)malloc(blobSt->nContigs * sizeof(int*));
            for (int contig = 0; contig < (int)blobSt->nContigs; contig++)
            {
                r_iblock_dataset[i][contig] = (int*)malloc(blobSt->contigNBlocks[contig] * sizeof(int));

				for (int block = 0; block < blobSt->contigNBlocks[contig]; block++)
				{

                    int r_ind = r_sample_ind(mS, 0);
                    r_iblock_dataset[i][contig][block] = r_ind;


        
					// int blockStart = blobSt->contigBlockStarts[contig][block];

					// shuffle genomic blocks among all individuals in given level

					// choose a random individual among the individual set in the shuffling level
					// then set the block to the chosen individual's block at the same position
					// iblock= block from an individual at this_block
                    // range [0, mS->nIndMetadata-1]




                    // randomly sample an individual to use the block from
                    // repeat sampling for each block
                    // store sampled individual indices in an array of size nBlocks






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
			}
		}
	}
	else
	{

		// more than one level; shuffle blocks within each level
		for (int level = 0; level < mS->nLevels; level++)
		{

			for (int i = 0; i < mS->nIndMetadata; i++)
			{
				// int chosen_iblock= rand() % mS->hierArr[level]->nIndPerStrata[i];
			}
		}
	}
}
