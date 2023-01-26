/*
 * ngsAMOVA
 */

#include <stdio.h>

#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

#include <inttypes.h>
#include <limits>
#include <math.h>
#include <time.h>

#include <pthread.h>

#include "main.h"
#include "shared.h"
#include "math_utils.h"
#include "vcf_utils.h"
#include "em.h"
#include "amova.h"

using size_t = decltype(sizeof(int));

void spawnThreads_pairEM_GL(argStruct *args, paramStruct *pars, DATA::pairStruct **pairSt, VCF::vcfData *VCF, IO::outFilesStruct *outSt, DATA::distanceMatrixStruct *distMatrix)
{

	pthread_t pairThreads[pars->nIndCmb];
	threadStruct *PTHREADS[pars->nIndCmb];

	for (int i = 0; i < pars->nIndCmb; i++)
	{
		PTHREADS[i] = new threadStruct(pairSt[i], VCF->lngl, outSt->out_sfs_fs, args, pars);
	}

	int nJobs_sent = 0;

	fprintf(outSt->out_sfs_fs->fp, "Method,Ind1,Ind2,A,D,G,B,E,H,C,F,I,n_em_iter,shared_nSites,Delta,Tole,Dij,Dij2\n");
	for (int pidx = 0; pidx < pars->nIndCmb; pidx++)
	{
		if (args->mThreads > 1)
		{
			if (nJobs_sent == args->mThreads)
			{
				int t = 0;
				while (nJobs_sent > 0)
				{
					t = pidx - nJobs_sent;

					if (pthread_join(pairThreads[t], NULL) != 0)
					{
						fprintf(stderr, "\n[ERROR] Problem with joining thread.\n");
						exit(1);
					}
					else
					{
						nJobs_sent--;
					}
				}
			}
			if (pthread_create(&pairThreads[pidx], NULL, t_EM_2DSFS_GL3, PTHREADS[pidx]) == 0)
			{
				nJobs_sent++;
			}
			else
			{
				fprintf(stderr, "\n[ERROR] Problem with spawning thread.\n");
				exit(1);
			}
		}
		else
		{
			pthread_t pairThread;
			if (pthread_create(&pairThread, NULL, t_EM_2DSFS_GL3, PTHREADS[pidx]) == 0)
			{
				if (pthread_join(pairThread, NULL) != 0)
				{
					fprintf(stderr, "\n[ERROR] Problem with joining thread.\n");
					exit(1);
				}
			}
			else
			{
				fprintf(stderr, "\n[ERROR] Problem with spawning thread.\n");
				exit(1);
			}
		}
	}

	// finished indPair loop
	int t = 0;
	while (nJobs_sent > 0)
	{
		t = pars->nIndCmb - nJobs_sent;
		if (pthread_join(pairThreads[t], NULL) != 0)
		{
			fprintf(stderr, "\n[ERROR] Problem with joining thread.\n");
			exit(1);
		}
		else
		{
			nJobs_sent--;
		}
	}

	for (int pidx = 0; pidx < pars->nIndCmb; pidx++)
	{

		DATA::pairStruct *pair = PTHREADS[pidx]->pair;
		ASSERT(pair->snSites > 0);

		if (args->do_square_distance == 1)
		{
			distMatrix->M[pidx] = (double)SQUARE((MATH::EST::Dij(pair->SFS)));
		}
		else
		{
			distMatrix->M[pidx] = (double)MATH::EST::Dij(pair->SFS);
		}

		IO::print::Sfs("gle", outSt->out_sfs_fs, pair, args, VCF->hdr->samples[pars->LUT_idx2inds[pidx][0]], VCF->hdr->samples[pars->LUT_idx2inds[pidx][1]]);

		delete PTHREADS[pidx];
	}
}

void DEV_spawnThreads_pairEM_GL(argStruct *args, paramStruct *pars, DATA::pairStruct **pairSt, VCF::vcfData *VCF, IO::outFilesStruct *outSt, DATA::distanceMatrixStruct *distMatrix)
{

	pthread_t pairThreads[pars->nIndCmb];
	threadStruct *PTHREADS[pars->nIndCmb];

	for (int i = 0; i < pars->nIndCmb; i++)
	{
		PTHREADS[i] = new threadStruct(pairSt[i], VCF->lngl, outSt->out_sfs_fs, args, pars);
	}

	int nJobs_sent = 0;

	fprintf(outSt->out_sfs_fs->fp, "Method,Ind1,Ind2,A,D,G,B,E,H,C,F,I,n_em_iter,shared_nSites,Delta,Tole,Dij,Dij2\n");
	for (int pidx = 0; pidx < pars->nIndCmb; pidx++)
	{
		if (args->mThreads > 1)
		{
			if (nJobs_sent == args->mThreads)
			{
				int t = 0;
				while (nJobs_sent > 0)
				{
					t = pidx - nJobs_sent;

					if (pthread_join(pairThreads[t], NULL) != 0)
					{
						fprintf(stderr, "\n[ERROR] Problem with joining thread.\n");
						exit(1);
					}
					else
					{
						nJobs_sent--;
					}
				}
			}
			if (pthread_create(&pairThreads[pidx], NULL, DEV_t_EM_2DSFS_GL3, PTHREADS[pidx]) == 0)
			{
				nJobs_sent++;
			}
			else
			{
				fprintf(stderr, "\n[ERROR] Problem with spawning thread.\n");
				exit(1);
			}
		}
		else
		{
			pthread_t pairThread;
			if (pthread_create(&pairThread, NULL, DEV_t_EM_2DSFS_GL3, PTHREADS[pidx]) == 0)
			{
				if (pthread_join(pairThread, NULL) != 0)
				{
					fprintf(stderr, "\n[ERROR] Problem with joining thread.\n");
					exit(1);
				}
			}
			else
			{
				fprintf(stderr, "\n[ERROR] Problem with spawning thread.\n");
				exit(1);
			}
		}
	}

	// finished indPair loop
	int t = 0;
	while (nJobs_sent > 0)
	{
		t = pars->nIndCmb - nJobs_sent;
		if (pthread_join(pairThreads[t], NULL) != 0)
		{
			fprintf(stderr, "\n[ERROR] Problem with joining thread.\n");
			exit(1);
		}
		else
		{
			nJobs_sent--;
		}
	}

	for (int pidx = 0; pidx < pars->nIndCmb; pidx++)
	{

		DATA::pairStruct *pair = PTHREADS[pidx]->pair;
		ASSERT(pair->snSites > 0);

		if (args->do_square_distance == 1)
		{
			distMatrix->M[pidx] = (double)SQUARE((MATH::EST::Dij(pair->SFS)));
		}
		else
		{
			distMatrix->M[pidx] = (double)MATH::EST::Dij(pair->SFS);
		}

		IO::print::Sfs("gle", outSt->out_sfs_fs, pair, args, VCF->hdr->samples[pars->LUT_idx2inds[pidx][0]], VCF->hdr->samples[pars->LUT_idx2inds[pidx][1]]);

		delete PTHREADS[pidx];
	}
}

void prepare_bootstrap_blocks(VCF::vcfData *VCF, paramStruct *pars, argStruct *args, DATA::distanceMatrixStruct *dMS, DATA::sampleStruct *sampleSt, DATA::metadataStruct *mS, DATA::formulaStruct *formulaSt, IO::outFilesStruct *outSt, DATA::blobStruct *blobSt)
{

	fprintf(stderr, "\n\n\n\n################ Bootstrap blocks");
	// Bootstrap genomic blocks among individuals 
	// while keeping the order of blocks the same (based on AMOVA levels)
	fprintf(stderr,"\nvcf: %d", VCF->nSites);


	if (mS->nLevels == 1)
	{
		// only one level; shuffle all individual blocks
		for(int i=0; i<mS->nIndMetadata; i++)
		{

			for(int contig=0; contig < (int) blobSt->nContigs; contig++){
				//loop through blocks and shuffle
				for(int block=0; block<blobSt->contigNBlocks[contig]; block++)
				{
					int blockStart = blobSt->contigBlockStarts[contig][block];

					// shuffle genomic blocks among all individuals in given level


					// choose a random individual among the individual set in the shuffling level
					// then set the block to the chosen individual's block at the same position
					// iblock= block from an individual at this_block
					int chosen_iblock =  1 + int(mS->nIndMetadata * (rand() / (RAND_MAX + 1.0)));
					fprintf(stderr, "\n\n\n\n######## chosen_block: %d", chosen_iblock);


					// to access the data
					// if this block is not the last block in the contig
					// 		loop: from this_blockStart to next_blockStart 
					// if this block is the last block in the contig
					// 		loop: from this_blockStart to the end of the contig
					for(int gti=0; gti<3; gti++){
						double x = VCF->lngl[0][3 * chosen_iblock + gti];
						fprintf(stderr, "\n\n\n\n######## x: %f", x);
					}

				}

			}
		}

	}else{
		
		// more than one level; shuffle blocks within each level
		for(int level=0; level<mS->nLevels; level++)
		{

			for(int i=0; i<mS->nIndMetadata; i++)
			{
					// int chosen_iblock= rand() % mS->hierArr[level]->nIndPerStrata[i];
			}
		}
	}
}


void DEV_prepare_dMS_orig(argStruct *args, paramStruct *pars,  DATA::distanceMatrixStruct *dMS_orig, VCF::vcfData *VCF, DATA::pairStruct **pairSt,  DATA::formulaStruct *formulaSt, IO::outFilesStruct *outSt, DATA::sampleStruct *sampleSt)
{

	switch (args->doAMOVA)
	{

		// --------------------------- doAMOVA 1 --------------------------- //
		case 1:{


			VCF->readSites_GL(args, pars, pairSt);

			DEV_spawnThreads_pairEM_GL(args, pars, pairSt, VCF, outSt, dMS_orig);

			break;

		}

		// --------------------------- doAMOVA 2 --------------------------- //
		case 2:{


			if(args->blockSize == 0){

				VCF->readSites_GT(args, pars, pairSt);

			}else{

				fprintf(stderr,"\n[ERROR]\t-> Not implemented yet. Please use -blockSize 0\n");
				exit(1);
				// DATA::blobStruct *blobSt = DATA::blobStruct_init(VCF->nContigs, args->blockSize, VCF->hdr);
				// VCF->readSites_GT(args, pars, pairSt, blobSt);
				// DATA::blobStruct_destroy(blobSt);
			}
			for (int pidx=0; pidx < pars->nIndCmb; pidx++)
			{
				int snSites = VCF->SFS_GT3[pidx][9];

				if (args->do_square_distance == 1)
				{
					dMS_orig->M[pidx] = (double)SQUARE(MATH::EST::Dij(VCF->SFS_GT3[pidx], snSites));
				}
				else
				{
					dMS_orig->M[pidx] = (double) MATH::EST::Dij(VCF->SFS_GT3[pidx], snSites);
				}
				
				IO::print::Sfs("gt", outSt->out_sfs_fs, args, VCF->SFS_GT3[pidx], snSites, VCF->hdr->samples[pars->LUT_idx2inds[pidx][0]], VCF->hdr->samples[pars->LUT_idx2inds[pidx][1]]);
			}

			break;

		}

		// --------------------------- doAMOVA 3 --------------------------- //
		case 3:{
			break;
		}

		// --------------------------- doAMOVA 0 --------------------------- //
		case 0:{
			// should never reach here; control is done above
			// but this ain't Rust, you can never be too safe with C/C++ :P
			ASSERT(0==1); 
			break;
		}


	}

	if (args->printMatrix == 1) dMS_orig->print(outSt->out_dm_fs->fp);

}

// prepare distance matrix using original data
void prepare_dMS_orig(argStruct *args, paramStruct *pars,  DATA::distanceMatrixStruct *dMS_orig, VCF::vcfData *VCF, DATA::pairStruct **pairSt,  DATA::formulaStruct *formulaSt, IO::outFilesStruct *outSt, DATA::blobStruct *blobSt, DATA::sampleStruct *sampleSt)
{


	switch (args->doAMOVA)
	{

		// --------------------------- doAMOVA 1 --------------------------- //
		case 1:{

			if(args->blockSize == 0){

				VCF->readSites_GL(args, pars, pairSt);

			}else{

				blobSt = DATA::blobStruct_init(VCF->nContigs, args->blockSize, VCF->hdr);
				VCF->readSites_GL(args, pars, pairSt, blobSt);

			}

			spawnThreads_pairEM_GL(args, pars, pairSt, VCF, outSt, dMS_orig);

			break;

		}

		// --------------------------- doAMOVA 2 --------------------------- //
		case 2:{


			if(args->blockSize == 0){

				VCF->readSites_GT(args, pars, pairSt);

			}else{

				fprintf(stderr,"\n[ERROR]\t-> Not implemented yet. Please use -blockSize 0\n");
				exit(1);
				// DATA::blobStruct *blobSt = DATA::blobStruct_init(VCF->nContigs, args->blockSize, VCF->hdr);
				// VCF->readSites_GT(args, pars, pairSt, blobSt);
				// DATA::blobStruct_destroy(blobSt);
			}
			for (int pidx=0; pidx < pars->nIndCmb; pidx++)
			{
				int snSites = VCF->SFS_GT3[pidx][9];

				if (args->do_square_distance == 1)
				{
					dMS_orig->M[pidx] = (double)SQUARE(MATH::EST::Dij(VCF->SFS_GT3[pidx], snSites));
				}
				else
				{
					dMS_orig->M[pidx] = (double) MATH::EST::Dij(VCF->SFS_GT3[pidx], snSites);
				}
				
				IO::print::Sfs("gt", outSt->out_sfs_fs, args, VCF->SFS_GT3[pidx], snSites, VCF->hdr->samples[pars->LUT_idx2inds[pidx][0]], VCF->hdr->samples[pars->LUT_idx2inds[pidx][1]]);
			}

			break;

		}

		// --------------------------- doAMOVA 3 --------------------------- //
		case 3:{
			break;
		}

		// --------------------------- doAMOVA 0 --------------------------- //
		case 0:{
			// should never reach here; control is done above
			// but this ain't Rust, you can never be too safe with C/C++ :P
			ASSERT(0==1); 
			break;
		}

		// ---------------------- doAMOVA NOT in {0,1,2,3} ---------------------- //
		default:{
			fprintf(stderr, "\n[ERROR]\t-> -doAMOVA %d not recognized\n", args->doAMOVA);
			exit(1);
			break;
		}

	}

	if (args->printMatrix == 1) dMS_orig->print(outSt->out_dm_fs->fp);

}

// --------------------------- INPUT: VCF/BCF --------------------------- //
void input_VCF(argStruct *args, paramStruct *pars, DATA::formulaStruct *formulaSt, IO::outFilesStruct *outSt){
	fprintf(stderr, "\n[INFO]\t-> Input file type: VCF/BCF\n");

	DATA::sampleStruct *sampleSt = new DATA::sampleStruct();
	VCF::vcfData* VCF = VCF::vcfData_init(args, pars, sampleSt);
	DATA::pairStruct **pairSt = new DATA::pairStruct *[pars->nIndCmb];

	for(int pidx=0; pidx < pars->nIndCmb; pidx++){
			pairSt[pidx] = new DATA::pairStruct(pars, pidx);
	}

	if(args->printDev == 1){
		DATA::distanceMatrixStruct **dMS = new DATA::distanceMatrixStruct *[1];
		dMS[0] = new DATA::distanceMatrixStruct(pars->nInd, pars->nIndCmb, args->do_square_distance);
		fprintf(stderr, "\n[INFO]\t-> -printDev 1; will print per EM iteration distance matrix\n");
		DEV_prepare_dMS_orig(args, pars, dMS[0], VCF, pairSt, formulaSt, outSt, sampleSt);
		

		fprintf(stderr, "Total number of sites processed: %lu\n", pars->totSites);
		fprintf(stderr, "Total number of sites skipped for all individual pairs: %lu\n", pars->totSites - pars->nSites);

		outSt->flushAll();

		for (int i = 0; i < pars->nIndCmb; i++)
		{
			delete pairSt[i];
		}
		delete[] pairSt;
		delete formulaSt;
		delete sampleSt;
		delete dMS[0];
		delete[] dMS;
		

		VCF::vcfData_destroy(VCF);

		return;
	}


	if(args->doAMOVA == 0){
		// --------------------------- doAMOVA 0 --------------------------- //
		// DO NOT run AMOVA
		fprintf(stderr, "\n[INFO]\t-> -doAMOVA 0; will not run AMOVA\n");
		fprintf(stderr, "\n\t-> Skipped running AMOVA\n");
		pars->nAmovaRuns = 0;

	}else{

		// --------------------------- doAMOVA!=0 --------------------------- //
		// will run AMOVA
		// 
		// prepare for AMOVA

		// [BEGIN] ----------------------- READ METADATA ------------------------- //
		// If VCF/BCF input, read metadata after VCF
		// so you can compare VCF records with metadata when reading metadata
		ASSERT(args->in_mtd_fn!=NULL);

		FILE *in_mtd_fp = IO::getFile(args->in_mtd_fn, "r");
		DATA::metadataStruct *metadataSt = DATA::metadataStruct_get(in_mtd_fp, sampleSt, formulaSt, args->hasColNames, pars);
		FCLOSE(in_mtd_fp);
		// [END] ----------------------- READ METADATA ------------------------- //


		pars->nAmovaRuns = 1;
		if(args->nBootstraps > 0) pars->nAmovaRuns += args->nBootstraps;
		
		DATA::blobStruct *blobSt = NULL;

		DATA::distanceMatrixStruct **dMS = new DATA::distanceMatrixStruct *[pars->nAmovaRuns];

		for (int r=0; r < pars->nAmovaRuns; r++){
			dMS[r] = new DATA::distanceMatrixStruct(pars->nInd, pars->nIndCmb, args->do_square_distance);
		}

		if(args->blockSize != 0){
			blobSt = DATA::blobStruct_init(VCF->nContigs, args->blockSize, VCF->hdr);
		}


		// prepare distance matrix dMS given analysis type (-doAMOVA)

		// dMS[0] is the original distance matrix (not bootstrapped)
		prepare_dMS_orig(args, pars, dMS[0], VCF, pairSt, formulaSt, outSt, blobSt, sampleSt);


		if(args->nBootstraps > 0){

			// dMS[0] is already filled with the original distance matrix
			// bootstrapped distance matrices start at dMS[1] and go to dMS[nBootstraps+1]
			// therefore bootstrap is 1-indexed
			int b=1;
			while(b < pars->nAmovaRuns){

				// fill dMS with bootstrapped distance matrices
				fprintf(stderr, "\n\t-> Bootstrapping %d/%d", b, args->nBootstraps);

				prepare_bootstrap_blocks(VCF, pars, args, dMS[b], sampleSt, metadataSt, formulaSt, outSt, blobSt);

				// // perform AMOVA using the bootstrapped dMS
				ASSERT(AMOVA::doAMOVA(dMS[b], metadataSt, sampleSt, outSt->out_amova_fs->fp, pars->LUT_inds2idx) == 0);

				fprintf(stderr, "\n\t-> Finished running AMOVA for bootstrap %d/%d", b, args->nBootstraps);
				b++;
			}
			DATA::blobStruct_destroy(blobSt);
		}

		for(int a=0; a < pars->nAmovaRuns; a++){
			ASSERT(AMOVA::doAMOVA(dMS[a], metadataSt, sampleSt, outSt->out_amova_fs->fp, pars->LUT_inds2idx) == 0);
		}

		for(int a=0; a < pars->nAmovaRuns; a++){
			delete dMS[a];
		}
		delete[] dMS;


		delete metadataSt;
	}
	
	fprintf(stderr, "Total number of sites processed: %lu\n", pars->totSites);
	fprintf(stderr, "Total number of sites skipped for all individual pairs: %lu\n", pars->totSites - pars->nSites);

	outSt->flushAll();

	for (int i = 0; i < pars->nIndCmb; i++)
	{
		delete pairSt[i];
	}
	delete[] pairSt;
	delete formulaSt;
	delete sampleSt;

	VCF::vcfData_destroy(VCF);



}

void input_DM(argStruct *args, paramStruct *pars, DATA::formulaStruct *formulaSt, IO::outFilesStruct *outSt)
{

	fprintf(stderr, "\n[INFO]\t-> Input file type: Distance Matrix\n");

	FILE *in_dm_fp = IO::getFile(args->in_dm_fn, "r");
	DATA::distanceMatrixStruct *dMS = DATA::distanceMatrixStruct_read(in_dm_fp, pars, args);

	FILE *in_mtd_fp = IO::getFile(args->in_mtd_fn, "r");
	DATA::sampleStruct *sampleSt = new DATA::sampleStruct();
	sampleSt->init(pars->nInd);

	DATA::metadataStruct *metadataSt = DATA::metadataStruct_get(in_mtd_fp, sampleSt, formulaSt, args->hasColNames, pars);
	FCLOSE(in_mtd_fp);
	ASSERT(AMOVA::doAMOVA(dMS, metadataSt, sampleSt, outSt->out_amova_fs->fp, pars->LUT_inds2idx) == 0);
	fprintf(stderr, "\n\t-> Finished running AMOVA\n");
}


int main(int argc, char **argv)
{

	if (argc == 1) usage(stderr);

	argStruct *args = argStruct_get(--argc, ++argv);
	paramStruct *pars = paramStruct_init(args);
	IO::outFilesStruct *outSt = new IO::outFilesStruct(args);

	char *DATETIME = pars->DATETIME;
	DATETIME = get_time();
	fprintf(stderr, "\n%s", DATETIME);

	argStruct_print(stderr, args);

	DATA::formulaStruct *formulaSt = NULL;
	if (args->formula != NULL) formulaSt = DATA::formulaStruct_get(args->formula);

	// SWITCH: input file type
	switch (pars->in_ft)
	{

		// --------------------------- INPUT: VCF/BCF --------------------------- //
		case IN_VCF:
		{
			input_VCF(args, pars, formulaSt, outSt);
			break;
		}

		// ---------------------------- INPUT: DISTANCE MATRIX ---------------------------- //
		case IN_DM:
		{
			input_DM(args, pars, formulaSt, outSt);
			break;
		}

		// ---------------------------- INPUT: UNRECOGNIZED ---------------------------- //
		default:
		{
			fprintf(stderr, "\n[ERROR]\t-> Input file type not recognized\n");
			exit(1);
			break;
		}

	}
	
	delete outSt;
	argStruct_destroy(args);
	paramStruct_destroy(pars);

	return 0;
}
