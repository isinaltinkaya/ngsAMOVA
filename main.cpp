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

void spawnThreads_pairEM_GL(argStruct *args, paramStruct *pars, DATA::pairStruct **pairSt, VCF::vcfData *VCF, IO::outFilesStruct *outSt, DATA::distanceMatrixStruct *dMS)
{

	pthread_t pairThreads[pars->n_ind_cmb];
	threadStruct *PTHREADS[pars->n_ind_cmb];

	for (int i = 0; i < pars->n_ind_cmb; i++)
	{
		PTHREADS[i] = new threadStruct(pairSt[i], VCF->lngl, pars->nSites, outSt->out_sfs_fs, args);
	}

	int nJobs_sent = 0;

	fprintf(outSt->out_sfs_fs->fp, "Method,Ind1,Ind2,A,D,G,B,E,H,C,F,I,n_em_iter,shared_nSites,Delta,Tole,Dij,Dij2\n");
	for (int pidx = 0; pidx < pars->n_ind_cmb; pidx++)
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
		t = pars->n_ind_cmb - nJobs_sent;
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

	for (int pidx = 0; pidx < pars->n_ind_cmb; pidx++)
	{

		DATA::pairStruct *pair = PTHREADS[pidx]->pair;
		ASSERT(pair->snSites > 0);

		if (args->do_square_distance == 1)
		{
			dMS->M[pidx] = (double)SQUARE((1 - MATH::EST::Sij(pair->SFS)));
		}
		else
		{
			dMS->M[pidx] = (double)(1 - MATH::EST::Sij(pair->SFS));
		}

		IO::print::Sfs("gle", outSt->out_sfs_fs, pair, args, VCF->hdr->samples[pair->i1], VCF->hdr->samples[pair->i2]);

		delete PTHREADS[pidx];
	}
}

void prepare_bootstrap_block(VCF::vcfData *VCF, paramStruct *pars, argStruct *args, DATA::distanceMatrixStruct *dMS, DATA::sampleStruct *sampleSt, DATA::metadataStruct *mS, DATA::formulaStruct *formulaSt, IO::outFilesStruct *outSt, DATA::contigStruct *contigSt)
{

	fprintf(stderr, "\n\n\n\n################ Bootstrap blocks");
	// Bootstrap genomic blocks among individuals 
	// while keeping the order of blocks the same (based on AMOVA levels)


	if (mS->nLevels == 1)
	{
		// only one level; shuffle all individual blocks
		for(int i=0; i<mS->nIndMetadata; i++)
		{

			for(int contig=0; contig < contigSt->nContigs; contig++){
				//loop through blocks and shuffle
				for(int block=0; block<contigSt->contigNBlocks[contig]; block++)
				{
					int blockStart = contigSt->contigBlockStarts[contig][block];
					// shuffle genomic blocks among all individuals in given level
					// save to a new simulated VCF data

					// choose a random individual among the individual set in the shuffling level
					// then set the block to the chosen individual's block at the same position
					int chosen_block = rand() % mS->nIndMetadata;
					fprintf(stderr, "\n\n\n\n######## chosen_block: %d", chosen_block);
					// bootVCF->lngl[blockStart][chosen_block] = VCF->lngl[blockStart][i];
				}

			}
		}

	}else{
		
		// more than one level; shuffle blocks within each level
		for(int level=0; level<mS->nLevels; level++)
		{
			for(int i=0; i<mS->nIndMetadata; i++)
			{
			}
		}
	}
}

int main(int argc, char **argv)
{

	if (argc == 1)
	{
		usage(stderr);
		exit(0);
	}

	argStruct *args = argStruct_get(--argc, ++argv);

	if (args != NULL)
	{

		paramStruct *pars = paramStruct_init(args);
		IO::outFilesStruct *outSt = new IO::outFilesStruct(args);

		char *DATETIME = pars->DATETIME;
		DATETIME = get_time();
		fprintf(stderr, "\n%s", DATETIME);

		argStruct_print(stderr, args);

		DATA::formulaStruct *formulaSt = NULL;
		if (args->formula != NULL)
		{
			fprintf(stderr, "\nAMOVA formula is defined as: %s", args->formula);
			formulaSt = DATA::formulaStruct_get(args->formula);
		}

		// input file type switch
		switch (pars->in_ft)
		{

		// --------------------------- INPUT: VCF/BCF --------------------------- //
		case IN_VCF:
		{
			fprintf(stderr, "\n[INFO]\t-> Input file type: VCF/BCF\n");


			DATA::sampleStruct *sampleSt = new DATA::sampleStruct();

			VCF::vcfData* VCF = VCF::vcfData_init(args, pars, sampleSt);

			ASSERT(pars->n_ind_cmb>0);
			DATA::pairStruct **pairSt = new DATA::pairStruct *[pars->n_ind_cmb];

			for (int i1 = 0; i1 < pars->nInd - 1; i1++)
			{
				for (int i2 = i1 + 1; i2 < pars->nInd; i2++)
				{
					int pidx = pars->LUT_indPairIdx[i1][i2];

					pairSt[pidx] = new DATA::pairStruct();
					pairSt[pidx]->i1 = i1;
					pairSt[pidx]->i2 = i2;
					pairSt[pidx]->idx = pidx;
				}
			}

			int nAmovaRuns=0;

			if(args->doAMOVA == 0){
				// --------------------------- doAMOVA 0 --------------------------- //
				// DO NOT run AMOVA
				fprintf(stderr, "\n[INFO]\t-> -doAMOVA 0; will not run AMOVA\n");
				fprintf(stderr, "\n\t-> Skipped running AMOVA\n");
				nAmovaRuns = 0;

			}else{

				// --------------------------- doAMOVA!=0 --------------------------- //
				// Will run AMOVA
				DATA::metadataStruct *metadataSt = NULL;

				nAmovaRuns = 1;
				if(args->nBootstraps > 0)
				{
					// bootstrap_0 is the original data
					// the rest are the real bootstraps
					nAmovaRuns += args->nBootstraps;
				}

				DATA::distanceMatrixStruct **dMS = new DATA::distanceMatrixStruct *[nAmovaRuns];
				DATA::contigStruct *contigSt = NULL;

				// prepare distance matrix dMS given analysis type (-doAMOVA)
				switch (args->doAMOVA)
				{

					case 0:{

						// should never reach here; control is done above
						// but this ain't Rust, you can never be too safe with C/C++ :P
						ASSERT(0==1); 
						break;
					}



					// --------------------------- doAMOVA 1 --------------------------- //
					case 1:{

						fprintf(stderr, "\n[INFO]\t-> -doAMOVA 1; will use 10 genotype likelihoods from VCF GL tag.\n");


						if(args->blockSize == 0){

							VCF->readSites_GL(args, pars, pairSt);

						}else{

							contigSt = DATA::contigStruct_init(VCF->nContigs, args->blockSize, VCF->hdr);
							VCF->readSites_GL(args, pars, pairSt, contigSt);

						}

						dMS[0] = new DATA::distanceMatrixStruct(pars->nInd, pars->n_ind_cmb, args->do_square_distance);


						// [BEGIN] ----------------------- READ METADATA ------------------------- //
						// If VCF/BCF input, read metadata after VCF
						// so you can compare VCF records with metadata when reading metadata
						if (args->in_mtd_fn!=NULL)
						{
							FILE *in_mtd_fp = IO::getFile(args->in_mtd_fn, "r");
							metadataSt = DATA::metadataStruct_get(in_mtd_fp, sampleSt, formulaSt, args->hasColNames, pars);
							FCLOSE(in_mtd_fp);
						}
						// [END] ----------------------- READ METADATA ------------------------- //

						spawnThreads_pairEM_GL(args, pars, pairSt, VCF, outSt, dMS[0]);


						break;

					}

					// --------------------------- doAMOVA 2 --------------------------- //
					case 2:{
						fprintf(stderr, "\n[INFO]\t-> -doAMOVA 2; will use genotypes from VCF GT tag.\n");

						if(args->blockSize == 0){

							VCF->readSites_GT(args, pars, pairSt);

						}else{

							fprintf(stderr,"\n[ERROR]\t-> Not implemented yet. Please use -blockSize 0\n");
							exit(1);
							// DATA::contigStruct *contigSt = DATA::contigStruct_init(VCF->nContigs, args->blockSize, VCF->hdr);
							// VCF->readSites_GT(args, pars, pairSt, contigSt);
							// DATA::contigStruct_destroy(contigSt);
						}

						dMS[0] = new DATA::distanceMatrixStruct(pars->nInd, pars->n_ind_cmb, args->do_square_distance);

						// [BEGIN] ----------------------- READ METADATA ------------------------- //
						// If VCF/BCF input, read metadata after VCF
						// so you can compare VCF records with metadata when reading metadata
						if (args->in_mtd_fn!=NULL)
						{
							FILE *in_mtd_fp = IO::getFile(args->in_mtd_fn, "r");
							metadataSt = DATA::metadataStruct_get(in_mtd_fp, sampleSt, formulaSt, args->hasColNames, pars);
							FCLOSE(in_mtd_fp);
						}
						// [END] ----------------------- READ METADATA ------------------------- //



						for (int i1 = 0; i1 < pars->nInd - 1; i1++)
						{
							for (int i2 = i1 + 1; i2 < pars->nInd; i2++)
							{

								int pidx = pars->LUT_indPairIdx[i1][i2];

								int snSites = VCF->SFS_GT3[pidx][9];
								if (args->do_square_distance == 1)
								{
									dMS[0]->M[pidx] = (double)SQUARE((1 - MATH::EST::Sij(VCF->SFS_GT3[pidx], snSites)));
								}
								else
								{
									dMS[0]->M[pidx] = (double)(1 - MATH::EST::Sij(VCF->SFS_GT3[pidx], snSites));
								}

								
								IO::print::Sfs("gt", outSt->out_sfs_fs, args, VCF->SFS_GT3[pidx], snSites, VCF->hdr->samples[i1], VCF->hdr->samples[i2]);
							}
						}
						break;

					}

					// --------------------------- doAMOVA 3 --------------------------- //
					case 3:{
						fprintf(stderr, "\n[INFO]\t-> -doAMOVA 3; will use both genotype likelihoods and genotypes from VCF.\n");
						break;
					}

					// ---------------------- doAMOVA NOT in {0,1,2,3} ---------------------- //
					default:{
						fprintf(stderr, "\n[ERROR]\t-> -doAMOVA %d not recognized\n", args->doAMOVA);
						exit(1);
						break;
					}
				}


				// dMS[0] now contains the original distance matrix
				if (args->printMatrix == 1)
				{
					dMS[0]->print(outSt->out_dm_fs->fp);
				}


				if(args->doAMOVA>0){

					if(args->nBootstraps > 0){



						// dMS[0] is already filled with the original distance matrix
						// bootstrapped distance matrices start at dMS[1] and go to dMS[nBootstraps+1]
						// therefore bootstrap is 1-indexed
						int b=0;
						while(b < args->nBootstraps){

							// fill dMS with bootstrapped distance matrices
							// and perform AMOVA on each

							fprintf(stderr, "\n\t-> Bootstrapping %d/%d", b+1, args->nBootstraps);

							dMS[b+1] = new DATA::distanceMatrixStruct(pars->nInd, pars->n_ind_cmb, args->do_square_distance);

							prepare_bootstrap_block(VCF, pars, args, dMS[b+1], sampleSt, metadataSt, formulaSt, outSt, contigSt);

							// // perform AMOVA using the bootstrapped dMS
							ASSERT(AMOVA::doAMOVA(dMS[b+1], metadataSt, sampleSt, outSt->out_amova_fs->fp, pars->LUT_indPairIdx) == 0);
							fprintf(stderr, "\n\t-> Finished running AMOVA for bootstrap %d/%d", b+1, args->nBootstraps);

							b++;
						}


						DATA::contigStruct_destroy(contigSt);

						for(int b=0; b < nAmovaRuns; b++){
							delete dMS[b];
						}
						delete[] dMS;

					}else if(args->nBootstraps == 0){
						// if no bootstrapping, then only use original data for AMOVA

						// perform AMOVA using the original dMS
						ASSERT(AMOVA::doAMOVA(dMS[args->nBootstraps], metadataSt, sampleSt, outSt->out_amova_fs->fp, pars->LUT_indPairIdx) == 0);
						fprintf(stderr, "\n\t-> Finished running AMOVA\n");

						delete dMS[0];
						delete[] dMS;

					}else{
						fprintf(stderr, "\n[ERROR]\t-> -nBootstraps %d not recognized\n", args->nBootstraps);
						exit(1);
					}

				}else{
					//no AMOVA
				}

				delete metadataSt;
			}
			

			fprintf(stderr, "Total number of sites processed: %lu\n", pars->totSites);
			fprintf(stderr, "Total number of sites skipped for all individual pairs: %lu\n", pars->totSites - pars->nSites);

			outSt->flushAll();

			for (int i = 0; i < pars->n_ind_cmb; i++)
			{
				delete pairSt[i];
			}
			delete[] pairSt;

			delete sampleSt;
			delete formulaSt;

			VCF::vcfData_destroy(VCF);
			delete outSt;

			break;
		}

		// ---------------------------- INPUT: DISTANCE MATRIX ---------------------------- //
		// case: input file is distance matrix
		case IN_DM:
		{
			fprintf(stderr, "\n[INFO]\t-> Input file type: Distance Matrix\n");

			FILE *in_dm_fp = IO::getFile(args->in_dm_fn, "r");
			DATA::distanceMatrixStruct *dMS = DATA::distanceMatrixStruct_read(in_dm_fp, pars, args);

			FILE *in_mtd_fp = IO::getFile(args->in_mtd_fn, "r");
			DATA::sampleStruct *sampleSt = new DATA::sampleStruct();
			sampleSt->init(pars->nInd);

			DATA::metadataStruct *metadataSt = DATA::metadataStruct_get(in_mtd_fp, sampleSt, formulaSt, args->hasColNames, pars);
			FCLOSE(in_mtd_fp);
			ASSERT(AMOVA::doAMOVA(dMS, metadataSt, sampleSt, outSt->out_amova_fs->fp, pars->LUT_indPairIdx) == 0);
			fprintf(stderr, "\n\t-> Finished running AMOVA\n");
			break;
		}

		// default: unrecognized input file type
		default:
		{
			fprintf(stderr, "\n[ERROR]\t-> Input file type not recognized\n");
			exit(1);
			break;
		}

			return 0;
		}
		argStruct_destroy(args);
		paramStruct_destroy(pars);
	}
	else
	{
		ASSERT(1 == 0);
	}

	return 0;
}
