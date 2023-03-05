/*
 * ngsAMOVA
 */

#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

#include <stdio.h>
#include <inttypes.h>
#include <limits>
#include <math.h>
#include <time.h>
#include <pthread.h>

#include "dataStructs.h"
#include "io.h"
#include "argStruct.h"
#include "paramStruct.h"
#include "shared.h"
#include "mathUtils.h"
#include "vcfUtils.h"
#include "em.h"
#include "amova.h"
#include "bootstrap.h"
#include "dev.h"

using size_t = decltype(sizeof(int));

void spawnThreads_pairEM_GL(argStruct *args, paramStruct *pars, pairStruct **pairSt, vcfData *vcfd, IO::outFilesStruct *outSt, distanceMatrixStruct *distMatrix)
{

	pthread_t pairThreads[pars->nIndCmb];
	threadStruct **PTHREADS = new threadStruct*[pars->nIndCmb];

	for (int i = 0; i < pars->nIndCmb; i++)
	{
		PTHREADS[i] = new threadStruct(pairSt[i], vcfd->lngl, args, pars);
	}

	int nJobs_sent = 0;

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
		pairStruct *pair = PTHREADS[pidx]->pair;
		ASSERT(pair->snSites > 0);

		for(int g=0; g<vcfd->nJointClasses; g++)
		{
			vcfd->JointGenoCountDistGL[pidx][g] = pair->optim_jointGenoCountDist[g];
			vcfd->JointGenoProbDistGL[pidx][g] = pair->optim_jointGenoProbDist[g];
		}
		vcfd->JointGenoCountDistGL[pidx][vcfd->nJointClasses] = pair->snSites;
		vcfd->JointGenoProbDistGL[pidx][vcfd->nJointClasses] = pair->snSites;

		if (args->do_square_distance == 1)
		{
			distMatrix->M[pidx] = (double)SQUARE((MATH::EST::Dij(vcfd->JointGenoProbDistGL[pidx])));
		}
		else
		{
			distMatrix->M[pidx] = (double)MATH::EST::Dij(vcfd->JointGenoProbDistGL[pidx]);
		}
		delete PTHREADS[pidx];
	}
	vcfd->print_JointGenoProbDist(outSt, args);
	vcfd->print_JointGenoCountDist(outSt, args);

	delete[] PTHREADS;

}


// prepare distance matrix using original data
void prepare_distanceMatrix(argStruct *args, paramStruct *pars, distanceMatrixStruct *dMS_orig, vcfData *vcfd, pairStruct **pairSt, formulaStruct *formulaSt, IO::outFilesStruct *outSt, blobStruct *blobSt, sampleStruct *sampleSt)
{

	switch (args->doAMOVA)
	{

	// --------------------------- doAMOVA 1 --------------------------- //
	case 1:
	{

		if (args->blockSize == 0)
		{
			readSites_GL(vcfd, args, pars, pairSt);
		}
		else
		{

			blobSt = blobStruct_init(vcfd->nContigs, args->blockSize, vcfd->hdr);
			readSites_GL(vcfd, args, pars, pairSt, blobSt);
		}

		spawnThreads_pairEM_GL(args, pars, pairSt, vcfd, outSt, dMS_orig);

		break;
	}

	// --------------------------- doAMOVA 2 --------------------------- //
	case 2:
	{

		if (args->blockSize == 0)
		{

			readSites_GT(vcfd, args, pars, pairSt);
		}
		else
		{

			fprintf(stderr, "\n[ERROR]\t-> Not implemented yet. Please use -blockSize 0\n");
			exit(1);
			// blobStruct *blobSt = blobStruct_init(vcfd->nContigs, args->blockSize, vcfd->hdr);
			// readSites_GT(vcfd,args, pars, pairSt, blobSt);
			// blobStruct_destroy(blobSt);
		}

		for (int pidx = 0; pidx < pars->nIndCmb; pidx++)
		{
			int snSites = vcfd->JointGenoCountDistGT[pidx][vcfd->nJointClasses];
			if (snSites==0)
			{
				fprintf(stderr,"\n[ERROR]\t-> No shared sites found for pair %d (snSites=%d). This is currently not allowed.\n",pidx,snSites);
				exit(1);
			}

			if (args->do_square_distance == 1)
			{
				dMS_orig->M[pidx] = (double)SQUARE(MATH::EST::Dij(vcfd->JointGenoCountDistGT[pidx], snSites));
			}
			else
			{
				dMS_orig->M[pidx] = (double)MATH::EST::Dij(vcfd->JointGenoCountDistGT[pidx], snSites);
			}
		}

		vcfd->print_JointGenoCountDist(outSt, args);

		break;
	}

	// --------------------------- doAMOVA 3 --------------------------- //
	case 3:
	{
		break;
	}

	// --------------------------- doAMOVA 0 --------------------------- //
	case 0:
	{

		// we do not want to calculate AMOVA, but we want to do EM and obtain the distance matrix
		if (args->doEM != 0)
		{

			ASSERT(args->blockSize == 0); // not implemented yet

			readSites_GL(vcfd, args, pars, pairSt);

			spawnThreads_pairEM_GL(args, pars, pairSt, vcfd, outSt, dMS_orig);
		}

		break;
	}

	// ---------------------- doAMOVA NOT in {0,1,2,3} ---------------------- //
	default:
	{
		fprintf(stderr, "\n[ERROR]\t-> -doAMOVA %d not recognized\n", args->doAMOVA);
		exit(1);
		break;
	}
	}

	if (args->printMatrix != 0)
	{
		dMS_orig->print(outSt->out_dm_fs);
	}
}

void input_VCF_noAMOVA(argStruct *args, paramStruct *pars, formulaStruct *formulaSt, IO::outFilesStruct *outSt, sampleStruct *sampleSt, vcfData *vcfd, pairStruct **pairSt)
{
	// --------------------------- doAMOVA 0 --------------------------- //
	// DO NOT run AMOVA
	fprintf(stderr, "\n[INFO]\t-> -doAMOVA 0; will not run AMOVA\n");
	fprintf(stderr, "\n\t-> Skipped running AMOVA\n");
	pars->nAmovaRuns = 0;

	if (args->doEM != 0)
	{
		distanceMatrixStruct *dMS = new distanceMatrixStruct(pars->nInd, pars->nIndCmb, args->do_square_distance);
		prepare_distanceMatrix(args, pars, dMS, vcfd, pairSt, formulaSt, outSt, NULL, sampleSt);
		if (args->printMatrix != 0)
		{
			dMS->print(outSt->out_dm_fs);
		}
		delete dMS;
	}
}

void input_VCF_doAMOVA(argStruct *args, paramStruct *pars, formulaStruct *formulaSt, IO::outFilesStruct *outSt, sampleStruct *sampleSt, vcfData *vcfd, pairStruct **pairSt)
{
	// --------------------------- doAMOVA!=0 --------------------------- //
	// will run AMOVA
	//
	// prepare for AMOVA

	// [BEGIN] ----------------------- READ METADATA ------------------------- //
	// If vcfd/BCF input, read metadata after vcfd
	// so you can compare VCF records with metadata when reading metadata
	ASSERT(args->in_mtd_fn != NULL);

	FILE *in_mtd_fp = IO::getFile(args->in_mtd_fn, "r");
	metadataStruct *metadataSt = metadataStruct_get(in_mtd_fp, sampleSt, formulaSt, args->hasColNames,pars);
	FCLOSE(in_mtd_fp);
	// [END] ----------------------- READ METADATA ------------------------- //

	pars->nAmovaRuns = 1;
	if (args->nBootstraps > 0)
		pars->nAmovaRuns += args->nBootstraps;

	blobStruct *blobSt = NULL;

	distanceMatrixStruct **dMS = new distanceMatrixStruct *[pars->nAmovaRuns];

	for (int r = 0; r < pars->nAmovaRuns; r++)
	{
		dMS[r] = new distanceMatrixStruct(pars->nInd, pars->nIndCmb, args->do_square_distance);
	}

	if (args->blockSize != 0)
	{
		blobSt = blobStruct_init(vcfd->nContigs, args->blockSize, vcfd->hdr);
	}

	// prepare distance matrix dMS given analysis type (-doAMOVA)

	// dMS[0] is the original distance matrix (not bootstrapped)

	prepare_distanceMatrix(args, pars, dMS[0], vcfd, pairSt, formulaSt, outSt, blobSt, sampleSt);

	if (args->nBootstraps > 0)
	{

		// dMS[0] is already filled with the original distance matrix
		// bootstrapped distance matrices start at dMS[1] and go to dMS[nBootstraps+1]
		// therefore bootstrap is 1-indexed
		int b = 1;
		while (b < pars->nAmovaRuns)
		{
			// fill dMS with bootstrapped distance matrices
			fprintf(stderr, "\n\t-> Bootstrapping %d/%d", b, args->nBootstraps);

			prepare_bootstrap_blocks(vcfd, pars, args, dMS[b], sampleSt, metadataSt, formulaSt, outSt, blobSt);

			// // perform AMOVA using the bootstrapped dMS
			// ASSERT(AMOVA::doAMOVA(dMS[b], metadataSt, sampleSt, outSt->out_amova_fs->fp, pars->lut_indsToIdx) == 0);

			fprintf(stderr, "\n\t-> Finished running AMOVA for bootstrap %d/%d", b, args->nBootstraps);
			b++;
		}
		blobStruct_destroy(blobSt);
	}

	AMOVA::amovaStruct **amv = new AMOVA::amovaStruct *[pars->nAmovaRuns];

	for (int a = 0; a < pars->nAmovaRuns; a++)
	{
		amv[a] = AMOVA::doAmova(dMS[a], metadataSt, sampleSt, pars);
		// ASSERT(AMOVA::doAMOVA(dMS[a], metadataSt, sampleSt, outSt->out_amova_fs->fp, pars->lut_indsToIdx) == 0);
		if (a==0 && args->printAmovaTable == 1)
		{
			amv[a]->print_as_table(stdout, metadataSt);
		}
		if (a==0){
			amv[a]->print_as_csv(outSt->out_amova_fs->fp, metadataSt);
		}
		//TODO print bootstrapped AMOVA tables and distance matrices too
	}



	for (int a = 0; a < pars->nAmovaRuns; a++)
	{
		delete dMS[a];
		delete amv[a];
	}
	delete[] dMS;
	delete[] amv;


	delete metadataSt;
}

// --------------------------- INPUT: VCF/BCF --------------------------- //
void input_VCF(argStruct *args, paramStruct *pars, formulaStruct *formulaSt, IO::outFilesStruct *outSt)
{
	sampleStruct *sampleSt = new sampleStruct();
	vcfData *vcfd = vcfData_init(args, pars, sampleSt);
	pairStruct **pairSt = new pairStruct *[pars->nIndCmb];

	for (int pidx = 0; pidx < pars->nIndCmb; pidx++)
	{
		pairSt[pidx] = new pairStruct(pars, pidx);
	}

	if (args->windowSize == 0)
	{
		// --------------------------- windowSize 0 --------------------------- //
		fprintf(stderr, "\n[INFO]\t-> -windowSize 0; will not use sliding window\n");
		// treat the whole genome as one window
		// args->windowSize = vcfd->nSites;
		// fprintf(stderr, "\n[INFO]\t-> -windowSize %d; will use the whole genome as one window. Therefore -windowSize is set to the total number of sites in the VCF/BCF file (nSites: %d)\n", args->windowSize, vcfd->nSites);
	}
	else
	{
		// --------------------------- windowSize > 0 --------------------------- //
		ASSERT(0==1);
	}

	// if (args->doEM != 0 && args->printMatrix == 1)
	// {
	// }

	if (args->doAMOVA == 0)
	{
		input_VCF_noAMOVA(args, pars, formulaSt, outSt, sampleSt, vcfd, pairSt);
	}
	else
	{
		input_VCF_doAMOVA(args, pars, formulaSt, outSt, sampleSt, vcfd, pairSt);
	}

	fprintf(stderr, "Total number of sites processed: %lu\n", pars->totSites);
	fprintf(stderr, "Total number of sites skipped for all individual pairs: %lu\n", pars->totSites - pars->nSites);

	// outSt->flushAll();

	for (int i = 0; i < pars->nIndCmb; i++)
	{
		delete pairSt[i];
	}
	delete[] pairSt;
	delete formulaSt;
	delete sampleSt;

	vcfData_destroy(vcfd);
}

// --------------------------- INPUT: DISTANCE MATRIX --------------------------- //
void input_DM(argStruct *args, paramStruct *pars, formulaStruct *formulaSt, IO::outFilesStruct *outSt)
{

	if (args->blockSize != 0)
	{
		fprintf(stderr, "\n[ERROR] blockSize is not supported for distance matrix input. Please set blockSize to 0.\n");
		exit(0);
	}

	pars->vprint(1, "input_DM is running\n");

	sampleStruct *sampleSt = new sampleStruct();

	if (args->doAMOVA == 0)
	{
		fprintf(stderr, "\n\t-> Nothing to do.\n");
		exit(1);
	}else{

		// will run AMOVA; prepare for AMOVA

		// [BEGIN] ----------------------- READ METADATA ------------------------- //
		//
		// Metadata reading
		// 	sets pars->nInd
		// 	sets pars->nIndCmb
		// 	sets pars->lut_indsToIdx
		// 	sets pars->lut_idxToInds
		ASSERT(args->in_mtd_fn != NULL);
		FILE *in_mtd_fp = IO::getFile(args->in_mtd_fn, "r");
		metadataStruct *metadataSt = metadataStruct_get(in_mtd_fp, sampleSt, formulaSt, args->hasColNames,pars);
		FCLOSE(in_mtd_fp);
		// [END] ----------------------- READ METADATA ------------------------- //

		distanceMatrixStruct *dMS = distanceMatrixStruct_read_csv(pars, args, metadataSt);

		pairStruct **pairSt = new pairStruct *[pars->nIndCmb];
		for (int pidx = 0; pidx < pars->nIndCmb; pidx++)
		{
			pairSt[pidx] = new pairStruct(pars, pidx);
		}

		AMOVA::amovaStruct *amv = AMOVA::doAmova(dMS, metadataSt, sampleSt, pars);

		if (args->printAmovaTable == 1)
		{
			amv->print_as_table(stdout, metadataSt);
		}
		amv->print_as_csv(outSt->out_amova_fs->fp, metadataSt);

		delete amv;

		// outSt->flushAll();
		delete metadataSt;

		for (int i = 0; i < pars->nIndCmb; i++)
		{
			delete pairSt[i];
		}
		delete[] pairSt;

		delete dMS;
	}

	delete formulaSt;
	delete sampleSt;
}

int main(int argc, char **argv)
{

	if (argc == 1)
		usage(stderr);

	argStruct *args = argStruct_get(--argc, ++argv);
	paramStruct *pars = paramStruct_init(args);
	IO::outFilesStruct *outSt = new IO::outFilesStruct(args);

	char *DATETIME = pars->DATETIME;
	DATETIME = get_time();
	fprintf(stderr, "\n%s", DATETIME);

	argStruct_print(stderr, args);

	formulaStruct *formulaSt = NULL;
	if (args->formula != NULL)
		formulaSt = formulaStruct_get(args->formula);

	// SWITCH: input file type
	switch (pars->in_ft)
	{

	// --------------------------- INPUT: VCF/BCF --------------------------- //
	case IN_VCF:
	{
		if (args->printDev == 1)
		{
			DEV_input_VCF(args, pars, formulaSt, outSt);
			break;
		}

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
