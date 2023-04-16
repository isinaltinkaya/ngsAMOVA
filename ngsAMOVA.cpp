/*
 * ngsAMOVA
 */

#include "dataStructs.h"
#include "io.h"
#include "argStruct.h"
#include "paramStruct.h"
#include "shared.h"
#include "mathUtils.h"
#include "vcfReader.h"
#include "em.h"
#include "amova.h"
#include "bootstrap.h"
#include "evaluation.h"
#include "dxy.h"
#include "neighborJoining.h"
#include "dev.h"

// TODO check size_t
using size_t = decltype(sizeof(int));

// prepare distance matrix using original data
void prepare_distanceMatrix(argStruct *args, paramStruct *pars, distanceMatrixStruct *dMS_orig, vcfData *vcfd, pairStruct **pairSt, formulaStruct *formulaSt, blobStruct *blobSt)
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

			NEVER;
			if (0 != args->nBootstraps)
			{
				blobSt = blobStruct_get(vcfd, args);
			}
			// readSites_GL(vcfd, args, pars, pairSt, blobSt);
		}

		spawnThreads_pairEM(args, pars, pairSt, vcfd, dMS_orig);

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

			fprintf(stderr, "\n[ERROR]\t-> Not implemented yet.\n");
			exit(1);
			if (0 != args->nBootstraps)
			{
				blobSt = blobStruct_get(vcfd, args);
			}
			// readSites_GT(vcfd,args, pars, pairSt, blobSt);
			// DELETE(blobSt);
		}

		for (int pidx = 0; pidx < pars->nIndCmb; pidx++)
		{
			int snSites = vcfd->JointGenoCountDistGT[pidx][vcfd->nJointClasses];
			if (snSites == 0)
			{
				fprintf(stderr, "\n[ERROR]\t-> No shared sites found for pair %d (snSites=%d). This is currently not allowed.\n", pidx, snSites);
				exit(1);
			}

			if (args->squareDistance == 1)
			{
				dMS_orig->M[pidx] = (double)SQUARE(MATH::EST::Dij(vcfd->JointGenoCountDistGT[pidx], snSites));
			}
			else
			{
				dMS_orig->M[pidx] = (double)MATH::EST::Dij(vcfd->JointGenoCountDistGT[pidx], snSites);
			}
		}

		vcfd->print_JointGenoCountDist(args);

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

			spawnThreads_pairEM(args, pars, pairSt, vcfd, dMS_orig);
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

	eval_distanceMatrixStruct(pars, dMS_orig);

	if (args->printMatrix != 0)
	{
		dMS_orig->print(outFiles->out_dm_fs);
	}
}

// --------------------------- INPUT: VCF/BCF --------------------------- //
void input_VCF(argStruct *args, paramStruct *pars, formulaStruct *formulaSt)
{

	vcfData *vcfd = vcfData_init(args, pars);
	pairStruct **pairSt = new pairStruct *[pars->nIndCmb];

	for (int i1 = 0; i1 < vcfd->nInd - 1; i1++)
	{
		for (int i2 = i1 + 1; i2 < vcfd->nInd; i2++)
		{
			int pidx = nCk_idx(vcfd->nInd, i1, i2);
			pairSt[pidx] = new pairStruct(pars, pidx, i1, i2);
		}
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
		ASSERT(0 == 1);
	}

	// --------------------------- doAMOVA 0 --------------------------- //
	// DO NOT run AMOVA
	if (args->doAMOVA == 0 && args->doEM == 1)
	{
		// do not run AMOVA, but do EM and get distance matrix
		distanceMatrixStruct *dMS = new distanceMatrixStruct(pars->nInd, pars->nIndCmb, args->squareDistance, NULL);
		prepare_distanceMatrix(args, pars, dMS, vcfd, pairSt, formulaSt, NULL);
		if (args->printMatrix != 0)
		{
			dMS->print(outFiles->out_dm_fs);
		}
		delete dMS;

		for (int i = 0; i < pars->nIndCmb; i++)
		{
			delete pairSt[i];
		}
		delete[] pairSt;

		vcfData_destroy(vcfd);
		return;
	}

	metadataStruct *metadataSt = metadataStruct_get(args, pars, formulaSt);

	pars->nAmovaRuns = 1;
	if (args->nBootstraps > 0)
	{
		pars->nAmovaRuns += args->nBootstraps;
	}

	blobStruct *blobSt = NULL;

	distanceMatrixStruct **dMS = new distanceMatrixStruct *[pars->nAmovaRuns];

	for (int r = 0; r < pars->nAmovaRuns; r++)
	{
		dMS[r] = new distanceMatrixStruct(pars->nInd, pars->nIndCmb, args->squareDistance, metadataSt->indNames);
	}

	if (0 != args->nBootstraps)
	{
		blobSt = blobStruct_get(vcfd, args);
	}

	// prepare distance matrix dMS given analysis type (-doAMOVA)

	// dMS[0] is the original distance matrix (not bootstrapped)

	prepare_distanceMatrix(args, pars, dMS[0], vcfd, pairSt, formulaSt, blobSt);

	// bootstrapDataset* bootstrap = bootstrapDataset_get(vcfd, pars, args, dMS[b], metadataSt, formulaSt, blobSt);

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


			// // perform AMOVA using the bootstrapped dMS

			fprintf(stderr, "\n\t-> Finished running AMOVA for bootstrap %d/%d", b, args->nBootstraps);
			b++;
		}
		DELETE(blobSt);
	}

	dxyStruct *dxySt = NULL;
	if (args->doDxy > 0)
	{
		if (args->in_dxy_fn == NULL)
		{
			// TODO: implement this as an independent function that does not require to go through amova stuff
			//  pars->in_ft=IN_DXY;
			dxySt = dxyStruct_get(args, pars, dMS[0], metadataSt);
		}
		else
		{
			dxySt = dxyStruct_read(args, pars, dMS[0], metadataSt);
		}
	}

	// TODO print in csv format
	//  for (int i1=0; i1<metadataSt->nInd-1; i1++)
	//  {
	//  	for (int i2=i1+1; i2<metadataSt->nInd; i2++)
	//  	{
	//  		int pidx = nCk_idx(metadataSt->nInd, i1, i2);
	//  		fprintf(stdout,"%d,%d,%d,", i1,i2,pidx);
	//  		fprintf(stdout,"%.*f\n",DBL_MAXDIG10, dMS[0]->M[pidx]);
	//  	}
	//  }

	njStruct *njSt = NULL;
	if (args->doNJ == 1)
	{
		ASSERT(dMS[0] != NULL);
		njSt = njStruct_get(args, pars, dMS[0]);
		njSt->print(outFiles->out_nj_fs);
	}
	else if (args->doNJ == 2)
	{
		ASSERT(dxySt != NULL);
		njSt = njStruct_get(args, pars, dxySt);
		njSt->print(outFiles->out_nj_fs);
	}

	AMOVA::amovaStruct **amv = new AMOVA::amovaStruct *[pars->nAmovaRuns];

	for (int a = 0; a < pars->nAmovaRuns; a++)
	{
		amv[a] = AMOVA::doAmova(dMS[a], metadataSt, pars);
		eval_amovaStruct(amv[a]);
		if (a == 0 && args->printAmovaTable == 1)
		{
			amv[a]->print_as_table(stdout, metadataSt);
		}
		if (a == 0)
		{
			amv[a]->print_as_csv(outFiles->out_amova_fs->fp, metadataSt);
		}
	}

	fprintf(stderr, "Total number of sites processed: %lu\n", pars->totSites);
	fprintf(stderr, "Total number of sites skipped for all individual pairs: %lu\n", pars->totSites - pars->nSites);

	for (int a = 0; a < pars->nAmovaRuns; a++)
	{
		delete dMS[a];
		delete amv[a];
	}
	delete[] dMS;
	delete[] amv;

	delete metadataSt;

	DELETE(njSt);
	DELETE(dxySt);

	for (int i = 0; i < pars->nIndCmb; i++)
	{
		delete pairSt[i];
	}
	delete[] pairSt;

	vcfData_destroy(vcfd);
}

// --------------------------- INPUT: DISTANCE MATRIX --------------------------- //
void input_DM(argStruct *args, paramStruct *pars, formulaStruct *formulaSt)
{

	if (args->blockSize != 0)
	{
		fprintf(stderr, "\n[ERROR] blockSize is not supported for distance matrix input. Please set blockSize to 0.\n");
		exit(0);
	}

	IO::vprint(1, "input_DM is running\n");

	distanceMatrixStruct *dMS = distanceMatrixStruct_read(pars, args);

	metadataStruct *metadataSt = metadataStruct_get(args, pars, formulaSt);

	dMS->set_item_labels(metadataSt->indNames);

	pairStruct **pairSt = new pairStruct *[pars->nIndCmb];

	for (int i1 = 0; i1 < pars->nInd - 1; i1++)
	{
		for (int i2 = i1 + 1; i2 < pars->nInd; i2++)
		{
			int pidx = nCk_idx(pars->nInd, i1, i2);
			pairSt[pidx] = new pairStruct(pars, pidx, i1, i2);
		}
	}

	if (args->doAMOVA != 0)
	{
		AMOVA::amovaStruct *amv = AMOVA::doAmova(dMS, metadataSt, pars);
		eval_amovaStruct(amv);
		if (args->printAmovaTable == 1)
		{
			amv->print_as_table(stdout, metadataSt);
		}
		amv->print_as_csv(outFiles->out_amova_fs->fp, metadataSt);
		delete amv;
	}

	dxyStruct *dxySt = NULL;
	if (args->doDxy > 0)
	{
		if (args->in_dxy_fn == NULL)
		{
			dxySt = dxyStruct_get(args, pars, dMS, metadataSt);
		}
		else
		{
			dxySt = dxyStruct_read(args, pars, dMS, metadataSt);
		}
	}

	njStruct *njSt = NULL;
	if (args->doNJ == 1)
	{
		ASSERT(dMS != NULL);
		njSt = njStruct_get(args, pars, dMS);
		njSt->print(outFiles->out_nj_fs);
	}
	else if (args->doNJ == 2)
	{
		ASSERT(dxySt != NULL);
		njSt = njStruct_get(args, pars, dxySt);
		njSt->print(outFiles->out_nj_fs);
	}

	delete metadataSt;

	for (int i = 0; i < pars->nIndCmb; i++)
	{
		delete pairSt[i];
	}
	delete[] pairSt;

	DELETE(dxySt);
	DELETE(njSt);

	delete dMS;
}

int main(int argc, char **argv)
{
	if (argc == 1)
	{
		print_help(stderr);
		exit(0);
	}

	argStruct *args = argStruct_get(--argc, ++argv);
	paramStruct *pars = paramStruct_init(args);
	IO::outFilesStruct_set(args, outFiles);

	char *DATETIME = pars->DATETIME;
	DATETIME = get_time();
	fprintf(stderr, "\n%s", DATETIME);

	argStruct_print(stderr, args);

	formulaStruct *formulaSt = formulaStruct_get(args->formula);
	formulaSt->print(stderr);

	// TODO put formulaStruct into pars
	// TODO use bitsetting to represent file types VCF/BCF, DM, etc.

	// == OBTAIN THE DISTANCE MATRIX ==

	// determine the method to obtain the distance matrix
	switch (pars->in_ft)
	{

	// --------------------------- INPUT: VCF/BCF --------------------------- //
	case IN_VCF:
		input_VCF(args, pars, formulaSt);
		break;

	// ---------------------------- INPUT: DISTANCE MATRIX ---------------------------- //
	case IN_DM:
		input_DM(args, pars, formulaSt);
		break;

	// ---------------------------- INPUT: JGCD ---------------------------- //
	case IN_JGCD:
		// input_JGCD(args, pars, formulaSt);
		break;

	// ---------------------------- INPUT: UNRECOGNIZED ---------------------------- //
	default:
		fprintf(stderr, "\n[ERROR]\t-> Input file type not recognized\n");
		exit(1);
	}

	// == CLEANUP ==
	formulaStruct_destroy(formulaSt);
	argStruct_destroy(args);
	paramStruct_destroy(pars);
	IO::outFilesStruct_destroy(outFiles);

	return 0;
}
