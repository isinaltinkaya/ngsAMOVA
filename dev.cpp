#include "dev.h"

/// @brief thread handler for EM_2DSFS_GL3 DEV version
/// @param p 
/// @return 
void* DEV_t_EM_2DSFS_GL3(void* p){

	threadStruct* THREAD= (threadStruct*) p;

	if (DEV_EM_2DSFS_GL3(THREAD) != 0){
		fprintf(stderr,"\n[ERROR] Problem with EM\n");
		exit(1);
	}
	return 0;
}


/// @brief EM algorithm for 3x3 PGC (pairwise genotype categories)
/// 		Version: Development (Print per iteration values)
/// @param THREAD 
/// @return 
///
/// outputs per-iteration values to stdout
/// columns in the output:
/// 	pair_index,n_em_iter,sq_dij,log10_d
//
int DEV_EM_2DSFS_GL3(threadStruct* THREAD){

	double **lngls=THREAD->lngls;
	DATA::pairStruct* pair=THREAD->pair;
	const double tole = THREAD->args->tole;
	const int mEmIter=THREAD->args->mEmIter;


	const int i1=pair->pars->LUT_idx2inds[pair->idx][0];
	const int i2=pair->pars->LUT_idx2inds[pair->idx][1];

	double temp;
	double sum;
	double d;


	// initial guess: 1/9 flat prior
	for (int i=0; i<9; i++){
		pair->SFS[i]=(double) 1 / (double) 9;
	}

	do{
		
		if(pair->n_em_iter >= mEmIter){
			break;
		}

		double TMP[3][3];
		// double *TMP;
		double ESFS[3][3];
		// double *ESFS;

		// for (int i=0; i<9; i++){
			// ESFS[i]=0.0;
		// }

			for(int i=0;i<3;i++){
				for(int j=0;j<3;j++){
					ESFS[i][j]=0.0;
				}
			}

		//loop through shared sites for pair
		for(size_t sn=0; sn<pair->snSites; sn++){
			size_t s=pair->sharedSites[sn];
			sum=0.0;

			// SFS * ind1 * ind2
			//lngls3 (anc,anc),(anc,der),(der,der)
			for(int i=0;i<3;i++){
				for(int j=0;j<3;j++){
					TMP[i][j] = pair->SFS[i*3+j] * exp( lngls[s][(3*i1)+i] + lngls[s][(3*i2)+j]);
					sum += TMP[i][j];
				}
			}

			for(int i=0;i<3;i++){
				for(int j=0;j<3;j++){
					ESFS[i][j] += TMP[i][j]/sum;
				}
			}
		}

		d=0.0;
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				temp=ESFS[i][j]/(double)pair->snSites;
				d += fabs(temp - pair->SFS[i*3+j]);
				pair->SFS[i*3+j]=temp;

			}
		}

		pair->n_em_iter++;

		// IO::print::Array(stdout,pair->SFS, 3, 3, ',');

// Column 1 contains the individual pair's index in individual pairs array.
// Column 2 contains the number of em iterations.
// Column 3 contains the squared distance between the individuals in the pair.
// Column 4 contains the difference in likelihood from the previous iteration.
#if 1
		fprintf(stdout,"%d,%d,",pair->idx,pair->n_em_iter);
		fprintf(stdout,"%.*f,",(int)DBL_MAXDIG10, (double) SQUARE(MATH::EST::Dij(pair->SFS)));
		fprintf(stdout,"%.*f",(int)DBL_MAXDIG10, log10(d));
		fprintf(stdout,"\n");
#endif

	}while(d>tole);

	pair->d=d;
	
	return 0;
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