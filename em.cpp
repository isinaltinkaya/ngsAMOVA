#include "em.h"

/// @brief spawnThreads_pairEM spawn threads for running EM algorithm for each individual pair 
/// @param args 
/// @param pars 
/// @param pairSt 
/// @param vcfd 
/// @param outSt 
/// @param distMatrix 
void spawnThreads_pairEM(argStruct *args, paramStruct *pars, pairStruct **pairSt, vcfData *vcfd,  distanceMatrixStruct *distMatrix)
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
			// vcfd->JointGenoCountDistGL[pidx][g] = pair->optim_jointGenoCountDist[g];
			vcfd->JointGenoCountDistGL[pidx][g] = pair->optim_jointGenoProbDist[g]*pair->snSites;
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
	vcfd->print_JointGenoProbDist( args);
	vcfd->print_JointGenoCountDist( args);

	delete[] PTHREADS;

}


/// @brief thread handler for EM_2DSFS_GL3
/// @param p 
/// @return 
void* t_EM_2DSFS_GL3(void* p){

	threadStruct* THREAD= (threadStruct*) p;

	if (EM_2DSFS_GL3(THREAD) != 0){
		fprintf(stderr,"\n[ERROR] Problem with EM\n");
		exit(1);
	}
	return 0;
}


/// @brief EM algorithm for 3x3 PGC (pairwise genotype categories)
/// 	2DSFS from 3 GL values
/// @param THREAD 
/// @return 
int EM_2DSFS_GL3(threadStruct* THREAD)
{

	double **lngls=THREAD->lngls;
	pairStruct* pair=THREAD->pair;
	const double tole = THREAD->tole;
	const int mEmIter=THREAD->maxEmIter;

	// const int i1=lut_idxToInds[pair->idx][0];
	// const int i2=lut_idxToInds[pair->idx][1];
	const int i1=pair->i1;
	const int i2=pair->i2;

	double sum=0.0;
	double d=0.0;

	double tmp_jointGenoProb=0.0;



	// set initial guess: 1/9 flat prior
	for (int i=0; i<9; i++){
		pair->optim_jointGenoProbDist[i]=FRAC_1_9;
	}

	do{
		
		if(pair->n_em_iter >= mEmIter){
			break;
		}

		// double TMP[9];
		// for (int i=0; i<9; i++){
		// 	TMP[i]=0.0;
		// }

		// double ESFS[9];
		// for (int i=0; i<9; i++){
		// 	ESFS[i]=0.0;
		// }

		// //loop through shared sites for pair
		// for(size_t sn=0; sn<pair->snSites; sn++){
		// 	size_t s=pair->sharedSites[sn];
		// 	sum=0.0;

		// 	// SFS * ind1 * ind2
		// 	//lngls3 (anc,anc),(anc,der),(der,der)
		// 	for(int i=0;i<9; i++){
		// 		TMP[i] = pair->optim_jointGenoProbDist[i] * exp( lngls[s][(3*i1)+
		// 		TMP[i] = pair->optim_jointGenoProbDist[i] * exp( lngls[s][(3*i1)+(i/3)] + lngls[s][(3*i2)+(i%3)]);
		// 		sum += TMP[i];
		// 	}


		// 	for(int i=0; i<9; i++){
		// 		ESFS[i] +=  TMP[i]/sum;
		// 	}

		// }

		// d=0.0;
		// for(int i=0;i<9;i++){
		// 	tmp_jointGenoProb = ESFS[i]/(double)pair->snSites;
		// 	d += fabs(tmp_jointGenoProb - pair->optim_jointGenoProbDist[i]);
		// 	pair->optim_jointGenoProbDist[i]=tmp_jointGenoProb;
		// 	pair->optim_jointGenoCountDist[i]=ESFS[i];
		// }
		// @@
		double TMP[3][3];
		double ESFS[3][3];
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
					TMP[i][j] = pair->optim_jointGenoProbDist[i*3+j] * exp( lngls[s][(3*i1)+i] + lngls[s][(3*i2)+j]);
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
				tmp_jointGenoProb=ESFS[i][j]/(double)pair->snSites;
				d += fabs(tmp_jointGenoProb - pair->optim_jointGenoProbDist[i*3+j]);
				pair->optim_jointGenoProbDist[i*3+j]=tmp_jointGenoProb;
			}
		}

		pair->n_em_iter++;

#if 0
		fprintf(stdout,"%d,%d,",pair->idx,pair->n_em_iter);
		fprintf(stdout,"%f,", log10(d));
		IO::print::Array(stdout,pair->optim_jointGenoProbDist, 3, 3, ',');
#endif

	}while(d>tole);

	pair->d=d;
	
	return 0;
}




