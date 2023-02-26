#include "em.h"

//TODO precalculate the flat prior 1/9 so we don't have to do it every time; const extern etc

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
int EM_2DSFS_GL3(threadStruct* THREAD){

	double **lngls=THREAD->lngls;
	pairStruct* pair=THREAD->pair;
	const double tole = THREAD->args->tole;
	const int mEmIter=THREAD->args->mEmIter;

	const int i1=pair->pars->lut_idxToInds[pair->idx][0];
	const int i2=pair->pars->lut_idxToInds[pair->idx][1];

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

		double TMP[9];
		for (int i=0; i<9; i++){
			TMP[i]=0.0;
		}

		double ESFS[9];
		for (int i=0; i<9; i++){
			ESFS[i]=0.0;
		}

		//loop through shared sites for pair
		for(size_t sn=0; sn<pair->snSites; sn++){
			size_t s=pair->sharedSites[sn];
			sum=0.0;

			// SFS * ind1 * ind2
			//lngls3 (anc,anc),(anc,der),(der,der)
			for(int i=0;i<9; i++){
				TMP[i] = pair->optim_jointGenoProbDist[i] * exp( lngls[s][(3*i1)+(i/3)] + lngls[s][(3*i2)+(i%3)]);
				sum += TMP[i];
			}


			for(int i=0; i<9; i++){
				ESFS[i] +=  TMP[i]/sum;
			}

		}

		d=0.0;
		for(int i=0;i<9;i++){
			tmp_jointGenoProb = ESFS[i]/(double)pair->snSites;
			d += fabs(tmp_jointGenoProb - pair->optim_jointGenoProbDist[i]);
			pair->optim_jointGenoProbDist[i]=tmp_jointGenoProb;
			pair->optim_jointGenoCountDist[i]=ESFS[i];
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





// /// @brief thread handler for block_EM_2DSFS_GL3
// /// @param p 
// /// @return 
// void* t_block_EM_2DSFS_GL3(void* p){

// 	threadStruct* THREAD= (threadStruct*) p;

// 	if (EM_2DSFS_GL3(THREAD) != 0){
// 		fprintf(stderr,"\n[ERROR] Problem with EM\n");
// 		exit(1);
// 	}
// 	return 0;
// }


// /// @brief Per-block EM algorithm for 2DSFS from 3 GL values
// /// @param THREAD 
// /// @return 
// int block_EM_2DSFS_GL3(threadStruct* THREAD){


// 	double **lngls=THREAD->lngls;
// 	pairStruct* pair=THREAD->pair;
// 	const double tole = THREAD->args->tole;
// 	const int mEmIter=THREAD->args->mEmIter;


// 	const int i1=pair->pars->lut_idxToInds[pair->idx][0];
// 	const int i2=pair->pars->lut_idxToInds[pair->idx][1];

// 	double temp;
// 	double sum;
// 	double d;


// 	// initial guess: 1/9 flat prior
// 	for (int i=0; i<9; i++){
// 		pair->optim_jointGenoProbDist[i]=FRAC_1_9;
// 	}

// 	do{
		
// 		if(pair->n_em_iter >= mEmIter){
// 			break;
// 		}


// 		double TMP[3][3];
// 		// double *TMP;
// 		double ESFS[3][3];
// 		// double *ESFS;

// 		// for (int i=0; i<9; i++){
// 			// ESFS[i]=0.0;
// 		// }

// 		for(int i=0;i<3;i++){
// 			for(int j=0;j<3;j++){
// 				ESFS[i][j]=0.0;
// 			}
// 		}

// 		//loop through the block
// 		for(size_t sn=0; sn<pair->snSites; sn++){
// 			size_t s=pair->sharedSites[sn];
// 			sum=0.0;

// 			// SFS * ind1 * ind2
// 			//lngls3 (anc,anc),(anc,der),(der,der)
// 			for(int i=0;i<3;i++){
// 				for(int j=0;j<3;j++){
// 					TMP[i][j] = pair->optim_jointGenoProbDist[i*3+j] * exp( lngls[s][(3*i1)+i] + lngls[s][(3*i2)+j]);
// 					sum += TMP[i][j];
// 				}
// 			}

// 			for(int i=0;i<3;i++){
// 				for(int j=0;j<3;j++){
// 					ESFS[i][j] += TMP[i][j]/sum;
// 				}
// 			}
// 		}

// 		d=0.0;
// 		for(int i=0;i<3;i++){
// 			for(int j=0;j<3;j++){
// 				pair->optim_jointGenoCountDist[i*3+j]=ESFS[i][j];
// 				temp=ESFS[i][j]/(double)pair->snSites;
// 				d += fabs(temp - pair->optim_jointGenoProbDist[i*3+j]);
// 				pair->optim_jointGenoProbDist[i*3+j]=temp;

// 			}
// 		}

// 		pair->n_em_iter++;

// #if 0
// 		fprintf(stdout,"%d,%d,",pair->idx,pair->n_em_iter);
// 		fprintf(stdout,"%f,", log10(d));
// 		IO::print::Array(stdout,pair->optim_jointGenoProbDist, 3, 3, ',');
// #endif

// 	}while(d>tole);

// 	pair->d=d;
	
// 	return 0;
// }


// // sum [ log [ sum [ M(g1,g2) * GL(g1) * GL(g2) ] ]
// void test_em(double **lngls3, double GLSFS[3][3],int **GTSFS,int i1,int i2,char *id1, char*id2, int pair_idx,size_t nSites, int snSites, FILE *outfile){
// // void test_em(double **lngls3, double *GLSFS,int **GTSFS,int i1,int i2,char *id1, char*id2, int pair_idx,size_t nSites, int snSites, FILE *outfile){

// 	double glres=0.0;
// 	double gtres=0.0;
// 	for(size_t s=0;s<nSites;s++){
// 		if ((lngls3[s][(3*i1)+0]==NEG_INF) && (lngls3[s][(3*i1)+1]==NEG_INF) && (lngls3[s][(3*i1)+2]==NEG_INF)){
// 			continue;
// 		}else if ((lngls3[s][(3*i2)+0]==NEG_INF) && (lngls3[s][(3*i2)+1]==NEG_INF) && (lngls3[s][(3*i2)+2]==NEG_INF)){
// 			continue;
// 		}else{
// 			double glresi=0.0;
// 			double gtresi=0.0;
// 			for(int i=0;i<3;i++){
// 				for(int j=0;j<3;j++){
// 					glresi += GLSFS[i][j] * exp(lngls3[s][(3*i1)+i]) * exp(lngls3[s][(3*i2)+j]);
// 					gtresi += ((double) GTSFS[pair_idx][(3*i)+j]/(double) snSites)  * exp(lngls3[s][(3*i1)+i]) * exp(lngls3[s][(3*i2)+j]);
// 					// fprintf(stderr,"->->->%d %d\n",GTSFS[pair_idx][(3*i)+j],snSites);
// 				}
// 			}
// 			// fprintf(stderr,"->->%f,%f\n",glresi,gtresi);
// 			glres+=log(glresi);
// 			gtres+=log(gtresi);
// 		}
// 	}

// 	fprintf(outfile,"%s,%s,%f,%f\n",id1,id2,glres,gtres);
// }

