#include "em.h"


void* t_EM_2DSFS_GL3(void* p){

	threadStruct* THREAD= (threadStruct*) p;

	if (EM_2DSFS_GL3(THREAD) != 0){
		fprintf(stderr,"\n[ERROR] Problem with EM\n");
		exit(1);
	}
	return NULL;

}


int EM_2DSFS_GL3(threadStruct* THREAD){


	double **lngls=THREAD->lngls;
	DATA::pairStruct* pair=THREAD->pair;
	FILE* out_sfs_ff=THREAD->out_sfs_ff;
	double tole=THREAD->tole;
	int mEmIter=THREAD->mEmIter;


	// fprintf(stderr,"\nEM begin for ind1:%d and ind2:%d \n",pair->i1,pair->i2);

	double temp;
	double sum;
	double d;


	// initial guess: 1/9 flat prior
	for (int i=0; i<9; i++){
		pair->SFS[i]=(double) 1/ (double) 9;
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

			size_t s=pair->sSites[sn];
			sum=0.0;

			// SFS * ind1 * ind2
			//lngls3 (anc,anc),(anc,der),(der,der)
			for(int i=0;i<3;i++){
				for(int j=0;j<3;j++){
					TMP[i][j] = pair->SFS[i*3+j] * exp( lngls[s][(3*pair->i1)+i] + lngls[s][(3*pair->i2)+j]);
					sum += TMP[i][j];
				}
			}
			// for(int g=0; g<9; g++){
				// TMP[g]
					//
				// pair->SFS[g]
//
				// pair->i1
			// }
//
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
#if 0
		fprintf(stdout,"%d,%d,",pair->idx,pair->n_em_iter);
		fprintf(stdout,"%f,", d);
		IO::print::Array(stdout,pair->SFS, 3, 3, ',');
#endif

	}while(d>tole);

	pair->d=d;
	
	return 0;
}


// sum [ log [ sum [ M(g1,g2) * GL(g1) * GL(g2) ] ]
// void test_em(double **lngls3, double GLSFS[3][3],int **GTSFS,int i1,int i2,char *id1, char*id2, int pair_idx,size_t nSites, int snSites, FILE *outfile){
// void test_em(double **lngls3, double *GLSFS,int **GTSFS,int i1,int i2,char *id1, char*id2, int pair_idx,size_t nSites, int snSites, FILE *outfile){
//
	// double glres=0.0;
	// double gtres=0.0;
	// for(size_t s=0;s<nSites;s++){
		// if ((lngls3[s][(3*i1)+0]==NEG_INF) && (lngls3[s][(3*i1)+1]==NEG_INF) && (lngls3[s][(3*i1)+2]==NEG_INF)){
			// continue;
		// }else if ((lngls3[s][(3*i2)+0]==NEG_INF) && (lngls3[s][(3*i2)+1]==NEG_INF) && (lngls3[s][(3*i2)+2]==NEG_INF)){
			// continue;
		// }else{
			// double glresi=0.0;
			// double gtresi=0.0;
			// for(int i=0;i<3;i++){
				// for(int j=0;j<3;j++){
					// glresi += GLSFS[i][j] * exp(lngls3[s][(3*i1)+i]) * exp(lngls3[s][(3*i2)+j]);
					// gtresi += ((double) GTSFS[pair_idx][(3*i)+j]/(double) snSites)  * exp(lngls3[s][(3*i1)+i]) * exp(lngls3[s][(3*i2)+j]);
					// // fprintf(stderr,"->->->%d %d\n",GTSFS[pair_idx][(3*i)+j],snSites);
				// }
			// }
			// // fprintf(stderr,"->->%f,%f\n",glresi,gtresi);
			// glres+=log(glresi);
			// gtres+=log(gtresi);
		// }
	// }
//
	// fprintf(outfile,"%s,%s,%f,%f\n",id1,id2,glres,gtres);
// }
//


//got some help from https://github.com/lz398/distAngsd/blob/main/vcftest.cpp
// double EM_2DSFS_GL10(double **lngl, double SFS[10][10], int i1, int i2, size_t nSites, double tole){
//
	// double temp;
//
	// double sum;
	// double d;
//
	// memset(SFS,0.01,10*10*sizeof(double));
//
	// do{
//
		// double TMP[10][10];
		// double ESFS[10][10];
		// memset(ESFS,0.0,10*10*sizeof(double));
//
		// for(size_t s=0; s<nSites; s++){
			// sum=0.0;
			// for(int i=0;i<10;i++){
				// for(int j=0;j<10;j++){
					// // SFS * ind1 * ind2
					// TMP[i][j]=SFS[i][j]*lngl[s][(10*i1)+i]*lngl[s][(10*i2)+j];
					// sum += TMP[i][j];
//
				// }
			// }
//
			// for(int i=0;i<10;i++){
				// for(int j=0;j<10;j++){
					// ESFS[i][j] += TMP[i][j]/sum;
				// }
			// }
//
		// }
		// d=0.0;
		// for(int i=0;i<10;i++){
			// for(int j=0;j<10;j++){
				// temp=ESFS[i][j]/(double)nSites;
				// d += fabs(temp-SFS[i][j]);
				// SFS[i][j]=temp;
//
			// }
		// }
//
	// }while(d>tole);
//
	// return d;
// }
