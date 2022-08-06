#include "em.h"

void* t_EM_2DSFS_GL3(void* p){

	IO::threadStruct* THREAD= (IO::threadStruct*) p;

	if (EM_2DSFS_GL3(THREAD) != 0){
		fprintf(stderr,"\n[ERROR] Problem with EM\n");
		exit(1);
	}

	// fprintf(stderr,"\n\nim a thread!\n\n");


	pthread_exit(NULL);
	

}


int EM_2DSFS_GL3(IO::threadStruct* THREAD){

	//TODO add keepSites if minInd set

	double **lngls=THREAD->lngls;
	size_t nSites=THREAD->total_nSites;
	DATA::pairStruct* pair=THREAD->pair;
	FILE* out_sfs_ff=THREAD->out_sfs_ff;

	double tole=THREAD->tole;
	int doDist=THREAD->doDist;
	double* M_PWD_GL_P=THREAD->M_PWD_GL_PAIR;



	//TODO check underflow
	// fprintf(stderr,"\nEM begin for ind1:%d and ind2:%d \n",pair->i1,pair->i2);

	double temp;
	double sum;
	double d;


	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			pair->SFS[i][j]=(double) 1/ (double) 9;
		}
	}

	do{
#if 0
		fprintf(stderr,"\n");
				for (int x=0;x<3;x++){
					for(int y=0;y<3;y++){
						fprintf(stderr,"%f ",pair->SFS[x][y]);
					}
				}
		fprintf(stderr,"\n");

#endif


		double TMP[3][3];
		double ESFS[3][3];

		for (int i=0; i<3; i++){
			for (int j=0; j<3; j++){
				ESFS[i][j]=0.0;
			}
		}

		// int tme=0;
		for(size_t s=0; s<nSites; s++){


			// skip the sites containing missing values for the individual pair
			if ((lngls[s][(3*pair->i1)+0]==NEG_INF) && (lngls[s][(3*pair->i1)+1]==NEG_INF) && (lngls[s][(3*pair->i1)+2]==NEG_INF)){
				continue;
			}
			if ((lngls[s][(3*pair->i2)+0]==NEG_INF) && (lngls[s][(3*pair->i2)+1]==NEG_INF) && (lngls[s][(3*pair->i2)+2]==NEG_INF)){
				continue;
			}

			sum=0.0;

			// SFS * ind1 * ind2
			//lngls3 (anc,anc),(anc,der),(der,der)
			for(int idx1=0;idx1<3;idx1++){
				for(int idx2=0;idx2<3;idx2++){
					TMP[idx1][idx2]= pair->SFS[idx1][idx2] * exp(lngls[s][(3*pair->i1)+idx1]) * exp(lngls[s][(3*pair->i2)+idx2]);
					sum += TMP[idx1][idx2];
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
				temp=ESFS[i][j]/(double)pair->shared_nSites;
				d += fabs(temp - pair->SFS[i][j]);
				pair->SFS[i][j]=temp;

			}
		}


		pair->n_em_iter++;
#if 0
		fprintf(stderr,"iter: %d d:%f \n",pair->n_em_iter,d);
		fprintf(stderr,"d: %e tole:%e \n",d,tole);
#endif
	}while(d>tole);

	pair->d=d;


	if(doDist==0){
		*M_PWD_GL_P=MATH::EST::Sij(pair->SFS);
	}else if(doDist==1){
		*M_PWD_GL_P=(double) (1-MATH::EST::Sij(pair->SFS));
	}else if(doDist==2){
		*M_PWD_GL_P=MATH::EST::Fij(pair->SFS);
	}else{
		exit(1);
	}


			fprintf(out_sfs_ff,"gle,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%e,%e,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
					pair->i1,
					pair->i2,
					pair->shared_nSites * pair->SFS[0][0],pair->shared_nSites * pair->SFS[0][1],pair->shared_nSites * pair->SFS[0][2],
					pair->shared_nSites * pair->SFS[1][0],pair->shared_nSites * pair->SFS[1][1],pair->shared_nSites * pair->SFS[1][2],
					pair->shared_nSites * pair->SFS[2][0],pair->shared_nSites * pair->SFS[2][1],pair->shared_nSites * pair->SFS[2][2],
					pair->n_em_iter,
					pair->shared_nSites,
					pair->d,
					tole,
					MATH::EST::Sij(pair->SFS),
					MATH::EST::Fij(pair->SFS),
					MATH::SQUARE(MATH::EST::Fij(pair->SFS)),
					MATH::EST::IBS0(pair->SFS),
					MATH::EST::IBS1(pair->SFS),
					MATH::EST::IBS2(pair->SFS),
					MATH::EST::R0(pair->SFS),
					MATH::EST::R1(pair->SFS),
					MATH::EST::Kin(pair->SFS));


	
	return 0;
}


// sum [ log [ sum [ M(g1,g2) * GL(g1) * GL(g2) ] ]
void test_em(double **lngls3, double GLSFS[3][3],int **GTSFS,int i1,int i2,char *id1, char*id2, int pair_idx,size_t nSites, int shared_nSites, FILE *outfile){

	double glres=0.0;
	double gtres=0.0;
	for(size_t s=0;s<nSites;s++){
		if ((lngls3[s][(3*i1)+0]==NEG_INF) && (lngls3[s][(3*i1)+1]==NEG_INF) && (lngls3[s][(3*i1)+2]==NEG_INF)){
			continue;
		}else if ((lngls3[s][(3*i2)+0]==NEG_INF) && (lngls3[s][(3*i2)+1]==NEG_INF) && (lngls3[s][(3*i2)+2]==NEG_INF)){
			continue;
		}else{
			double glresi=0.0;
			double gtresi=0.0;
			for(int i=0;i<3;i++){
				for(int j=0;j<3;j++){
					glresi += GLSFS[i][j] * exp(lngls3[s][(3*i1)+i]) * exp(lngls3[s][(3*i2)+j]);
					gtresi += ((double) GTSFS[pair_idx][(3*i)+j]/(double) shared_nSites)  * exp(lngls3[s][(3*i1)+i]) * exp(lngls3[s][(3*i2)+j]);
					// fprintf(stderr,"->->->%d %d\n",GTSFS[pair_idx][(3*i)+j],shared_nSites);
				}
			}
			// fprintf(stderr,"->->%f,%f\n",glresi,gtresi);
			glres+=log(glresi);
			gtres+=log(gtresi);
		}
	}

	fprintf(outfile,"%s,%s,%f,%f\n",id1,id2,glres,gtres);
	// fprintf(stderr,"%s,%s,%f,%f\n",id1,id2,glres,gtres);
}



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
