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
#include "em.h"
#include "math_utils.h"
#include "vcf_utils.h"
#include "amova.h"

using size_t=decltype(sizeof(int));


int main(int argc, char **argv) {


	if(argc==1){
		usage(stderr);
		exit(0);
	}

	argStruct *args=argStruct_get(--argc,++argv);

	if(args!=0){



		FILE *in_mtd_ff=NULL;

		paramStruct *pars = paramStruct_init(args);

		IO::outFilesStruct *OUTS = NULL;
		OUTS = new IO::outFilesStruct(args);

		IO::outputStruct* out_emtest_fs=NULL;
		IO::outputStruct* out_dm_fs=NULL;


		if(args->doTest==1){
			out_emtest_fs=OUTS->out_emtest_fs;
		}
		if(args->printMatrix==1){
			out_dm_fs=OUTS->out_dm_fs;
		}

		IO::outputStruct* out_sfs_fs=OUTS->out_sfs_fs;
		IO::outputStruct* out_amova_fs=OUTS->out_amova_fs;


		char *DATETIME=pars->DATETIME;
		DATETIME=get_time();
		fprintf(stderr,"\n%s",DATETIME);

		//TODO print based on analysis type, and to args file
		fprintf(stderr,"\nngsAMOVA -doAMOVA %d -doTest %d -in %s -out %s -isSim %d -minInd %d -printMatrix %d -m %s -doDist %d -maxIter %d -nThreads %d",args->doAMOVA,args->doTest,args->in_fn,args->out_fp,args->isSim,args->minInd,args->printMatrix,args->in_mtd_fn,args->doDist, args->mEmIter, args->mThreads);


		if(args->doAMOVA == -1 || args->doAMOVA == 3 || args->doAMOVA==1){
			fprintf(stderr," -tole %e ",args->tole);
		}
		fprintf(stderr,"\n");

		vcfFile *in_ff = bcf_open(args->in_fn, "r");

		ASSERT(in_ff);

		bcf_hdr_t *hdr = bcf_hdr_read(in_ff);
		bcf1_t *bcf = bcf_init();

		ASSERT(bcf);


		fprintf(stderr, "\n\t-> Reading file: %s\n", args->in_fn);

		pars->nInd=bcf_hdr_nsamples(hdr);
		fprintf(stderr, "\nNumber of samples: %i", pars->nInd);

		const int nContigs=hdr->n[BCF_DT_CTG];
		fprintf(stderr,	"\nNumber of contigs: %d",nContigs);



		for(int ci=0; ci < nContigs; ci++){
			// const bcf_idpair_t& idPair=hdr->id[BCF_DT_CTG][ci];
			// const int contigSize=idPair.val->info[0];
			const int contigSize=hdr->id[BCF_DT_CTG][ci].val->info[0];
			fprintf(stderr,"\nContig %d length:%d\n",ci,contigSize);

			int nBlocks=0;
			if(args->blockSize!=0){
				if(args->blockSize < contigSize)
					nBlocks=(contigSize / args->blockSize) + 1;
				fprintf(stderr,"\nContig %d length:%d nBlocks: %d\n",ci,contigSize,nBlocks);
				for(int bi=0; bi<nBlocks; bi++){
					int blockStart=bi * args->blockSize;
					fprintf(stderr,"\nBlock %d starts at %d\n",bi,blockStart);
					//TODO collect these as pointers to start locations in lngl
				}
			}

		}


		//TODO do below properly
		DATA::sampleStruct *SAMPLES;
		DATA::sampleStruct Samples;
		SAMPLES=&Samples;
		DATA::metadataStruct *MTD;
		DATA::metadataStruct Metadata;
		MTD=&Metadata;


		if(args->in_mtd_fn!=NULL){
			//// BEGIN Read metadata

			in_mtd_ff=IO::getFILE(args->in_mtd_fn,"r");

			//sep can be \t \whitespace or comma
			const char* delims="\t ,\n";


			ASSERT(IO::readFILE::METADATA(MTD, in_mtd_ff, args->whichCol, delims,SAMPLES)==0);
			//// END Read metadata


			if( pars->nInd != MTD->nInds_total){
				fprintf(stderr,"\n[ERROR]: Number of samples in input file (%i) is not equal to number of samples in metadata file (%i); will exit!\n\n", pars->nInd, MTD->nInds_total);
				exit(1);
			}

		}



		if(args->minInd==pars->nInd){
			fprintf(stderr,"\n\t-> -minInd %d is equal to the number of individuals found in file: %d. Setting -minInd to 0 (all).\n",args->minInd, pars->nInd);
			args->minInd=0;
		}

		if(pars->nInd==1){
			fprintf(stderr,"\n\n[ERROR]\tOnly one sample; will exit\n\n");
			free(args->in_fn);
			args->in_fn=NULL;
			free(args);

			paramStruct_destroy(pars);

			exit(1);
		}

		if(pars->nInd<args->minInd){
			fprintf(stderr,"\n\n[ERROR]\tMinimum number of individuals -minInd is set to %d, but input file contains %d individuals; will exit!\n\n",args->minInd,pars->nInd);

			//TODO nice exit?
			// free(args->in_fn);
			// args->in_fn=NULL;
			// free(args);
			// paramStruct_destroy(pars);

			exit(1);
		}

		//TODO
		//print args nicely to an arg file and use -out
		//but smart -out. if output prefix is given detect and truncate
		//so it would play nicely with smk pipelines



		int buf_size=1024;


		/*
		 * lngls[nSites][nInd*10*double]
		 */


		double **lngl=0;
		if(args->doAMOVA==-1 ||args->doAMOVA==1||args->doAMOVA==3){
			lngl=(double**) malloc(buf_size*sizeof(double*));
		}

		//TODO how to make lookup table better
		pars->n_ind_cmb=nCk(pars->nInd,2);



		pars->LUT_indPair_idx=(int **)malloc(pars->nInd * sizeof (int*));

		for (int mi=0;mi<pars->nInd;mi++){
			pars->LUT_indPair_idx[mi]=(int*) malloc( pars->nInd * sizeof(int)) ;
		}

		prepare_LUT_indPair_idx(pars->nInd, pars->LUT_indPair_idx);


		fprintf(stderr,"\nNumber of individual pairs: %d\n",pars->n_ind_cmb);

		//TODO only create if doGeno etc
		/*
		 * SFS_GT3[n_pairs][9+1]
		 * last element contains total number of sites shared
		 */
		int **SFS_GT3;
		if(args->doAMOVA==2 || args->doAMOVA==3){
			SFS_GT3=(int**) malloc(pars->n_ind_cmb*sizeof(int*));
			for (int i=0;i<pars->n_ind_cmb;i++){
				//9 categories per individual pair
				SFS_GT3[i]=(int*)malloc(10*sizeof(int));
				for (int y=0; y<10; y++){
					SFS_GT3[i][y]=0;
				}
			}
		}


		//Pairwise distance matrix
		double *M_PWD_GL=NULL;
		double *M_PWD_GT=NULL;


		// TODO use a lookup table instead
		if(args->doDist!=-1){

			switch(args->doAMOVA){
				case -1:
					M_PWD_GL=(double*) malloc(pars->n_ind_cmb*sizeof(double));
					for(int i=0;i<pars->n_ind_cmb;i++){
						M_PWD_GL[i]=0.0;
					}
					break;
				case 1:
					M_PWD_GL=(double*) malloc(pars->n_ind_cmb*sizeof(double));
					for(int i=0;i<pars->n_ind_cmb;i++){
						M_PWD_GL[i]=0.0;
					}
					break;
				case 2:
					M_PWD_GT=(double*) malloc(pars->n_ind_cmb*sizeof(double));
					for(int i=0;i<pars->n_ind_cmb;i++){
						M_PWD_GT[i]=0.0;
					}
					break;
				case 3:
					M_PWD_GL=(double*) malloc(pars->n_ind_cmb*sizeof(double));
					for(int i=0;i<pars->n_ind_cmb;i++){
						M_PWD_GL[i]=0.0;
					}
					M_PWD_GT=(double*) malloc(pars->n_ind_cmb*sizeof(double));
					for(int i=0;i<pars->n_ind_cmb;i++){
						M_PWD_GT[i]=0.0;
					}
					break;
			}
		}


		/*
		 * [START] Reading sites
		 *
		 * hdr	header
		 * bcf	struct for storing each record
		 *
		 *
		 *
		 * nSites=0
		 * totSites=0
		 *
		 *
		 * if minInd is set
		 * nSites==where minInd threshold can be passed
		 * totSites==all sites processed
		 *
		 *
		 * lngl
		 *
		 *
		 */



		DATA::pairStruct* PAIRS[pars->n_ind_cmb];

		for(int i1=0;i1<pars->nInd-1;i1++){
			for(int i2=i1+1;i2<pars->nInd;i2++){
				int pidx=pars->LUT_indPair_idx[i1][i2];

				PAIRS[pidx] = new DATA::pairStruct(i1,i2,pidx);
			}
		}


		while (bcf_read(in_ff, hdr, bcf) == 0) {



			if(bcf->rlen != 1){
				fprintf(stderr,"\n[ERROR](File reading)\tVCF file REF allele with length of %ld is currently not supported, will exit!\n\n", bcf->rlen);
				exit(1);
			}



			while((int) pars->nSites >= buf_size){

#if 0
				fprintf(stderr,"\n\nrealloc %d at site %d\n",buf_size,pars->nSites);
				fprintf(stderr,"new size: %d\n",buf_size);
#endif
				buf_size=buf_size*2;


				if(args->doAMOVA==-1 ||args->doAMOVA==1||args->doAMOVA==3){
					lngl=(double**) realloc(lngl, buf_size * sizeof(*lngl) );
				}

			}




			if(args->doAMOVA==-1 ||args->doAMOVA==1||args->doAMOVA==3){

				if(VCF::read_GL10_to_GL3(hdr,bcf,lngl,pars,args, pars->nSites, PAIRS)!=0){
					fprintf(stderr,"\n->\tSkipping site %lu for all individuals\n\n",pars->totSites);


					pars->totSites++;

					//next loop will skip this and use the same site_i
					free(lngl[pars->nSites]);
					lngl[pars->nSites]=NULL;

					continue;
				}
			}

			//if doAMOVA==3 and site is skipped for gle; it will be skipped for gt, too

			if(args->doAMOVA==2 || args->doAMOVA ==3 ){

				if(args->gl2gt==1){
					ASSERT(VCF::GT_to_i2i_SFS(hdr,bcf,SFS_GT3,pars,args)==0);
				}else if(args->gl2gt<0){
					ASSERT(VCF::GT_to_i2i_SFS(hdr,bcf,SFS_GT3,pars,args)==0);
				}else{
					exit(1);
				}

			}


			pars->nSites++;
			pars->totSites++;




#if 0
			fprintf(stderr,"\n\n\t-> Printing at site %d\n\n",pars->nSites);
			fprintf(stderr,"%d\n\n",pars->nSites);
			fprintf(stderr,"\r\t-> Printing at site %lu",pars->nSites);
			fprintf(stderr,"\nPrinting at (idx: %lu, pos: %lu 1based %lu) totSites:%d\n\n",pars->nSites,bcf->pos,bcf->pos+1,pars->totSites);
#endif


		}
		fprintf(stderr,"\n\t-> Finished reading sites\n");
		// [END] Reading sites


		pthread_t pairThreads[pars->n_ind_cmb];
		threadStruct* PTHREADS[pars->n_ind_cmb];


		for(int i1=0;i1<pars->nInd-1;i1++){
			for(int i2=i1+1;i2<pars->nInd;i2++){
				int pidx=pars->LUT_indPair_idx[i1][i2];

				PTHREADS[pidx] = new threadStruct(PAIRS[pidx], lngl, pars->nSites, out_sfs_fs, args->tole, args->mEmIter);
			}
		}


		if(args->doAMOVA!=-1)
			fprintf(out_sfs_fs->ff,"Method,Ind1,Ind2,A,D,G,B,E,H,C,F,I,n_em_iter,shared_nSites,Delta,Tole,Sij,Fij,Fij2,IBS0,IBS1,IBS2,R0,R1,Kin\n");


		int nJobs_sent=0;

		if(args->doAMOVA==-1 ||args->doAMOVA==1||args->doAMOVA==3){


			for(int i1=0;i1<pars->nInd-1;i1++){
				for(int i2=i1+1;i2<pars->nInd;i2++){


					int pidx=pars->LUT_indPair_idx[i1][i2];

					if(args->mThreads>1){

						if(nJobs_sent==args->mThreads){
							int t=0;
							while(nJobs_sent>0){
								t=pidx-nJobs_sent;

								if(pthread_join( pairThreads[t], NULL) != 0){
									fprintf(stderr,"\n[ERROR] Problem with joining thread.\n");
									exit(1);
								}else{
									nJobs_sent--;
								}
							}
						}
					}


					if(args->mThreads>1){
						if(pthread_create( &pairThreads[pidx], NULL,t_EM_2DSFS_GL3, PTHREADS[pidx]) == 0){
							nJobs_sent++;
						}else{
							fprintf(stderr,"\n[ERROR] Problem with spawning thread.\n");
							exit(1);
						}
					}else{
						pthread_t pairThread;
						if(pthread_create( &pairThread, NULL,t_EM_2DSFS_GL3, PTHREADS[pidx]) == 0){
							if(pthread_join( pairThread, NULL) != 0){
								fprintf(stderr,"\n[ERROR] Problem with joining thread.\n");
								exit(1);
							}
						}else{
							fprintf(stderr,"\n[ERROR] Problem with spawning thread.\n");
							exit(1);
						}

					}
				}

			}
			//finished indPair loop

			int t=0;
			while(nJobs_sent>0){
				t=pars->n_ind_cmb-nJobs_sent;
				if(pthread_join( pairThreads[t], NULL) != 0){
					fprintf(stderr,"\n[ERROR] Problem with joining thread.\n");
					exit(1);
				}else{
					nJobs_sent--;
				}
			}


			for (int pidx=0; pidx<pars->n_ind_cmb; pidx++){

				DATA::pairStruct* pair=PTHREADS[pidx]->pair;

				ASSERT(pair->snSites>0);

				if(args->doDist==0){
					// M_PWD_GL[pidx]=MATH::EST::Sij(pair->SFS);
					exit(1);
				}else if(args->doDist==1){
					M_PWD_GL[pidx]=(double) (1-MATH::EST::Sij(pair->SFS));
				}else if(args->doDist==2){
					exit(1);
					// M_PWD_GL[pidx]=MATH::EST::Fij(pair->SFS);
				}else{
					exit(1);
				}



				fprintf(out_sfs_fs->ff,"gle,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%ld,%e,%e,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
						hdr->samples[pair->i1],
						hdr->samples[pair->i2],
						pair->snSites * pair->SFS[0],pair->snSites * pair->SFS[1],pair->snSites * pair->SFS[2],
						pair->snSites * pair->SFS[3],pair->snSites * pair->SFS[4],pair->snSites * pair->SFS[5],
						pair->snSites * pair->SFS[6],pair->snSites * pair->SFS[7],pair->snSites * pair->SFS[8],
						pair->n_em_iter,
						pair->snSites,
						pair->d,
						args->tole,
						MATH::EST::Sij(pair->SFS),
						MATH::EST::Fij(pair->SFS),
						MATH::SQUARE(MATH::EST::Fij(pair->SFS)),
						MATH::EST::IBS0(pair->SFS),
						MATH::EST::IBS1(pair->SFS),
						MATH::EST::IBS2(pair->SFS),
						MATH::EST::R0(pair->SFS),
						MATH::EST::R1(pair->SFS),
						MATH::EST::Kin(pair->SFS));

			}



			//TODO maybe join this loop 

		}
		if(args->doAMOVA==2||args->doAMOVA==3){

			for(int i1=0;i1<pars->nInd-1;i1++){
				for(int i2=i1+1;i2<pars->nInd;i2++){

					int snSites=0;

					int pidx=pars->LUT_indPair_idx[i1][i2];

					if(args->doAMOVA==3){

						if((int) PAIRS[pidx]->snSites!=SFS_GT3[pidx][9]){
							fprintf(stderr,"\n[ERROR] Number of shared sites in gle analysis is different than in gt, will exit!\n\n");
							exit(1);
						}
					}

					snSites=SFS_GT3[pidx][9];


					if(args->doDist==0){
						// M_PWD_GT[pidx]=(double) MATH::EST::Sij(SFS_GT3[pidx], snSites);
						exit(1);
					}else if(args->doDist==1){
						M_PWD_GT[pidx]=(double) (1-MATH::EST::Sij(SFS_GT3[pidx], snSites));
					}else if(args->doDist==2){
						exit(1);
						// M_PWD_GT[pidx]=MATH::EST::Fij(SFS_GT3[pidx], snSites);
					}else{
						exit(1);
					}

#if 1
					fprintf(out_sfs_fs->ff,"gt,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
							hdr->samples[i1],
							hdr->samples[i2],
							SFS_GT3[pidx][0],SFS_GT3[pidx][1],SFS_GT3[pidx][2],
							SFS_GT3[pidx][3],SFS_GT3[pidx][4],SFS_GT3[pidx][5],
							SFS_GT3[pidx][6],SFS_GT3[pidx][7],SFS_GT3[pidx][8],
							"gt",
							snSites,
							"gt",
							"gt",
							MATH::EST::Sij(SFS_GT3[pidx], snSites),
							MATH::EST::Fij(SFS_GT3[pidx], snSites),
							MATH::SQUARE(MATH::EST::Fij(SFS_GT3[pidx], snSites)),
							MATH::EST::IBS0(SFS_GT3[pidx], snSites),
							MATH::EST::IBS1(SFS_GT3[pidx], snSites),
							MATH::EST::IBS2(SFS_GT3[pidx], snSites),
							MATH::EST::R0(SFS_GT3[pidx], snSites),
							MATH::EST::R1(SFS_GT3[pidx], snSites),
							MATH::EST::Kin(SFS_GT3[pidx], snSites));

#endif

				}
			}
		}




#if 0
		//print lookup table
		for(int i1=0;i1<pars->nInd-1;i1++){
			for(int i2=i1+1;i2<pars->nInd;i2++){
				fprintf(stderr,"\n%i %i %i\n",LUT_indPair_idx[i1][i2],i1,i2);
			}
		}
#endif

		if(args->printMatrix==1){

			if(args->doAMOVA==-1 ||args->doAMOVA==1||args->doAMOVA==3){

				fprintf(out_dm_fs->ff,"gl,");
				for (int px=0;px<pars->n_ind_cmb;px++){
					fprintf(out_dm_fs->ff,"%f",M_PWD_GL[px]);
					if(px!=pars->n_ind_cmb-1){
						fprintf(out_dm_fs->ff,",");
					}else{
						fprintf(out_dm_fs->ff,"\n");
					}
				}

			}else if(args->doAMOVA==2||args->doAMOVA==3){

				fprintf(out_dm_fs->ff,"gt,");
				for (int px=0;px<pars->n_ind_cmb;px++){
					fprintf(out_dm_fs->ff,"%f",M_PWD_GT[px]);
					if(px!=pars->n_ind_cmb-1){
						fprintf(out_dm_fs->ff,",");
					}else{
						fprintf(out_dm_fs->ff,"\n");
					}
				}
			}
		}

		switch (args->doAMOVA){

			case 1:
				ASSERT(doAMOVA(pars->n_ind_cmb, pars->nInd, MTD, SAMPLES, out_amova_fs->ff, args->sqDist, M_PWD_GL, pars->LUT_indPair_idx, "gl")==0);
				fprintf(stderr, "\n\t-> Finished running AMOVA\n");
				break;
			case 2:
				ASSERT(doAMOVA(pars->n_ind_cmb, pars->nInd, MTD, SAMPLES, out_amova_fs->ff, args->sqDist, M_PWD_GT, pars->LUT_indPair_idx, "gt")==0);
				fprintf(stderr, "\n\t-> Finished running AMOVA\n");
				break;
			case 3:
				ASSERT(doAMOVA(pars->n_ind_cmb, pars->nInd, MTD, SAMPLES, out_amova_fs->ff, args->sqDist, M_PWD_GL, pars->LUT_indPair_idx, "gl")==0);
				ASSERT(doAMOVA(pars->n_ind_cmb, pars->nInd, MTD, SAMPLES, out_amova_fs->ff, args->sqDist, M_PWD_GT, pars->LUT_indPair_idx, "gt")==0);
				fprintf(stderr, "\n\t-> Finished running AMOVA\n");
				break;
			case -1:
				fprintf(stderr, "\n\t-> Skipped running AMOVA\n");
				break;
			default:
				exit(1);
		}


		fprintf(stderr, "Total number of sites processed: %lu\n", pars->totSites);
		fprintf(stderr, "\n");
		fprintf(stderr, "Total number of sites skipped for all individual pairs: %lu\n", pars->totSites-pars->nSites);

		//TODO output is not sorted when threads are used
		//and ind indexes are written instead of ind ids
		fflush(out_sfs_fs->ff);


		bcf_hdr_destroy(hdr);
		bcf_destroy(bcf);

		int BCF_CLOSE;
		if ( (BCF_CLOSE=bcf_close(in_ff))){
			fprintf(stderr,"bcf_close(%s): non-zero status %d\n",args->in_fn,BCF_CLOSE);
			exit(BCF_CLOSE);
		}

		for (int i=0;i<pars->n_ind_cmb;i++){
			delete PAIRS[i];
			delete PTHREADS[i];
		}

		if(args->doAMOVA==1 || args->doAMOVA==3){

			for (size_t s=0;s<pars->nSites;s++){
				free(lngl[s]);
				lngl[s]=NULL;
			}

			free(lngl);


		}

		if(args->doAMOVA==2 || args->doAMOVA==3){

			for (int i=0;i<pars->n_ind_cmb;i++){
				free(SFS_GT3[i]);
				SFS_GT3[i]=NULL;
			}
			free(SFS_GT3);

		}


		if(args->doDist!=-1){
			switch(args->doAMOVA){
				case -1:
					free(M_PWD_GL);
					break;
				case 1:
					free(M_PWD_GL);
					break;
				case 2:
					free(M_PWD_GT);
					break;
				case 3:
					free(M_PWD_GL);
					free(M_PWD_GT);
					break;
			}
		}


		free(args->in_fn);
		args->in_fn=NULL;

		free(args->out_fp);
		args->out_fp=NULL;

		free(args->in_mtd_fn);
		args->in_mtd_fn=NULL;

		free(args);

		paramStruct_destroy(pars);

		if(in_mtd_ff!=NULL){
			fclose(in_mtd_ff);
		}

		delete OUTS;

	}else{
		exit(1);
	}

	return 0;

}

