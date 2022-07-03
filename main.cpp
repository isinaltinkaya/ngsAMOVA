/*
 *
 */


#include <stdio.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

#include <inttypes.h>
#include <math.h>

#include <limits>

#include <time.h>

#include "main.h"
#include "em.h"
#include "math_utils.h"
#include "vcf_utils.h"
#include "io.h"
#include "shared.h"


using size_t=decltype(sizeof(int));

const double NEG_INF = -std::numeric_limits<double>::infinity();

FILE *getFILE(const char*fname,const char* mode){
	FILE *fp;
	if(NULL==(fp=fopen(fname,mode))){
		fprintf(stderr,"[%s:%s()]\t->Error opening FILE handle for file:%s exiting\n",__FILE__,__FUNCTION__,fname);
		exit(1);
	}
	return fp;
}


void get_gt_sfs( int* gtdata, int **sfs, int32_t *ptr1, int32_t *ptr2, int pair_idx){

	int gti1=0;
	int gti2=0;

	//using binary input genotypes from VCF GT tag
	//assume ploidy=2
	for (int i=0; i<2;i++){
		gti1 += bcf_gt_allele(ptr1[i]);
		gti2 += bcf_gt_allele(ptr2[i]);
	}
	// fprintf(stderr,"\n%d:%d %d:%d -> %d\n",i1,gti1,i2,gti2,get_3x3_idx[gti1][gti2]);
	sfs[pair_idx][get_3x3_idx[gti1][gti2]]++;
}



int main(int argc, char **argv) {


	if(argc==1){
		usage();
		exit(0);
	}

	argStruct *args=argStruct_get(--argc,++argv);

	if(args!=NULL){

		paramStruct *pars = paramStruct_init();

		int nInd=pars->nInd;
		size_t nSites=pars->nSites;
		char *anc=pars->anc;
		char *der=pars->der;
		char *in_fn=args->in_fn;

		vcfFile * in_ff = bcf_open(in_fn, "r");

		if (in_ff == NULL) {
			exit(1);
		}

		if (bcf == 0) {
			exit(1);
		}

		bcf_hdr_t *hdr = bcf_hdr_read(in_ff);
		bcf1_t *bcf = bcf_init();

		nInd=bcf_hdr_nsamples(hdr);

		if(nInd==1){
			fprintf(stderr,"\n\n[ERROR]\tOnly one sample; will exit\n\n");
			exit(0);
		}

		//TODO
		//print args nicely to an arg file and use -out
		//but smart -out. if output prefix is given detect and truncate
		//so it would play nicely with smk pipelines
		nSites=0;


		char *DATETIME=pars->DATETIME;
		DATETIME=get_time();

		fprintf(stderr,"\n%s",DATETIME);
		fprintf(stderr,"\n./ngsAMOVA -in %s -tole %e -isSim %d -onlyShared %d\n",args->in_fn,args->tole,args->isSim,args->onlyShared);

		fprintf(stderr, "\nReading file: \"%s\"", in_fn);
		fprintf(stderr, "\nNumber of samples: %i", bcf_hdr_nsamples(hdr));
		fprintf(stderr,	"\nNumber of chromosomes: %d",hdr->n[BCF_DT_CTG]);



		int buf_size=1024;

		/*
		 * lngls[nSites][nInd*10*double]
		 */
		double **lngls=0;

		lngls=(double**) malloc(buf_size*sizeof(double*));
		// for (int i=0;i<buf_size;i++){
			// lngls[i]=(double*)malloc(nInd*10*sizeof(double));
		// }

		//TODO build lookup table

		int n_ind_cmb=nCk(nInd,2);
		
		fprintf(stderr,"\nNumber of individual pairs: %d\n",n_ind_cmb);

		//TODO only create if doGeno etc
		/*
		 * gt_sfs[n_pairs][9]
		 */
		int **gt_sfs;
		gt_sfs=(int**) malloc(n_ind_cmb*sizeof(int*));

		for (int i=0;i<n_ind_cmb;i++){
			//9 categories per individual pair
			gt_sfs[i]=(int*)malloc(9*sizeof(int));
		}

		for(int x=0;x<n_ind_cmb;x++){
			for (int y=0; y<9; y++){
				gt_sfs[x][y]=0;
			}
		}


		if(args->isSim==1){
			anc=(char*)malloc(buf_size*sizeof(char));
			der=(char*)malloc(buf_size*sizeof(char));
		}

		const char* TAG=NULL;

		/*
		 * [START] Reading sites
		 *
		 * hdr	header
		 * bcf	struct for storing each record
		 */

		while (bcf_read(in_ff, hdr, bcf) == 0) {
//
			// TAG="DP";
			// get_data<int32_t> dp;
			// dp.n = bcf_get_format_int32(hdr,bcf,TAG,&dp.data,&dp.size_e);
			// if(dp.n<0){
				// fprintf(stderr,"\n[ERROR](File reading)\tVCF tag \"%s\" does not exist; will exit!\n\n",TAG);
				// exit(1);
			// }
//
			// //TODO use only sites non-missing for all individuals for now
			// for(int indi=0; indi<nInd; indi++){
				// if(dp.data[indi]==0){
					// // fprintf(stderr,"\nDepth: %d at (site_i: %d, pos: %d) in (ind_i: %d, ind_id: %s)\n",dp.data[indi],nSites,bcf->pos,indi,hdr->samples[indi]);
					// dp.n_missing++;
				// }
			// }
			// if (dp.n_missing>0) continue;
//
//
			get_data<float> lgl;

			TAG="GL";
			lgl.n = bcf_get_format_float(hdr,bcf,TAG,&lgl.data,&lgl.size_e);

			if(lgl.n<0){
				fprintf(stderr,"\n[ERROR](File reading)\tVCF tag \"%s\" does not exist; will exit!\n\n",TAG);
				exit(1);
			}



			while(nSites>=buf_size){

				// fprintf(stderr,"\n\nrealloc %d at site %d\n",buf_size,nSites);
				buf_size=buf_size*2;
				// fprintf(stderr,"new size: %d\n",buf_size);

				lngls=(double**) realloc(lngls, buf_size * sizeof(*lngls) );

				if(args->isSim==1){
					anc=(char*)realloc(anc,buf_size*sizeof(char));
					der=(char*)realloc(der,buf_size*sizeof(char));
				}

			}

			lngls[nSites]=(double*)malloc(nInd*10*sizeof(double));

			//TODO
			//- only store lngls3 if majmin 3x3 will be used
			//- check below what are you


			for(int indi=0; indi<nInd; indi++){
				int n_missing_gl_indi=0;
				for(int i=0;i<10;i++){
					lngls[nSites][(10*indi)+i]=NEG_INF;
					if(isnan(lgl.data[(10*indi)+i])){
						// fprintf(stderr,"\nMissing data\n");
						n_missing_gl_indi++;
					}else if (bcf_float_is_vector_end(lgl.data[(10*indi)+i])){
						fprintf(stderr,"\nwhat are you\n");
						// lgl.n_missing++;
					// }else if(bcf_float_is_missing(lgl.data[i])){
						// fprintf(stderr,"\nMissing data\n");
						// lgl.n_missing++;
					// }else if(std::isnan(lgl.data[(10*indi)+i])){
					}else{
						// fprintf(stderr,"\n->->i:%d ind:%d %s %f -> %f\n",i,indi,hdr->samples[indi],lgl.data[(10*indi)+i],(double)log2ln(lgl.data[(10*indi)+i]));
						lngls[nSites][(10*indi)+i]=(double)log2ln(lgl.data[(10*indi)+i]);
					}
				}
				
				//if all 10 gls missing for ind
				if (n_missing_gl_indi==10){

					lgl.n_missing_ind++;
					
				}
			}

			if(args->onlyShared==1){
				//if all 10 gls missing for at least 1 ind
				if (lgl.n_missing_ind>0){
					free(lngls[nSites]);
					lngls[nSites]=NULL;
					continue;
				}
			}

			//skip site if missing for all individuals
			if (lgl.n_missing_ind == nInd) {
				// fprintf(stderr,"->->n_missing_ind: %d, nInd: %d\n",lgl.n_missing_ind,nInd);
				free(lngls[nSites]);
				lngls[nSites]=NULL;
				continue;

			//if there are only 2 individuals, skip site regardless of onlyShared val
			}else if (nInd==2){
				if(lgl.n_missing_ind>0){
					free(lngls[nSites]);
					lngls[nSites]=NULL;
					continue;
				}
			}

			// fprintf(stderr,"->->->n_missing_ind: %d, nInd: %d\n",lgl.n_missing_ind,nInd);


			get_data<int32_t> gt;
			gt.n = bcf_get_genotypes(hdr,bcf,&gt.data,&gt.size_e);
			gt.ploidy=gt.n/nInd;

			if(gt.n<0){
				fprintf(stderr,"\n[ERROR](File reading)\tProblem with reading GT; will exit!\n\n");
				exit(1);
			}

			if (gt.ploidy!=2){
				fprintf(stderr,"ERROR:\n\nploidy: %d not supported\n",gt.ploidy);
				return 1;
			}



			for(int i1=0;i1<nInd-1;i1++){
				for(int i2=i1+1;i2<nInd;i2++){
					
					int32_t *ptr1 = gt.data + i1*gt.ploidy;
					int32_t *ptr2 = gt.data + i2*gt.ploidy;
					int pair_idx=nCk_idx(nInd,i1,i2);
					// fprintf(stderr,"\n->i1: %d, i2: %d, pair: %d, nInd: %d\n\n", i1, i2, pair_idx,nInd);
					// fprintf(stderr,"%d\t%d\t%d\n", i1, i2, pair_idx);

					get_gt_sfs(gt.data,gt_sfs,ptr1,ptr2,pair_idx);
					
				}
			}


			if(bcf_is_snp(bcf)){
				if(bcf->rlen == 1){

					if(args->isSim==1){
						anc[nSites]=bcf_allele_charToInt[bcf->d.allele[0][0]];
						der[nSites]=bcf_allele_charToInt[bcf->d.allele[1][0]];
						// fprintf(stderr,"\n->->gt %d \n",bcf_alleles_get_gtidx(bcf->d.allele[0][0],bcf->d.allele[1][0]));
						// fprintf(stderr,"\n->->gt %d \n",bcf_alleles_get_gtidx(bcf_allele_charToInt[bcf->d.allele[0][0]],bcf_allele_charToInt[bcf->d.allele[1][0]]));
					}

				}else{
					fprintf(stderr,"\n[ERROR](File reading)\tVCF file REF allele with length of %d is currently not supported, will exit!\n\n", bcf->rlen);
					exit(1);
				}
			}

			//TODO bcf_hdr_set_samples efficient sample parsing
			// if (args->doInd==1){
			// int i1=args->ind1;
			// int i2=args->ind2;
			// get_data<int32_t> gt;
			// gt.n = bcf_get_genotypes(hdr,bcf,&gt.data,&gt.size_e);
			// }

			// fprintf(stderr,"\n\n\t-> Printing at site %d\n\n",nSites);
			// fprintf(stderr,"%d\n\n",nSites);


			nSites++;

		}
		// [END] Reading sites


#if 0
		for(int s=0;s<nSites;s++){
			int i1=0;
			fprintf(stderr,"\n-> site: %d anc:%d der:%d gtidx ancanc:%d ancder:%d derder:%d",s,anc[s],der[s],bcf_alleles_get_gtidx(anc[s],anc[s]),bcf_alleles_get_gtidx(anc[s],der[s]),bcf_alleles_get_gtidx(der[s],der[s]));
			fprintf(stderr,"\n-> ind:%s (%f",hdr->samples[i1],lngls[s][(10*i1)+bcf_alleles_get_gtidx(anc[s],anc[s])]);
			fprintf(stderr," %f",lngls[s][(10*i1)+bcf_alleles_get_gtidx(anc[s],der[s])]);
			fprintf(stderr," %f",lngls[s][(10*i1)+bcf_alleles_get_gtidx(der[s],der[s])]);
			fprintf(stderr,")\n");
		}
#endif

		
		fprintf(stdout,"Method,Ind1,Ind2,A,D,G,B,E,H,C,F,I,N_EM_ITER,nSites\n");

		for(int i1=0;i1<nInd-1;i1++){
			for(int i2=i1+1;i2<nInd;i2++){


				int shared_nSites=0;

				if(args->onlyShared==1){
					shared_nSites=nSites;
				}else{


					for(size_t s=0; s<nSites; s++){

						int n_neg_inf_i1=0;
						int n_neg_inf_i2=0;

						if (lngls[s][(10*i1)+s]==NEG_INF){
							n_neg_inf_i1++;
							// fprintf(stderr,"\nNEG INF\n");
						}
						if (lngls[s][(10*i2)+s]==NEG_INF){
							n_neg_inf_i2++;
							// fprintf(stderr,"\nNEG INF2\n");
						}
						//skip if all -inf == all missing for i1 or i2
						if(n_neg_inf_i1==10){
							continue;
						}else if(n_neg_inf_i2==10){
							continue;
						}else{
							shared_nSites++;
						}
					}
					
				}
				if(shared_nSites==0){
					continue;
				}

				int pair_idx=nCk_idx(nInd,i1,i2);

				fprintf(stdout,"gt,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%s,%d\n",
						hdr->samples[i1],
						hdr->samples[i2],
						gt_sfs[pair_idx][0],gt_sfs[pair_idx][1],gt_sfs[pair_idx][2],
						gt_sfs[pair_idx][3],gt_sfs[pair_idx][4],gt_sfs[pair_idx][5],
						gt_sfs[pair_idx][6],gt_sfs[pair_idx][7],gt_sfs[pair_idx][8],
						"gt",
						shared_nSites);

				// double SFS2D[10][10];
				// EM_2DSFS_GL10(lngls,SFS2D,i1,i2,nSites,args->tole);
				double SFS2D3[3][3];
				int n_em_iter;
				n_em_iter=EM_2DSFS_GL3(lngls,SFS2D3,i1,i2,nSites,shared_nSites,args->tole,anc,der);
				// print_2DM(3,3,*SFS2D3);
				// fprintf(stdout,"gl,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
						// hdr->samples[i1],
						// hdr->samples[i2],
						// SFS2D3[0][0],SFS2D3[0][1],SFS2D3[0][2],
						// SFS2D3[1][0],SFS2D3[1][1],SFS2D3[1][2],
						// SFS2D3[2][0],SFS2D3[2][1],SFS2D3[2][2]);
//
				// fprintf(stdout,"gle,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
						// hdr->samples[i1],
						// hdr->samples[i2],
						// shared_nSites*SFS2D3[0][0],shared_nSites*SFS2D3[0][1],shared_nSites*SFS2D3[0][2],
						// shared_nSites*SFS2D3[1][0],shared_nSites*SFS2D3[1][1],shared_nSites*SFS2D3[1][2],
						// shared_nSites*SFS2D3[2][0],shared_nSites*SFS2D3[2][1],shared_nSites*SFS2D3[2][2]);
//
				fprintf(stdout,"gle,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d\n",
						hdr->samples[i1],
						hdr->samples[i2],
						shared_nSites*SFS2D3[0][0],shared_nSites*SFS2D3[0][1],shared_nSites*SFS2D3[0][2],
						shared_nSites*SFS2D3[1][0],shared_nSites*SFS2D3[1][1],shared_nSites*SFS2D3[1][2],
						shared_nSites*SFS2D3[2][0],shared_nSites*SFS2D3[2][1],shared_nSites*SFS2D3[2][2],
						n_em_iter,
						shared_nSites);
			}
		}

		fprintf(stderr, "Total number of sites processed: %i\n", nSites);

		bcf_hdr_destroy(hdr);
		bcf_destroy(bcf);

		int BCF_CLOSE;
		if ( (BCF_CLOSE=bcf_close(in_ff))){
			fprintf(stderr,"bcf_close(%s): non-zero status %d\n",in_fn,BCF_CLOSE);
			exit(BCF_CLOSE);
		}

		if(args->isSim==1){
			free(anc);
			anc=NULL;
			free(der);
			der=NULL;
		}

		for (int i=0;i<n_ind_cmb;i++){
			free(gt_sfs[i]);
			gt_sfs[i]=NULL;
		}
		free(gt_sfs);

		for (size_t s=0;s<nSites;s++){
			free(lngls[s]);
			lngls[s]=NULL;
		}

		free(lngls);

		paramStruct_destroy(pars);

	}
	free(args->in_fn);
	free(args->out_fp);
	free(args);

	return 0;

}
