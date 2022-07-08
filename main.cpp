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


int set_lngls3(double **lngls3, float *lgldata,argStruct *args, int nInd,char a1, char a2, int site){

	int n_missing_ind=0;
	for(int indi=0; indi<nInd; indi++){

		for(int ix=0;ix<3;ix++){
			lngls3[site][(3*indi)+ix]=NEG_INF;
		}

		//TODO only checking the first for now
		//what is the expectation in real cases?
		//should we skip sites where at least one is set to missing?
		//
		if(isnan(lgldata[(10*indi)+0])){
			if(args->onlyShared==1){
				return 1;
			}else{
				n_missing_ind++;
			}
		}else{
			// fprintf(stderr,"\n->->gt %d %d\n",bcf_alleles_get_gtidx(bcf->d.allele[0][0],bcf->d.allele[1][0]));
			// fprintf(stderr,"\n->->gt %d \n",bcf_alleles_get_gtidx(bcf_allele_charToInt[bcf->d.allele[0][0]],bcf_allele_charToInt[bcf->d.allele[1][0]]));
			lngls3[site][(3*indi)+0]=(double) log2ln(lgldata[(10*indi)+bcf_alleles_get_gtidx(a1,a1)]);
			lngls3[site][(3*indi)+1]=(double) log2ln(lgldata[(10*indi)+bcf_alleles_get_gtidx(a1,a2)]);
			lngls3[site][(3*indi)+2]=(double) log2ln(lgldata[(10*indi)+bcf_alleles_get_gtidx(a2,a2)]);
			// fprintf(stderr,"\n->->lgl %f %f %f\n",lgldata[(10*indi)+bcf_alleles_get_gtidx(a1,a1)],lgldata[(10*indi)+bcf_alleles_get_gtidx(a1,a2)],lgldata[(10*indi)+bcf_alleles_get_gtidx(a2,a2)]);
			// fprintf(stderr,"\n->->lngls3 %f %f %f\n",lngls3[site][(3*indi)+0],lngls3[site][(3*indi)+1],lngls3[site][(3*indi)+2]);
		}
	}

	//skip site if missing for all individuals
	if (nInd==n_missing_ind){
		return 1;
	}

	// //if there are only 2 individuals, skip site regardless of onlyShared val
	if (nInd==2){
		if(n_missing_ind>0){
			return 1;
		}
	}

	if(args->minInd>0){
		//skip site if minInd is defined and #non-missing inds=<nInd
		if( (nInd - n_missing_ind) < args->minInd ){
			fprintf(stderr,"\n\nMinimum number of individuals -minInd is set to %d, but nInd-n_missing_ind==n_nonmissing_ind is %d at site %d\n\n",args->minInd,nInd-n_missing_ind,site);
			return 1;
		}else{
			return 0;
		}
	}else{
		return 0;
	}

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
		//TODO already exit(0) in usage
		exit(0);
	}

	argStruct *args=argStruct_get(--argc,++argv);

	if(args!=0){

		paramStruct *pars = paramStruct_init();

		int nInd=pars->nInd;
		size_t nSites=pars->nSites;
		char *in_fn=args->in_fn;

		size_t totSites=0;

		vcfFile *in_ff = bcf_open(in_fn, "r");

		if (in_ff == NULL) {
			free(args);
			exit(1);
		}

		bcf_hdr_t *hdr = bcf_hdr_read(in_ff);

		bcf1_t *bcf = bcf_init();

		if (bcf == 0) {
			free(args);
			exit(1);
		}

		nInd=bcf_hdr_nsamples(hdr);
		//
		// if(nInd==1){
		// fprintf(stderr,"\n\n[ERROR]\tOnly one sample; will exit\n\n");
		// free(args->in_fn);
		// args->in_fn=NULL;
		// free(args);
		//
		// paramStruct_destroy(pars);
		//
		// exit(1);
		// }
		if(nInd<args->minInd){
			fprintf(stderr,"\n\n[ERROR]\tMinimum number of individuals -minInd is set to %d, but input file contains %d individuals; will exit!\n\n",args->minInd,nInd);

			free(args->in_fn);
			args->in_fn=NULL;
			free(args);

			paramStruct_destroy(pars);

			exit(1);
		}
		//
		//TODO
		//print args nicely to an arg file and use -out
		//but smart -out. if output prefix is given detect and truncate
		//so it would play nicely with smk pipelines
		nSites=0;


		char *DATETIME=pars->DATETIME;
		DATETIME=get_time();

		fprintf(stderr,"\n%s",DATETIME);
		fprintf(stderr,"\n./ngsAMOVA -in %s -tole %e -isSim %d -onlyShared %d -minInd %d\n",args->in_fn,args->tole,args->isSim,args->onlyShared,args->minInd);

		fprintf(stderr, "\nReading file: \"%s\"", in_fn);
		fprintf(stderr, "\nNumber of samples: %i", bcf_hdr_nsamples(hdr));
		fprintf(stderr,	"\nNumber of chromosomes: %d",hdr->n[BCF_DT_CTG]);



		int buf_size=1024;
		// int buf2_size=1024;


		/*
		 * lngls[nSites][nInd*10*double]
		 */
		// double **lngls=0;
		//
		// lngls=(double**) malloc(buf_size*sizeof(double*));
		// for (int i=0;i<buf_size;i++){
		// lngls[i]=(double*)malloc(nInd*10*sizeof(double));
		// }


		double **lngls3=0;
		lngls3=(double**) malloc(buf_size*sizeof(double*));

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

			totSites++;


			while((int) nSites >= buf_size){

				// fprintf(stderr,"\n\nrealloc %d at site %d\n",buf_size,nSites);
				// fprintf(stderr,"new size: %d\n",buf_size);
				buf_size=buf_size*2;

				// lngls=(double**) realloc(lngls, buf_size * sizeof(*lngls) );
				lngls3=(double**) realloc(lngls3, buf_size * sizeof(*lngls3) );

				// if(args->isSim==1){
				// anc=(char*)realloc(anc,buf_size*sizeof(char));
				// der=(char*)realloc(der,buf_size*sizeof(char));
				// }
			}

			//
			// if(bcf_is_snp(bcf)){
			// if(args->isSim==1){
			// anc[nSites]=bcf_allele_charToInt[(unsigned char) bcf->d.allele[0][0]];
			// der[nSites]=bcf_allele_charToInt[(unsigned char) bcf->d.allele[1][0]];
			// // fprintf(stderr,"\n->->gt %d \n",bcf_alleles_get_gtidx(bcf->d.allele[0][0],bcf->d.allele[1][0]));
			// // fprintf(stderr,"\n->->gt %d \n",bcf_alleles_get_gtidx(bcf_allele_charToInt[bcf->d.allele[0][0]],bcf_allele_charToInt[bcf->d.allele[1][0]]));
			// }
			if(bcf->rlen != 1){
				fprintf(stderr,"\n[ERROR](File reading)\tVCF file REF allele with length of %ld is currently not supported, will exit!\n\n", bcf->rlen);
				exit(1);
			}

			// lngls[nSites]=(double*)malloc(nInd*10*sizeof(double));
			lngls3[nSites]=(double*)malloc(nInd*3*sizeof(double));


			if(bcf_is_snp(bcf)){
				char a1=bcf_allele_charToInt[(unsigned char) bcf->d.allele[0][0]];
				char a2=bcf_allele_charToInt[(unsigned char) bcf->d.allele[1][0]];

				if (set_lngls3(lngls3, lgl.data,args, nInd,a1,a2, nSites)==1){
					free(lngls3[nSites]);
					lngls3[nSites]=NULL;
// fprintf(stderr,"\nset_lngls3 returns 1 at site (idx: %lu, pos: %lu 1based: %lu)\n\n",nSites,bcf->pos,bcf->pos+1);
					continue;
				// }else{
				}
			}else{
				//TODO check
				fprintf(stderr,"\n\nHERE BCF_IS_SNP==0!!!\n\n");
				free(lngls3[nSites]);
				lngls3[nSites]=NULL;
				continue;
			}


			get_data<int32_t> gt;
			gt.n = bcf_get_genotypes(hdr,bcf,&gt.data,&gt.size_e);
			gt.ploidy=gt.n/nInd;

			if(gt.n<0){
				fprintf(stderr,"\n[ERROR](File reading)\tProblem with reading GT; will exit!\n\n");
				exit(1);
			}

			if (gt.ploidy!=2){
				fprintf(stderr,"ERROR:\n\nploidy: %d not supported\n",gt.ploidy);
				exit(1);
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

			//TODO bcf_hdr_set_samples efficient sample parsing
			// if (args->doInd==1){
			// int i1=args->ind1;
			// int i2=args->ind2;
			// get_data<int32_t> gt;
			// gt.n = bcf_get_genotypes(hdr,bcf,&gt.data,&gt.size_e);
			// }

			// fprintf(stderr,"\n\n\t-> Printing at site %d\n\n",nSites);
			// fprintf(stderr,"%d\n\n",nSites);
			// fprintf(stderr,"\r\t-> Printing at site %lu",nSites);
			// fprintf(stderr,"\nPrinting at (idx: %lu, pos: %lu 1based %lu)\n\n",nSites,bcf->pos,bcf->pos+1);

			nSites++;

		}
		fprintf(stderr,"\n");
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


		fprintf(stdout,"Method,Ind1,Ind2,A,D,G,B,E,H,C,F,I,N_EM_ITER,nSites,tole\n");


		for(int i1=0;i1<nInd-1;i1++){
			for(int i2=i1+1;i2<nInd;i2++){

				//TODO avoid checking shared sites multiple times for ind pairs
				//
				int shared_nSites=0;

				if(args->onlyShared==1){

					shared_nSites=nSites;

				}else{
				
					for(size_t s=0; s<nSites; s++){

						if ((lngls3[s][(3*i1)+0]==NEG_INF) && (lngls3[s][(3*i1)+1]==NEG_INF) && (lngls3[s][(3*i1)+2]==NEG_INF)){
	// fprintf(stderr,"\n\nHERE!%d !!\n\n",shared_nSites);
							// fprintf(stderr,"\nNEG INF\n");
							continue;
						}else if ((lngls3[s][(3*i2)+0]==NEG_INF) && (lngls3[s][(3*i2)+1]==NEG_INF) && (lngls3[s][(3*i2)+2]==NEG_INF)){
	// fprintf(stderr,"\n\nHERE2!%d !!\n\n",shared_nSites);
							// fprintf(stderr,"\nNEG INF2\n");
							continue;
						}else{
							shared_nSites++;
						}
					}
				}

				if(shared_nSites==0){
					fprintf(stderr,"\n->No shared sites found for i1:%d i2:%d\n",i1,i2);
					fprintf(stdout,"gt,%s,%s,NA,NA,NA,NA,NA,NA,NA,NA,NA,%s,%d,%e\n",
						hdr->samples[i1],
						hdr->samples[i2],
						"gt",
						shared_nSites,
						args->tole);
					//tole is not used with gt but print to indicate that result is from this specific run

					fprintf(stdout,"gle,%s,%s,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,%d,%e\n",
							hdr->samples[i1],
							hdr->samples[i2],
							shared_nSites,
							args->tole);
					continue;
				}

				int pair_idx=nCk_idx(nInd,i1,i2);

				fprintf(stdout,"gt,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%s,%d,%e\n",
						hdr->samples[i1],
						hdr->samples[i2],
						gt_sfs[pair_idx][0],gt_sfs[pair_idx][1],gt_sfs[pair_idx][2],
						gt_sfs[pair_idx][3],gt_sfs[pair_idx][4],gt_sfs[pair_idx][5],
						gt_sfs[pair_idx][6],gt_sfs[pair_idx][7],gt_sfs[pair_idx][8],
						"gt",
						shared_nSites,
						args->tole);


				double SFS2D3[3][3];
				int n_em_iter;
				n_em_iter=EM_2DSFS_GL3(lngls3,SFS2D3,i1,i2,nSites,shared_nSites,args->tole);

				fprintf(stdout,"gle,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%e\n",
						hdr->samples[i1],
						hdr->samples[i2],
						shared_nSites*SFS2D3[0][0],shared_nSites*SFS2D3[0][1],shared_nSites*SFS2D3[0][2],
						shared_nSites*SFS2D3[1][0],shared_nSites*SFS2D3[1][1],shared_nSites*SFS2D3[1][2],
						shared_nSites*SFS2D3[2][0],shared_nSites*SFS2D3[2][1],shared_nSites*SFS2D3[2][2],
						n_em_iter,
						shared_nSites,
						args->tole);
			}
		}

		fprintf(stderr, "Total number of sites processed: %lu\n", totSites);
		fprintf(stderr, "Total number of sites skipped for all individual pairs: %lu\n", totSites-nSites);

		bcf_hdr_destroy(hdr);
		bcf_destroy(bcf);

		int BCF_CLOSE;
		if ( (BCF_CLOSE=bcf_close(in_ff))){
			fprintf(stderr,"bcf_close(%s): non-zero status %d\n",in_fn,BCF_CLOSE);
			exit(BCF_CLOSE);
		}

		for (int i=0;i<n_ind_cmb;i++){
			free(gt_sfs[i]);
			gt_sfs[i]=NULL;
		}
		free(gt_sfs);

		for (size_t s=0;s<nSites;s++){
			// free(lngls[s]);
			// lngls[s]=NULL;
			free(lngls3[s]);
			lngls3[s]=NULL;
		}

		// free(lngls);
		free(lngls3);

		free(args->in_fn);
		args->in_fn=NULL;
		free(args);
		// free(args->out_fp);


		paramStruct_destroy(pars);

	}else{
		//if args==NULL
		//already freed in io.cpp
		exit(1);
	}

	return 0;

}
