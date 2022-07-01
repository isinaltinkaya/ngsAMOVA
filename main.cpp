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
#include "random_generator.h"
#include "estimator.h"
#include "vcf_utils.h"
#include "io.h"
#include "shared.h"


using size_t=decltype(sizeof(int));

const double NEG_INF = -std::numeric_limits<double>::infinity();


//print 2D matrix
void print_2DM(size_t x1, size_t x2, double *M){
	printf("\n");
	int m,n;
	for (m=0;m<x1;m++){
		for (n=0;n<x2;n++){
			printf("A[%d,%d]=%f\n",m,n,M[(m*x1)+n]);
		}
	}
};

/*
 * Binomial coefficient
 * n choose k
 */
int nCk(int n, int k){
	int res = 1;
	if (k > n-k){
		k = n-k;
	}
	// [ n * (n-1) * ... * (n-k+1) ] / [ k * (k-1) * ... * 1 ]
	for (int i = 0; i < k; ++i) {
		res = res * (n-i) / (i+1) ;
	}
	return res;
}

/*
 * Maps a given pair of objects to their index in
 * the lexicographically ordered binomial coefficients
 * a.k.a. array of pairs
 */
int nCk_idx(int nInd, int i1, int i2){
	// int ret=(nInd - nCk((nInd-i1),2))+(i2-i1)+1;
	int ret=(nCk(nInd,2) - nCk((nInd-i1),2))+(i2-i1)-1;
	// fprintf(stderr,"\n->ret: %d, nInd:%d, %d + %d + 1\n",ret,nInd,nCk((nInd-i1),2),(i2-i1));
	return ret;
	// return (nInd - nCk((nInd-i1),2))+(i2-i1)+1;
}


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

	// int idx=get_3x3_idx[gti1][gti2];
	// sfs[pair_idx][idx]++;
	sfs[pair_idx][get_3x3_idx[gti1][gti2]]++;

}



int main(int argc, char **argv) {


	if(argc==1){
		usage();
		exit(0);
	}

	argStruct *args=args_get(--argc,++argv);
	paramStruct *pars = pars_init();

	int nInd=pars->nInd;
	size_t nSites=pars->nSites;

	char *anc=pars->anc;
	char *der=pars->der;


	if(args!=NULL){

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

		// if(nInd==1){
			// fprintf(stderr,"\n\n[ERROR]\tOnly one sample; will exit\n\n");
			// exit(0);
		// }


		nSites=0;


		fprintf(stderr, "\nReading file: \"%s\"", in_fn);
		fprintf(stderr, "\nNumber of samples: %i", bcf_hdr_nsamples(hdr));
		fprintf(stderr,	"\nNumber of chromosomes: %d",hdr->n[BCF_DT_CTG]);


		//todo first use index of bcf etc to determine nsites
		//then start analyses
		//

		int buf_size=1024;
		// int buf_size=2;

		/*
		 * lngls[nSites][nInd*10*double]
		 */
		double **lngls=0;

		lngls=(double**) malloc(buf_size*sizeof(double*));
		for (int i=0;i<buf_size;i++){
			lngls[i]=(double*)malloc(nInd*10*sizeof(double));
		}


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

		char* TAG;


		//TODO not filling with -INF? 
		// lngls.n=buf_size*10*nInd;
		// for(int i=0;i<lngls.n;i++){
		// lngls.data[i]=-INFINITY;
		// }





		/*
		 * [START] Reading sites
		 *
		 * hdr	header
		 * bcf	struct for storing each record
		 */

		while (bcf_read(in_ff, hdr, bcf) == 0) {

			get_data<float> lgl;
			//TODO check
			// lgl.data=(float*)malloc(buf_size*10*nInd*sizeof(float));

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
			//- if site is not shared between two inds
			//- only store lngls3 if majmin 3x3 will be used

			for(int indi=0; indi<nInd; indi++){
				for(int i=0;i<10;i++){
					lngls[nSites][(10*indi)+i]=NEG_INF;
					if(isnan(lgl.data[i])){
						fprintf(stderr,"\nMissing data\n");
						break;
					}else if (lgl.data[i]==bcf_float_vector_end){
						fprintf(stderr,"\n???\n");
					}else{
						// fprintf(stderr,"\n->->i:%d ind:%d %s %f -> %f\n",i,indi,hdr->samples[indi],lgl.data[(10*indi)+i],(double)log2ln(lgl.data[(10*indi)+i]));
						lngls[nSites][(10*indi)+i]=(double)log2ln(lgl.data[(10*indi)+i]);
					}
				}
			}

			get_data<int32_t> gt;

			gt.n = bcf_get_genotypes(hdr,bcf,&gt.data,&gt.size_e);
			gt.ploidy=gt.n/nInd;

			if(gt.n<0){
				fprintf(stderr,"\n[ERROR](File reading)\tVCF tag \"%s\" does not exist; will exit!\n\n",TAG);
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

		for(int i1=0;i1<nInd-1;i1++){
			for(int i2=i1+1;i2<nInd;i2++){

				int pair_idx=nCk_idx(nInd,i1,i2);

				// double SFS2D[10][10];
				// EM_2DSFS_GL10(lngls,SFS2D,i1,i2,nSites,args->tole);
				double SFS2D3[3][3];
				double dx;
				dx=EM_2DSFS_GL3(lngls,SFS2D3,i1,i2,nSites,args->tole,anc,der);
			// fprintf(stderr,"), ind1: (%f,%f,%f)\n",lngls[0][0],lngls[0][1],lngls[0][2]);
			// fprintf(stderr,"), ind2: (%f,%f,%f)\n",lngls[0][10],lngls[0][11],lngls[0][12]);
				// print_2DM(3,3,*SFS2D3);
				fprintf(stdout,"gl,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
						hdr->samples[i1],
						hdr->samples[i2],
						SFS2D3[0][0],SFS2D3[0][1],SFS2D3[0][2],
						SFS2D3[1][0],SFS2D3[1][1],SFS2D3[1][2],
						SFS2D3[2][0],SFS2D3[2][1],SFS2D3[2][2]);

				fprintf(stdout,"gle,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
						hdr->samples[i1],
						hdr->samples[i2],
						nSites*SFS2D3[0][0],nSites*SFS2D3[0][1],nSites*SFS2D3[0][2],
						nSites*SFS2D3[1][0],nSites*SFS2D3[1][1],nSites*SFS2D3[1][2],
						nSites*SFS2D3[2][0],nSites*SFS2D3[2][1],nSites*SFS2D3[2][2]);
				// fprintf(stderr,"\n\nHERE,%f\n\n",SFS2D3[0][0]);

				fprintf(stdout,"gt,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
						hdr->samples[i1],
						hdr->samples[i2],
						gt_sfs[pair_idx][0],gt_sfs[pair_idx][1],gt_sfs[pair_idx][2],
						gt_sfs[pair_idx][3],gt_sfs[pair_idx][4],gt_sfs[pair_idx][5],
						gt_sfs[pair_idx][6],gt_sfs[pair_idx][7],gt_sfs[pair_idx][8]);

			}
		}

		fprintf(stderr, "Total number of sites analysed: %i\n", nSites);

		bcf_hdr_destroy(hdr);
		bcf_destroy(bcf);

		int BCF_CLOSE;
		if ( (BCF_CLOSE=bcf_close(in_ff))){
			fprintf(stderr,"bcf_close(%s): non-zero status %d\n",in_fn,BCF_CLOSE);
			exit(BCF_CLOSE);
		}

		free(args->in_fn);
		free(args);


	}

	return 0;

}
