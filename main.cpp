/*
 *
 */


#include <stdio.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

#include <inttypes.h>
#include <math.h>

#include <time.h>

#include "main.h"
#include "random_generator.h"
#include "estimator.h"
#include "utils.h"
#include "io.h"
#include "shared.h"


using size_t=decltype(sizeof(int));

extern char bcf_allele_charToInt[256];

const int vcf_gl_order_idx[10]={0,1,4,2,5,7,3,6,8,9};


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

FILE *getFILE(const char*fname,const char* mode){
	FILE *fp;
	if(NULL==(fp=fopen(fname,mode))){
		fprintf(stderr,"[%s:%s()]\t->Error opening FILE handle for file:%s exiting\n",__FILE__,__FUNCTION__,fname);
		exit(1);
	}
	return fp;
}



int do_gt_sfs(bcf1_t *bcf,bcf_hdr_t *hdr,get_data<int32_t> gt ,int pairM[],char *TAG, int i1, int i2){

	gt.n = bcf_get_genotypes(hdr,bcf,&gt.data,&gt.size_e);

	if(gt.n<0){
		fprintf(stderr,"\n[ERROR](File reading)\tVCF tag \"%s\" does not exist; will exit!\n\n",TAG);
		exit(1);
	}



	// for(int i1=0;i1<nInd-1;i1++)
	// for(int i2=i1+1;i2<nInd;i2++)


	int32_t *ptr1 = gt.data + i1*gt.ploidy;
	int32_t *ptr2 = gt.data + i2*gt.ploidy;


	int gti1=0;
	int gti2=0;

	//binary input genotypes from simulated input
	for (int i=0; i<gt.ploidy;i++){
		gti1 += bcf_gt_allele(ptr1[i]);
		gti2 += bcf_gt_allele(ptr2[i]);
	}
	// fprintf(stderr,"\n%d:%d %d:%d\n",i1,gti1,i2,gti2);

	//0 1 2
	//00 01 02
	//MMMM MMMm MMmm
	//
	//3 4 5
	//10 11 12
	//MmMM MmMm Mmmm
	//
	//6 7 8
	//20 21 22
	//mmMM mmMm mmmm

	switch(gti1){
		case 0:
			switch(gti2){
				case 0:
					// pairM[i1][i2][0]++;
					pairM[0]++;
					break;
				case 1:
					// pairM[i1][i2][1]++;
					pairM[1]++;
					break;
				case 2:
					pairM[2]++;
					// pairM[i1][i2][2]++;
					break;
			}
			break;
		case 1:

			switch(gti2){
				case 0:
					// pairM[i1][i2][3]++;
					pairM[3]++;
					break;
				case 1:
					// pairM[i1][i2][4]++;
					pairM[4]++;
					break;
				case 2:
					// pairM[i1][i2][5]++;
					pairM[5]++;
					break;
			}
			break;
		case 2:

			switch(gti2){
				case 0:
					pairM[6]++;
					// pairM[i1][i2][6]++;
					break;
				case 1:
					pairM[7]++;
					// pairM[i1][i2][7]++;
					break;
			}
	}
	return 0;
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
	char *der;



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

		// int nInd=bcf_hdr_nsamples(hdr);
		nInd=bcf_hdr_nsamples(hdr);

		if(nInd==1){
			fprintf(stderr,"\n\n[ERROR]\tOnly one sample; will exit\n\n");
			exit(0);
		}
		//print sample ids
		// for( size_t i = 0; i < nInd; i++ ) {
		// fprintf(stderr,"\n->->%s\n",hdr->samples[i]);
		// }



		nSites=0;


		fprintf(stderr, "\nReading file:\t\"%s\"\n", in_fn);
		fprintf(stderr, "Number of samples: %i\n", bcf_hdr_nsamples(hdr));
		fprintf(stderr,	"Number of chromosomes: %d\n",hdr->n[BCF_DT_CTG]);


		//todo first use index of bcf etc to determine nsites
		//then start analyses
		//

		int buf_size=1024;
		// int buf_size=10;

		double *lngls[buf_size];

		for (int i=0;i<buf_size;i++){
			lngls[i]=(double*)malloc(nInd*10*sizeof(double));
		}



		int32_t *gts[buf_size];	
		for (int i=0;i<buf_size;i++){
			gts[i]=(int32_t*)malloc(nInd*10*sizeof(int32_t));
		}


		// lngls.data=(double*)malloc(buf_size*10*nInd*sizeof(double));

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


		// if (args->doGeno==1){
		int nDim=9;
		// int pairM[nInd][nInd*nDim] = (int*)malloc(nDim*nInd*sizeof(int));
		// int pairM[nInd][nInd][9] ={};
		//todo only create if dogeno 1
		int pairM[9]={};
		// }


		/*
		 * hdr	header
		 * bcf	struct for storing each record
		 */

		while (bcf_read(in_ff, hdr, bcf) == 0) {

			get_data<float> lgl;
			lgl.data=(float*)malloc(buf_size*10*nInd*sizeof(float));

			TAG="GL";
			lgl.n = bcf_get_format_float(hdr,bcf,TAG,&lgl.data,&lgl.size_e);

			if(lgl.n<0){
				fprintf(stderr,"\n[ERROR](File reading)\tVCF tag \"%s\" does not exist; will exit!\n\n",TAG);
				exit(1);
			}

			while(nSites>=buf_size){

				fprintf(stderr,"\n\nrealloc %d at site %d\n\n",buf_size,nSites);
				buf_size=buf_size*2;
				//TODO realloc lngls
				// lngls.data=(double*)realloc(lngls.data,buf_size*10*nInd*sizeof(double));

				// GL=(double**)realloc(GL,sizeof(double *[buf_size]));


				if(args->isSim==1){
					anc=(char*)realloc(anc,buf_size*sizeof(char));
					der=(char*)realloc(der,buf_size*sizeof(char));
				}

			}

			//TODO
			//-if site is not shared between two inds

			for(int indi=0; indi<nInd; indi++){
				for(int i=0;i<10;i++){
					if(isnan(lgl.data[i])){
						fprintf(stderr,"\nMissing data\n");
					}else if (lgl.data[i]==bcf_float_vector_end){
					}else{
						// lngls.data[(nSites*10*nInd)+(10*indi)+i]=(double)log2ln(lgl.data[i]);
						// GL[nSites][(10*indi)+i]=(double)log2ln(lgl.data[i]);
						lngls[nSites][(10*indi)+i]=(double)log2ln(lgl.data[i]);
						//TODO fix gt sfs 
						// gts[nSites][indi+i]=(int)
//
	// //binary input genotypes from simulated input
	// for (int i=0; i<gt.ploidy;i++){
		// gti1 += bcf_gt_allele(ptr1[i]);
		// gti2 += bcf_gt_allele(ptr2[i]);
	// }
					}
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

			if (args->doGeno==1){
				if (args->doInd==1){
					int i1=args->ind1;
					int i2=args->ind2;
					get_data<int32_t> gt;

					// gt.n = bcf_get_genotypes(hdr,bcf,&gt.data,&gt.size_e);

					// if(gt.n<0){
					// fprintf(stderr,"\n[ERROR](File reading)\tVCF tag \"%s\" does not exist; will exit!\n\n",TAG);
					// exit(1);
					// }
					do_gt_sfs(bcf,hdr,gt,pairM,"GT",i1,i2);
				}
				// fprintf(stderr,"\n\n\t-> Printing at site %d\n\n",nSites);
			}
			nSites++;
		}



		//TODO bcf_hdr_set_samples efficient sample parsing
		if (args->doGeno==1){

			if (args->doInd==1){
				int i1=args->ind1;
				int i2=args->ind2;
				fprintf(stdout,"ind%d,ind%d,",i1,i2);
				//TODO print sample ids,below. not sure if works
				// fprintf(stdout,"%s,%s,",hdr->samples[i1],hdr->samples[i2]);
				for (int k=0; k<nDim;k++){
					// fprintf(stderr,"%d",pairM[i1][i2][k]);
					fprintf(stdout,"%d",pairM[k]);
					if(k<nDim-1){
						fprintf(stdout,",");
					}else{
						fprintf(stdout,"\n");
					}
				}
			}

			// fprintf(stdout,"\ngt,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d",
					// hdr->samples[i1],
					// hdr->samples[i2],
					// pairM[0],
					// pairM[1],
					// pairM[2],
					// pairM[3],
					// pairM[4],
					// pairM[5],
					// pairM[6],
					// pairM[7],
					// pairM[8]);
		}


		// if(args->isSim==1){
		// // fprintf(stderr,"\n->here %d \n",anc[site]);
		// // fprintf(stderr,"\n->here %d \n",der[site]);
		// }
		//


		for(int i1=0;i1<nInd-1;i1++){
			for(int i2=i1+1;i2<nInd;i2++){

				// double SFS2D[10][10];
				// EM_2DSFS_GL10(lngls,SFS2D,i1,i2,nSites,args->tole);
				double SFS2D3[3][3];
				double dx;
				dx=EM_2DSFS_GL3(lngls,SFS2D3,i1,i2,nSites,args->tole,anc,der);
				// print_2DM(3,3,*SFS2D3);
				fprintf(stdout,"\ngl,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f",
						hdr->samples[i1],
						hdr->samples[i2],
						SFS2D3[0][0],SFS2D3[0][1],SFS2D3[0][2],
						SFS2D3[1][0],SFS2D3[1][1],SFS2D3[1][2],
						SFS2D3[2][0],SFS2D3[2][1],SFS2D3[2][2]);


			}
		}
		fprintf(stdout,"\n");

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
