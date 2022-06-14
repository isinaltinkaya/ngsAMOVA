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
#include "io.h"


double log2ln(float ivar){
	return (double) ivar/M_LOG10E;
}

const int vcf_gl_order_idx[10]={0,1,4,2,5,7,3,6,8,9};
//
// typedef struct VCF{
//
// VCF(char *ptr);
// ~VCF();
//
// void print();
//
// struct gt_data{
// uint8_t allele1;
// uint8_t allele2;
// bool is_phased;
// };
// gt_data gt;
//
//
// // double **gl;
//
//
// private:
// char *p;
// int len;
//
// }VCF;
//
//
// VCF::VCF(char *ptr){
// len=strlen(ptr);
// p=(char *)malloc(len+1);
// if(!p){
// fprintf(stderr, "\nAllocation error\n");
// exit(1);
// }
// strcpy(p,ptr);
// }
//
// VCF::~VCF(){
// fprintf(stderr, "\nDestructor is deallocating memory\n");
// free(p);
// }
//
//
// void VCF::print(){
// fprintf(stderr, "\nLength of %s is %d\n",p,len);
// }
//
//
FILE *getFILE(const char*fname,const char* mode){
	FILE *fp;
	if(NULL==(fp=fopen(fname,mode))){
		fprintf(stderr,"[%s:%s()]\t->Error opening FILE handle for file:%s exiting\n",__FILE__,__FUNCTION__,fname);
		exit(1);
	}
	return fp;
}

int main(int argc, char **argv) {


	if(argc==1){
		usage();
		exit(0);
	}
	argStruct *args=args_get(--argc,++argv);

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

		int nSamples=bcf_hdr_nsamples(hdr);


		fprintf(stderr, "\nReading file:\t\"%s\"\n", in_fn);
		fprintf(stderr, "Number of samples: %i\n", bcf_hdr_nsamples(hdr));
		fprintf(stderr,	"Number of chromosomes: %d\n",hdr->n[BCF_DT_CTG]);

		int nSites=0;


		//todo first use index of bcf etc to determine nsites
		//then start analyses
		//

		// int buf_size=1024;
		int buf_size=1;

		get_data<float> lgl;
		get_data<double> lngl;
		get_data<int32_t> gt;

		lngl.data=(double*)malloc(buf_size*10*nSamples*sizeof(double));
		lgl.data=(float*)malloc(buf_size*10*nSamples*sizeof(float));

		char* TAG;

		// lgl.data=new float[10*nSamples];

		
// fprintf(stderr,"\n\nHERE nsamples%d!!!\n\n",nSamples);
		lngl.n=buf_size*10*nSamples;
		for(int i=0;i<lngl.n;i++){
			lngl.data[i]=-INFINITY;
		}


		// if (args->doGeno==1){
			int nDim=9;
			// int pairM[nSamples][nSamples*nDim] = (int*)malloc(nDim*nSamples*sizeof(int));
			int pairM[nSamples][nSamples][9] ={}; 
		// }


		//TODO
		//todo read all sites to memory using n_sites_size as dynamic reallocated memory for these
		//keep track of nsites
		//then implement em
		while (bcf_read(in_ff, hdr, bcf) == 0) {

			TAG="GL";
			lgl.n = bcf_get_format_float(hdr,bcf,TAG,&lgl.data,&lgl.size_e);

			while(nSites>=buf_size){

				// fprintf(stderr,"\n\nrealloc %d at site %d\n\n",buf_size,nSites);
				buf_size=buf_size*2;
				lngl.data=(double*)realloc(lngl.data,buf_size*10*nSamples*sizeof(double));

			}
				// fprintf(stderr,"\n\nbufsize %d at site %d\n\n",buf_size,nSites);


					

			if(lgl.n<0){
				fprintf(stderr,"\n[ERROR](File reading)\tVCF tag \"%s\" does not exist; will exit!\n\n",TAG);
				exit(1);
			}

			// for(int ind=0;ind<nSamples;ind++){

// fprintf(stderr,"\n\nHERE!!!%d\n\n",lgl.n);
					fprintf(stderr,"\n\n");
			// for(int i=0; i<lgl.n; i++){
			for(int indi=0; indi<nSamples; indi++){
				for(int i=0;i<10;i++){
				if(isnan(lgl.data[i])){
					fprintf(stderr,"\nMissing data\n");
				}else if (lgl.data[i]==bcf_float_vector_end){
					// fprintf(stderr,"\n\nHERE!!!\n\n");
				}else{
					fprintf(stderr, "\n->lgl %lf",lgl.data[i]);
					lngl.data[(nSites*10*nSamples)+(10*indi)+i]=(double)log2ln(lgl.data[i]);
					fprintf(stderr, "->lngl %lf",lngl.data[i]);
				}
				}
			}
			// }


			// for(int ind=0;ind<10;ind++){
				// int32_t *ptr = gt.data + ind*gt.ploidy;
				// for(int i=0;i<gt.ploidy;i++){
// // fprintf(stderr,"\n\nHERE!!!\n\n");
					// // fprintf(stderr,"\n->\t%d",bcf_gt_allele(ptr[i]));
				// }
			// }

			if (args->doGeno==1){

			// TAG="GT";
			// gt.n = bcf_get_genotypes(hdr,bcf,&gt.data,&gt.size_e);
			// if(gt.n<0){
				// fprintf(stderr,"\n[ERROR](File reading)\tVCF tag \"%s\" does not exist; will exit!\n\n",TAG);
				// exit(1);
			// }
//
			// gt.ploidy=gt.n/nSamples;

				int32_t *gt_arr=NULL;
				int32_t ngt_arr=0;
				int ngt=bcf_get_genotypes(hdr, bcf, &gt_arr, &ngt_arr);
				if ( ngt<=0 ){
					fprintf(stderr,"\nGT not present\n");
				}

				int gt_ploidy=ngt/nSamples;

				if (gt_ploidy!=2){
					fprintf(stderr,"ERROR:\n\nploidy: %d not supported\n",gt_ploidy);
					return 1;
				}


				for(int i1=0;i1<nSamples-1;i1++){
					for(int i2=i1+1;i2<nSamples;i2++){


						int32_t *ptr1 = gt_arr + i1*gt_ploidy;
						int32_t *ptr2 = gt_arr + i2*gt_ploidy;


						int gti1=0;
						int gti2=0;

						//binary input genotypes from simulated input
						for (int i=0; i<gt_ploidy;i++){
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
										pairM[i1][i2][0]++;
										break;
									case 1:
										pairM[i1][i2][1]++;
										break;
									case 2:
										pairM[i1][i2][2]++;
										break;
								}
								break;
							case 1:

								switch(gti2){
									case 0:
										pairM[i1][i2][3]++;
										break;
									case 1:
										pairM[i1][i2][4]++;
										break;
									case 2:
										pairM[i1][i2][5]++;
										break;
								}
								break;
							case 2:

								switch(gti2){
									case 0:
										pairM[i1][i2][6]++;
										break;
									case 1:
										pairM[i1][i2][7]++;
										break;
									case 2:
										pairM[i1][i2][8]++;
										break;
								}
								break;
						}
					}
				}
			}



		nSites++;
		}


		if (args->doGeno==1){

			for(int i1=0;i1<nSamples-1;i1++){
				for(int i2=i1+1;i2<nSamples;i2++){
					fprintf(stderr,"ind%d,ind%d,",i1,i2);
					for (int k=0; k<nDim;k++){
						fprintf(stderr,"%d",pairM[i1][i2][k]);
						if(k<nDim-1){
							fprintf(stderr,",");
						}else{
							fprintf(stderr,"\n");
						}
					}

				}

			}
		}

		for(int ind=0; ind<nSamples; ind++){
			for (int site=0; site<nSites;site++){
				for (int i=0;i<10;i++){
				fprintf(stderr, "\nafter lngl->\t%lf\n",lngl.data[(site*10*nSamples)+(10*ind)+i]);
				}
			}
		}

		// free(gt_arr);
		//
		// for(int i1=0;i1<nSamples-1;i1++){
		// for(int i2=i1+1;i2<nSamples;i2++){
		// fprintf(stderr,"ind%d,ind%d,",i1,i2);
		// for (int k=0; k<nDim;k++){
		// fprintf(stderr,"%d",pairM[i1][i2][k]);
		// if(k<nDim-1){
		// fprintf(stderr,",");
		// }else{
		// fprintf(stderr,"\n");
		// }
		// }
		// }
		// }

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
