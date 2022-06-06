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
		exit(0);
	}
	return fp;
}


//
// FILE *openFile(const char* a,const char* b){
// char *c = (char*)malloc(strlen(a)+strlen(b)+1);
// strcpy(c,a);
// strcat(c,b);
// // fprintf(stderr,"\t-> Dumping file: %s\n",c);
// FILE *fp = getFILE(c,"w");
// free(c);
// return fp;
// }
//

int main(int argc, char **argv) {


	if(argc==1){
		usage();
		return 0;
	}
	argStruct *args=args_get(--argc,++argv);

	if(args!=NULL){

		char *in_fn=args->in_fn;

		vcfFile * in_ff = bcf_open(in_fn, "r");



		if (in_ff == NULL) {
			return 1;
		}

		if (bcf == 0) {
			return 1; 
		}

		bcf_hdr_t *hdr = bcf_hdr_read(in_ff);



		bcf1_t *bcf = bcf_init();


		int nSamples=bcf_hdr_nsamples(hdr);


		fprintf(stderr, "\nReading file:\t\"%s\"\n", in_fn);
		fprintf(stderr, "Number of samples: %i\n", bcf_hdr_nsamples(hdr));
		fprintf(stderr,	"Number of chromosomes: %d\n",hdr->n[BCF_DT_CTG]);



		int nSites=0;

		while (bcf_read(in_ff, hdr, bcf) == 0) {


			get_data<float> lgl;
			get_data<double> lngl;

			// lngl.data=new double[10*nSamples];
			lngl.data=(double*)malloc(10*nSamples*sizeof(double));

			// lgl.data=new float[10*nSamples];

			for(int i=0;i<10*nSamples;i++){
				lngl.data[i]=-INFINITY;
			}

			const char* TAG;
			TAG="GL";
			lgl.n = bcf_get_format_float(hdr,bcf,TAG,&lgl.data,&lgl.size_e);


			if(lgl.n<0){
				fprintf(stderr,"\n[ERROR](File reading)\tVCF tag \"%s\" does not exist; will exit!\n\n",TAG);
				exit(1);
			}

			fprintf(stderr,"\n");
			for(int i=0; i<lgl.n; i++){
				if(isnan(lgl.data[i])){
				//missing does not work??
				// if(lgl.data[i]==bcf_float_missing){
					fprintf(stderr,"\n\nHERE!!!\n\n");
				}else if (lgl.data[i]==bcf_float_vector_end){
					fprintf(stderr,"\n\nHERE!!!\n\n");
				}else{
					// fprintf(stderr, "\n->\t%lf\n",lgl.data[i]);
					// fprintf(stderr, "\n->\t%lf\n",lngl.data[i]);
					lngl.data[i]=(double)log2ln(lgl.data[i]);
					fprintf(stderr, "%lf ",lngl.data[i]);

					fprintf(stderr,"\n");
				}
			}

			for(int ind=0;ind<10;ind++){
				for(int j=0;j<10;j++){
				}
			}

			// for(int i1=0;i1<nSamples-1;i1++){
			// for(int i2=i1+1;i2<nSamples;i2++){
			// }
			// } //end sample loop


			nSites++;
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
