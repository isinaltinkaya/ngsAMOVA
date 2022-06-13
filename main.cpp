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

		int n_sites_size=100;
fprintf(stderr,"\nXXX!!%u\n\n",bcf->n_allele);

		get_data<float> lgl;
		get_data<double> lngl;
		get_data<int32_t> gt;

		// lngl.data=new double[10*nSamples];
		// lngl.data=(double*)malloc(10*nSamples*sizeof(double));
		// lgl.data=(float*)malloc(10*nSamples*sizeof(float));
		// gt.data=(double*)malloc(10*nSamples*sizeof(double));


		lngl.data=(double*)malloc(n_sites_size*10*nSamples*sizeof(double));
		lgl.data=(float*)malloc(n_sites_size*10*nSamples*sizeof(float));

		char* TAG;

		// lgl.data=new float[10*nSamples];

		
fprintf(stderr,"\n\nHERE nsamples%d!!!\n\n",nSamples);
		lngl.n=10*nSamples;
		for(int i=0;i<lngl.n;i++){
			lngl.data[i]=-INFINITY;
		}


		//TODO
		//todo read all sites to memory using n_sites_size as dynamic reallocated memory for these
		//keep track of nsites
		//then implement em
		while (bcf_read(in_ff, hdr, bcf) == 0) {

			TAG="GL";
			lgl.n = bcf_get_format_float(hdr,bcf,TAG,&lgl.data,&lgl.size_e);



			if(lgl.n<0){
				fprintf(stderr,"\n[ERROR](File reading)\tVCF tag \"%s\" does not exist; will exit!\n\n",TAG);
				exit(1);
			}

		fprintf(stderr,"\n->loggl%d\n",lgl.n);
		fprintf(stderr,"\n->lngl%d\n",lngl.n);
			fprintf(stderr,"\n");
			for(int i=0; i<lgl.n; i++){
				if(isnan(lgl.data[i])){
					// fprintf(stderr,"\nMissing data\n");
				}else if (lgl.data[i]==bcf_float_vector_end){
					// fprintf(stderr,"\n\nHERE!!!\n\n");
				}else{
					// fprintf(stderr, "\n->\t%lf\n",lgl.data[i]);
					// fprintf(stderr, "\n->\t%lf\n",lngl.data[i]);
					lngl.data[i]=(double)log2ln(lgl.data[i]);
					// fprintf(stderr, "%lf ",lngl.data[i]);
					// fprintf(stderr,"\n");
				}
			}



			TAG="GT";
			gt.n = bcf_get_genotypes(hdr,bcf,&gt.data,&gt.size_e);
			if(gt.n<0){
				fprintf(stderr,"\n[ERROR](File reading)\tVCF tag \"%s\" does not exist; will exit!\n\n",TAG);
				exit(1);
			}

			gt.ploidy=gt.n/nSamples;
			for(int ind=0;ind<10;ind++){
				int32_t *ptr = gt.data + ind*gt.ploidy;
				for(int i=0;i<gt.ploidy;i++){
// fprintf(stderr,"\n\nHERE!!!\n\n");
					// fprintf(stderr,"\n->\t%d",bcf_gt_allele(ptr[i]));
				}
			}

			for(int i1=0;i1<nSamples-1;i1++){
			for(int i2=i1+1;i2<nSamples;i2++){
			}
			} //end sample loop

		fprintf(stderr,"\n->loggl%d\n",lgl.n);
		fprintf(stderr,"\n->lngl%d\n",lngl.n);

			nSites++;
		}

		fprintf(stderr,"\n->loggl%d\n",lgl.n);
		fprintf(stderr,"\n->lngl%d\n",lngl.n);
		for(int i=0; i<lngl.n; i++){
			// if(isnan(lgl.data[i])){
				// fprintf(stderr,"\nMissing data\n");
			// }else if (lgl.data[i]==bcf_float_vector_end){
				// fprintf(stderr,"\n\nHERE!!!\n\n");
			// }else{
				// fprintf(stderr, "\n->\t%lf\n",lgl.data[i]);
				// fprintf(stderr, "\n->\t%lf\n",lngl.data[i]);
				fprintf(stderr, "%lf ",lngl.data[i]);
				fprintf(stderr,"\n");
			// }
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
