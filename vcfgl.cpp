/*
 *
 * VCF parser
 *
 *
 *
 * [[ MSPRIME SIMULATED VCF INPUT ]]
 * generated by tskit.write_vcf()
 *
 *
 * ** [POSITIONS]
 * Positions are 0-based
 * in contrast to regular vcf specs=1-based
 *
 * ** [ALLELES]
 * * REF(=ANC)=0
 * Reference allele is the same as ancestral allele
 * Encoded as 0
 * * ALT=1
 * Only one ALT allele is present
 * Encoded as 1
 * 
 *
 * [GL order in VCF format]
 *
 * for P=ploidy and N=number of alternate alleles;
 * for a_p in 0:N; for a_p-1 in 0:a_p; print(a1,a2);
 *
 * For P=2 N=3
 *
 * 0,1,2,3
 * A,C,G,T
 * 00,01,11,02,12,22,03,13,23,33
 * AA,AC,CC,AG,CG,GG,AT,CT,GT,TT
 *
 *
 * ** GLF
 * A, C, G, T
 *
 */

#include <stdio.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

#include <inttypes.h>
#include <math.h>

#include <time.h>

#include "random_generator.h"
#include "estimator.h"
#include "io.h"



const int vcf_gl_order_idx[10]={0,1,4,2,5,7,3,6,8,9};



/*
 * @template struct
 * @abstract	wrapper for bcf_get_data_*
 *
 * @field data
 * @field size_e	watermark for number of elements
 * @field n			number of returned values
 * 					if <0; error
 * 					else; number of written values
 * 					used for accessing entries
 *
 */
template<typename T> struct get_data{

	T *data = NULL;

	int size_e=0;
	int n=0;
	

	T& operator[](unsigned i){
		return data[i];
	}

	bool is_empty() const {
		return data == NULL;
	}

	void destroy(){
		free(data);
		data=NULL;
	}

	~get_data(){
		free(data);
	}
};


typedef struct VCF{

	VCF(char *ptr);
	~VCF();

	void print();

	struct gt_data{
		uint8_t allele1;
		uint8_t allele2;
		bool is_phased;
	};
	gt_data gt;

	
	// double **gl;


	private:
	char *p;
	int len;

}VCF;


VCF::VCF(char *ptr){
	len=strlen(ptr);
	p=(char *)malloc(len+1);
	if(!p){
		fprintf(stderr, "\nAllocation error\n");
		exit(1);
	}
	strcpy(p,ptr);
}

VCF::~VCF(){
	fprintf(stderr, "\nDestructor is deallocating memory\n");
	free(p);
}


void VCF::print(){
	fprintf(stderr, "\nLength of %s is %d\n",p,len);
}


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
		fprintf(stderr,"\n\nhelp\n\n");
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


		// fprintf(stderr, "\nReading file:\t\"%s\"\n", in_fn);
		// fprintf(stderr, "Number of samples: %i\n", bcf_hdr_nsamples(hdr));
		// fprintf(stderr,	"Number of chromosomes: %d\n",hdr->n[BCF_DT_CTG]);
		//
		// fprintf(stderr, "\n\n\n");
		//
		//

		double **gl =NULL;
		gl=new double *[10*nSamples];



		int nSites=0;
		while (bcf_read(in_ff, hdr, bcf) == 0) {


			get_data<float> lgl;

			lgl.n = bcf_get_format_float(hdr,bcf,"GL",&lgl.data,&lgl.size_e);

			
			if(lgl.n<0){
				fprintf(stderr,"\n\nn %d\n\n",lgl.n);
			}
			
			// double ***pgl=NULL;
			
			// get_data<double> lngl;
			// lngl.data=(double*)malloc(10*nSamples*sizeof(double));

			// for(int i=0; i<lgl.n; i++){
				// fprintf(stderr,"%f ",lgl.data[i]);
				// lngl.data[i]=(double)lgl.data[i];
				// fprintf(stderr,"%lf\n",lngl.data[i]/M_LOG10E);
			// }
			for(int ind=0;ind<10;ind++){
				for(int j=0;j<10;j++){
					// gl[nSites][ind*10+j]=(double)lgl.data[ind];
					// *gl[ind*10+j]=(double)lgl.data[ind];
					// gl[ind]=(double)lgl.data[ind];
					// gl[ind]=(double*)&lgl.data[ind];
					*gl[ind]=(double)lgl.data[ind];
					fprintf(stderr,"%f ",lgl.data[ind]);
					// fprintf(stderr,"%lf\n",*gl[ind]);
					// fprintf(stderr,"%f", gl[ind*10+j]);
					// gl[nSites][ind*10+j]=(double)lgl.data[ind]/M_LOG10E;
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

		// fprintf(stderr, "Total number of sites: %i\n", nSites);

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
