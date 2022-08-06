#ifndef __PARAM_STRUCT__
#define __PARAM_STRUCT__


#include <limits>
#include <cstddef>
#include <time.h>

#include <cstdlib>
#include <cstdio>

#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/stat.h>

#include <math.h>

const double NEG_INF = -std::numeric_limits<double>::infinity();

using size_t=decltype(sizeof(int));


namespace DATA{

	typedef struct pairStruct{

		int i1;
		int i2;
		int idx;

		int shared_nSites=0;

		// int nDim;

		double d;
		int n_em_iter;

		// double SFS[3][3]={NEG_INF};
		double SFS[3][3]={{NEG_INF,NEG_INF,NEG_INF},
			{NEG_INF,NEG_INF,NEG_INF},
			{NEG_INF,NEG_INF,NEG_INF}};
		// double **SFS;
		// double SFS[nDim][nDim];

		pairStruct(int ind1, int ind2, int pair_idx){
			i1=ind1;
			i2=ind2;
			idx=pair_idx;
			d=0.0;
			n_em_iter=0;
		}

		
	}pairStruct;


	typedef struct sampleStruct{

		int *sampleArr;
		int _size_sampleArr=200;

		sampleStruct(){
			sampleArr=new int[_size_sampleArr];
			for(int i=0;i<_size_sampleArr;i++){
				sampleArr[i]=0;
			}
		}
		~sampleStruct(){
			delete [] sampleArr;
		}

	}sampleStruct;

	typedef struct strataStruct{

		int nInds=0;
		char *id;

	//TODO Associate hierarchical levels
	//
		int assoc=0;

		strataStruct(){
			id=NULL;
		}
		~strataStruct(){
			free(id);
		}

	}strataStruct;

	typedef struct metadataStruct{

		strataStruct *strataArr;
		int _size_strataArr=3;

		int nInds_total=0;

		int nStrata=0;

		metadataStruct(){
			strataArr=new strataStruct[_size_strataArr];
		};
		~metadataStruct(){
			delete [] strataArr;
		}

	}metadataStruct;



};
/*
 * @typedef
 * @abstract argStruct - argument structure
 *
 * @field *in_fn		pointer to input file name
 * @field *in_mtd_fn	pointer to input metadata file name
 * @field *out_fp		pointer to output file prefix [angsdput]
 * @field seed			random seed
 *
 * @field isSim			input is vcfgl simulation output
 * 						anc=ref and der=alt[0]
 *
 *
 * @field minInd		[-1 = not set]
 * 						minimum number of individuals needed
 * 						for site to be included in analyses
 *
 * @field whichCol		[-1] defines the index (1-based) of the
 * 						column in metadata file 'in_mtd_fn'
 * 						to use to define stratification levels
 *
 * @field doAMOVA		[0]
 * 						1 use 10 genotype likelihoods (GL)
 * 						2 use genotypes (GT) (NOTE: Only for benchmark purposes)
 * 
 * @field doDist		[0] use Sij similarity index
 * 						[1] use Dij (1-Sij) dissimilarity index
 * 						[2] use Fij F statistic
 *
 * @field sqDist		[0] use absolute value of distance measure (|dist_ij|)
 * 						[1] use squared distance measure (dist_ij^2)
 *
 * @field doInd			do ind pairs
 * @field ind1			ind1 id
 * @field ind2			ind2 id
 *
 *
 * @field doTest		test for convergence
 *
 *
 * @field mThreads		maximum number of threads defined by user
 *
 */


typedef struct {


	char* in_fn;
	char* in_mtd_fn;
	char* out_fp;

	int whichCol;
	int doAMOVA;
	int printMatrix;
	int isSim;
	int doDist;
	int sqDist;
	int minInd;

	int seed;

	double tole;

	int doInd;
	int ind1;
	int ind2;

	int doTest;


	int mThreads;
	

}argStruct;


argStruct *argStruct_init();

argStruct *argStruct_get(int argc, char **argv);







namespace IO {


	char *setFileName(const char* a,const char* b);

	
	FILE *getFILE(const char*fname,const char* mode);
	FILE *openFILE(const char* a,const char* b);
	FILE* openFILE(char* c);

	namespace readFILE{
		int METADATA(DATA::metadataStruct * MTD, FILE* in_mtd_ff, int whichCol, const char* delims, DATA::sampleStruct *SAMPLES);
	};

	namespace inspectFILE{
		int count_nColumns(char* line, const char* delims);
	};


	typedef struct outputStruct{
		char* fn=NULL;
		FILE* ff=NULL;

		outputStruct(const char* fp, const char* suffix){
			fn=setFileName(fp, suffix);
			ff=openFILE(fn);
		}
		~outputStruct(){
			fclose(ff);
			free(fn);
			fn=NULL;
		}

	}outputStruct;

	typedef struct outFilesStruct{
		outputStruct* out_emtest_fs=NULL;
		outputStruct* out_sfs_fs=NULL;
		outputStruct* out_dm_fs=NULL;
		outputStruct* out_amova_fs=NULL;

		outFilesStruct(argStruct* args){
			if(args->printMatrix==1){
				out_dm_fs= new outputStruct(args->out_fp,".dm.csv");
			}
			if(args->doTest==1){
				out_emtest_fs= new outputStruct(args->out_fp,".emtest.csv");
			}

			out_amova_fs= new outputStruct(args->out_fp,".amova.csv");
			out_sfs_fs= new outputStruct(args->out_fp,".sfs.csv");
		}

		~outFilesStruct(){
			delete out_emtest_fs;
			delete out_dm_fs;
			delete out_sfs_fs;
			delete out_amova_fs;
		}

	}outFilesStruct;

	typedef struct threadStruct{

		DATA::pairStruct* pair;
		double **lngls;

		size_t total_nSites;


		FILE* out_sfs_ff;

		double* M_PWD_GL_PAIR;

		int doDist;
		double tole;

		threadStruct(DATA::pairStruct* tPair, double **lngl, size_t nSites, outputStruct* out_sfs_fs, double* M_PWD_GL_P, double toleArg, int doDistArg){
			pair=tPair;
			lngls=lngl;
			total_nSites=nSites;
			out_sfs_ff=out_sfs_fs->ff;
			M_PWD_GL_PAIR=M_PWD_GL_P;
			tole=toleArg;
			doDist=doDistArg;
			// if(args->minInd!=0){
				// pair->shared_nSites=nSites;
			// }
		}

	}threadStruct;
}



/*
 * @typedef
 * @abstract paramStruct - parameter structure
 *
 * @field nSites	number of sites
 * @field nInd		number of individuals
 * @field pos		position
 * @field keepSites	int array with values indicating if a site will be kept
 * @field major		major allele
 * @field minor		minor allele
 * @field ref		reference allele
 * @field anc		ancestral allele
 * @field der		derived allele
 */


typedef struct{

	size_t nSites;
	int nInd;

	int *pos;

	int *keepSites;

	char *major;
	char *minor;
	char *ref;
	char *anc;
	char *der;

	char *DATETIME;

	
}paramStruct;



paramStruct *paramStruct_init(argStruct *args);
void paramStruct_destroy(paramStruct *p);

extern const int get_3x3_idx[3][3];

char *get_time();
// void *argStruct_destroy(argStruct *arg);

void usage();

#endif
