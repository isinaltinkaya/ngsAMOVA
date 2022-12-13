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

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))



/*
 * Macro:[DBL_MAXDIG10]
 * Defines the maximum number of decimal digits that can be represented by a
 * double-precision floating-point number on a 64-bit system
 *
 * Calculates the max number of decimal digits that a double-precision
 * floating-point number can have on a 64-bit system using the value of DBL_MANT_DIG
 * The UL suffix ensures that our constant values are unsigned long 
 * to prevent any overflow issues.
 *
 */
#define DBL_MAXDIG10 (2+ (DBL_MANT_DIG * 30103UL)/100000UL)


/*
 * Macro:[AT]
 * Injects the file and line info as string
 */

#define STRINGIFY(x) #x
#define ASSTR(x) STRINGIFY(x)
#define AT __FILE__ ":" ASSTR(__LINE__)

/*
 * Macro:[ASSERT]
 * shortcut to evaluate an expression, works the same way as the C-macro assert
 */
#define ASSERT(expr) if (!(expr)) {fprintf(stderr,"\n\n*******\n[ERROR](%s:%d) %s\n*******\n",__FILE__,__LINE__,#expr);exit(1);}


extern const int get_3x3_idx[3][3];

const double NEG_INF = -std::numeric_limits<double>::infinity();

using size_t=decltype(sizeof(int));


namespace DATA{


	typedef struct contigsStruct{
		int nContigs=0;
		char **contigNames;
		size_t _contigNames=24;
		int *contigLengths;
		size_t _contigLengths=24;

		int **contigBlockStarts = NULL;

		//realloc if needed	
		//expand
		// void realloc_contigBlockStarts(int **contigBlockStarts, int ci, int nBlocks){
		// 	contigBlockStarts[ci] = (int *)realloc(contigBlockStarts[ci], nBlocks * sizeof(int));
		// }

		contigsStruct(){
			contigNames=new char*[_contigNames];
			contigLengths=new int[_contigLengths];
			contigBlockStarts = (int **)malloc(nContigs * sizeof(int *));
			for (int ci = 0; ci < nContigs; ci++)
			{
				contigBlockStarts[ci] = NULL;
			}
		}
		~contigsStruct(){
			for (size_t i=0;i<_contigNames;i++){
				free(contigNames[i]);
				contigNames[i]=NULL;
			}
			delete [] contigNames;
			delete [] contigLengths;
			for(int i = 0; i < nContigs; i++)
			{
				free(contigBlockStarts[i]);
				contigBlockStarts[i] = NULL;
			}
			free(contigBlockStarts);
			contigBlockStarts = NULL;

		}


	} contigsStruct;

	typedef struct pairStruct{

		int i1;
		int i2;
		int idx;

		size_t snSites=0;

		int *sSites;
		// TODO init better
		size_t _sSites=100;


		// int nDim;

		double d;
		int n_em_iter;

		//TODO should not always create below
		// double SFS[3][3]={{NEG_INF,NEG_INF,NEG_INF},
			// {NEG_INF,NEG_INF,NEG_INF},
			// {NEG_INF,NEG_INF,NEG_INF}};

		double *SFS = NULL;


		pairStruct(int ind1, int ind2, int pair_idx){
			i1=ind1;
			i2=ind2;
			idx=pair_idx;
			d=0.0;
			n_em_iter=0;
			sSites=(int*) malloc(_sSites*sizeof(int));
			for (size_t i=0;i<_sSites;i++){
				sSites[i]=-1;
			}

			SFS = (double *) malloc (9*sizeof(double));
			for (int i = 0; i < 9; i++) {
				SFS[i] = NEG_INF;
			}
		}
		~pairStruct(){
			free(sSites);
			sSites=NULL;

			free(SFS);
			SFS=NULL;
		}

	}pairStruct;


	typedef struct samplesStruct{

		int *sampleArr;
		size_t _sampleArr=200;
		//TODO increase if overflow

		samplesStruct(){
			// sampleArr=new int[_sampleArr];
			sampleArr=(int*) malloc(_sampleArr*sizeof(int));
			for (size_t i=0;i<_sampleArr;i++){
				sampleArr[i]=0;
			}
		}
		~samplesStruct(){
			free(sampleArr);
			sampleArr=NULL;
			// for (size_t i=0;i<_sampleArr;i++){
			// free(sampleArr[i]);
			// sampleArr[i]=NULL;
			// }
			// delete [] sampleArr;
		}

	}samplesStruct;

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
			id=NULL;
		}

	}strataStruct;

	typedef struct metadataStruct{

		strataStruct *strataArr;
		size_t _strataArr=10;

		int nInds_total=0;

		int nStrata=0;

		metadataStruct(){
			strataArr=new strataStruct[_strataArr];
			fprintf(stderr,"->Creating metadataStruct");
		};
		~metadataStruct(){
			delete [] strataArr;
			strataArr=NULL;
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
 * 							anc=ref and der=alt[0]
 *
 * @field minInd		[-1 = not set]
 * 						minimum number of individuals needed
 * 						for site to be included in analyses
 * 						
 * 						if minInd not set; set minInd=2
 * 						== no filter, include all sites that exist in pair
 *
 * 						if minInd==nInd; set minInd=0
 * 						== site should be nonmissing for all individuals
 * 						to be included
 *
 * @field blockSize		[0 = not set]
 *						Block size to be used in block bootstrapping
 *
 *
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
 * 						[2] use Fij F statistic [DEPRECATED]
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
 * @field mEmIter		maximum number of iterations allowed for em
 *
 */


typedef struct {


	char* in_fn;
	char* in_sfs_fn;
	char* in_mtd_fn;
	char* out_fp;

	int whichCol;
	int blockSize;
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
	int mEmIter;


	int gl2gt;


}argStruct;


argStruct *argStruct_init();

argStruct *argStruct_get(int argc, char **argv);







namespace IO {


	char *setFileName(const char* a,const char* b);


	FILE *getFILE(const char*fname,const char* mode);
	FILE *openFILE(const char* a,const char* b);
	FILE* openFILE(char* c);

	namespace readFILE{
		int METADATA(DATA::metadataStruct * MTD, FILE* in_mtd_ff, int whichCol, const char* delims, DATA::samplesStruct *SAMPLES);
		int SFS(FILE* in_sfs_ff, const char* delims, DATA::samplesStruct *SAMPLES);
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
				out_dm_fs= new outputStruct(args->out_fp,".distance_matrix.csv");
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

	namespace print{


		/* IO::print::Array
		 * Prints the elements of an array to a file or stream.
		 *
		 * @param arr: the array to print
		 * @param N: the number of rows, dimension 1 
		 * @param M: the number of columns, dimension 2
		 * @param out: the file or stream to print to
		 * @param sep: the character to use as a separator
		 *
		 * @example
		 * IO::print::Array(stdout, myArray, DOUBLE, 3, 3, ',');
		 */
		void Array(FILE *fp, double *arr, size_t N, size_t M, char sep);

		// IO::print::Array
		// :overload: int array
		void Array(FILE *fp, int *arr, size_t N, size_t M, char sep);
	}

}


typedef struct threadStruct{

	DATA::pairStruct* pair;
	double **lngls;



	FILE* out_sfs_ff;

	// double* M_PWD_GL_PAIR;

	//TODO use them globally?
	double tole;
	int mEmIter;


	size_t nSites;

	threadStruct(DATA::pairStruct* tPair, double **lngl, size_t nSites_t, IO::outputStruct* out_sfs_fs, double toleArg, int mEmIterArg){
		pair=tPair;
		lngls=lngl;
		nSites=nSites_t;
		out_sfs_ff=out_sfs_fs->ff;
		// M_PWD_GL_PAIR=M_PWD_GL_P;
		tole=toleArg;
		// doDist=doDistArg;
		mEmIter=mEmIterArg;
	}

}threadStruct;


/*
 * @typedef
 * @abstract paramStruct - parameter structure
 *
 * @field nSites				number of sites
 * @field nInd					number of individuals
 * @field pos					position
 *
 * @field LUT_indPair_idx		lookup table for mapping two individuals
 * 								to their pair index
 * @field n_ind_cmb				number of unique pairwise individual combinations
 *
 * @field major					major allele
 * @field minor					minor allele
 * @field ref					reference allele
 * @field anc					ancestral allele
 * @field der					derived allele
 */


typedef struct{

	size_t nSites;
	size_t totSites;

	int nInd;

	int *pos;

	int **LUT_indPair_idx;
	int n_ind_cmb;


	char *major;
	char *minor;
	char *ref;
	char *anc;
	char *der;

	char *DATETIME;


}paramStruct;



paramStruct *paramStruct_init(argStruct *args);
void paramStruct_destroy(paramStruct *p);


char *get_time();
// void *argStruct_destroy(argStruct *arg);

void usage(FILE *fp);

#endif
