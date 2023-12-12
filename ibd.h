#ifndef __IBD__
#define __IBD__

#include "dataStructs.h"
#include "io.h"
// #include "mathUtils.h"
#include "shared.h"
#include "paramStruct.h"
#include "vcfReader.h"

struct paramStruct;

// void readSites_ibdseq(vcfData *vcfd, paramStruct *pars, ibdStruct* ibd);
void readSites_doIbd(vcfData *vcfd, paramStruct *pars);

int get_ibdData(const int contig_i, const int site_i, vcfData *vcfd, paramStruct *pars);

typedef struct ibdStruct {

	int nSites=0;
	int nInd=0;
	// void add_doseData(const int site_i,const int ind_i, const int dose);

	// doseData[nSites][nInd]
	// int** doseData=NULL;

	double* maxErrorArray=NULL;
	double* ibdScores = NULL;
	double* hbdScores = NULL;

	double **pairScores=NULL;
	// double **selfScores=NULL;


	// ----------
	// // nIndices=6 from ibdseq
	// pars->ibdScores = (double**) malloc(6 * sizeof(double*));
	// pars->ibdScores[0] = (double*) malloc(sizeof(double));
	// pars->ibdScores[1] = (double*) malloc(sizeof(double));
	// pars->ibdScores[2] = (double*) malloc(sizeof(double));
	// pars->ibdScores[3] = (double*) malloc(sizeof(double));
	// pars->ibdScores[4] = (double*) malloc(sizeof(double));
	// pars->ibdScores[5] = (double*) malloc(sizeof(double));
	// // ----------
//


	double ibdLike(int dose1, int dose2, double err, double fB);
	double ibdLike_GL(const int idx, vcfData* vcfd, paramStruct* pars,double fm);
	double nullLike(int dose1, int dose2, double fB);
	double nullLike_GL(const int idx,vcfData* vcfd, paramStruct* pars,double fm);
	double errorRate(double fB);
	double* errorArray(double e);
	double ibdScore(const int dose1, const int dose2, const double fB);
	double ibdScore_GL(const int i1, const int i2,vcfData* vcfd, paramStruct* pars, const double fa, const size_t site_i,double* lngl);
	double hbdScore(int dose, double fB);

	ibdStruct(vcfData* vcfd, paramStruct* pars);
	// ibdStruct(vcfData* vcfd);
	~ibdStruct();

} ibdStruct;

ibdStruct *ibdStruct_get(paramStruct *pars, vcfData *vcfd);

// int get_ibdScores(const int contig_i, const int site_i, vcfData *vcfd, paramStruct *pars);
int site_doIbd(const int contig_i, const int site_i, vcfData *vcfd, paramStruct *pars,ibdStruct* ibd);
// int site_doIbd_GL(const int contig_i, const int site_i, vcfData *vcfd, paramStruct *pars,ibdStruct* ibd,double* lngl);
double score(paramStruct* pars, vcfData* vcfd, const int sampleA, const int sampleB, const int marker,const  bool isCor);

int trimmedEnd(paramStruct* pars, const int i1, const int i2, const int start, const int end);
int trimmedStart(paramStruct* pars, const int i1, const int i2, const int start, const int end);
int get_ibd_segments(const char* contigName, paramStruct *pars, int* site2pos);

#endif // __IBD__
