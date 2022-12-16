#include "shared.h"
#include "math_utils.h"

// int doAMOVA(int nInd, int n_ind_cmb, DATA::metadataStruct *MTD, DATA::sampleStruct* SAMPLES, FILE *out_amova_ff, int sqDist, double *M_PWD, int **LUT_indPair_idx);

int doAMOVA(int n_ind_cmb, int nInd, DATA::metadataStruct *MTD, DATA::samplesStruct *SAMPLES, FILE *out_amova_ff, int sqDist, double *M_PWD, int **LUT_indPair_idx, const char *type);
