// developer functions

#include "vcfUtils.h"
#include "shared.h"
#include "em.h"


// // print the binary representation of a char
// #define PRINT_CHAR_BITS(a) 
// 	printf("%d%d%d%d%d%d%d%d", !!((a << 0) & 0x80), !!((a << 1) & 0x80), !!((a << 2) & 0x80), !!((a << 3) & 0x80), !!((a << 4) & 0x80), !!((a << 5) & 0x80), !!((a << 6) & 0x80), !!((a << 7) & 0x80));

// print the binary representation of a char with additional info
//      with info about where this function is called from,name of the variable
//       and the value of the variable
#define PRINT_CHAR_BITS(a) \
    printf("\n\n[DEV][PRINT_CHAR_BITS]\t%s:%d: var:%s, val:%d, 0b:%d%d%d%d%d%d%d%d\n\n", __FILE__, __LINE__, #a, a, !!((a << 0) & 0x80), !!((a << 1) & 0x80), !!((a << 2) & 0x80), !!((a << 3) & 0x80), !!((a << 4) & 0x80), !!((a << 5) & 0x80), !!((a << 6) & 0x80), !!((a << 7) & 0x80));

// print the binary representation of a char with additional info
//      [DEV][PRINT_CHAR_BITS]    dev.cpp:123: a = 00000000
    // fprintf(stderr, "\n\n[DEV][PRINT_CHAR_BITS]\t%s:%d: %s = %d%d%d%d%d%d%d%d\n", __FILE__, __LINE__, #a, !!((a << 0) & 0x80), !!((a << 1) & 0x80), !!((a << 2) & 0x80), !!((a << 3) & 0x80), !!((a << 4) & 0x80), !!((a << 5) & 0x80), !!((a << 6) & 0x80), !!((a << 7) & 0x80));

void DEV_input_VCF(argStruct *args, paramStruct *pars, formulaStruct *formulaSt, IO::outFilesStruct *outSt);

void *DEV_t_EM_2DSFS_GL3(void *p);

int DEV_EM_2DSFS_GL3(threadStruct *THREAD);

void DEV_prepare_distanceMatrix_originalData(argStruct *args, paramStruct *pars, distanceMatrixStruct *dMS_orig, vcfData *vcfd, pairStruct **pairSt, formulaStruct *formulaSt, IO::outFilesStruct *outSt, sampleStruct *sampleSt);

void DEV_spawnThreads_pairEM_GL(argStruct *args, paramStruct *pars, pairStruct **pairSt, vcfData *vcfd, IO::outFilesStruct *outSt, distanceMatrixStruct *distMatrix);
