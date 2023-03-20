/*
 *
 * [Parameters]
 *
 * idea from angsd funkyPars
 * at https://github.com/ANGSD/angsd/blob/master/shared.cpp
 *
 */

#include "shared.h"
#include "mathUtils.h"
#include <stdlib.h>

#include "argStruct.h"
#include "paramStruct.h"



//
// 0 1 2
// 00 01 02
// MMMM MMMm MMmm
//
// 3 4 5
// 10 11 12
// MmMM MmMm Mmmm
//
// 6 7 8
// 20 21 22
// mmMM mmMm mmmm
// TODO rename and consider change format
extern const int get_3x3_idx[3][3] = {
	{0, 1, 2},
	{3, 4, 5},
	{6, 7, 8}};

// TODO check this
using size_t = decltype(sizeof(int));

int find_n_given_nC2(int nC2_res)
{
	int n = 0;
	while (NC2_LUT[n] < nC2_res)
	{
		n++;
	}
	if (NC2_LUT[n] != nC2_res)
	{
		fprintf(stderr, "[%s:%s()]\t->Error: nC2_res:%d not found in NC2_LUT[]\n", __FILE__, __FUNCTION__, nC2_res);
		exit(1);
	}
	return n;
}

// TODO
//  extract digits using bit masking
//  int extractDigits(int num, int digits)
//  {
//  	uint32_t x = 100101;
//  	uint32_t y = 100102;
//  	for(int i=0;i<sizeof(uint32_t);i++){
//  		fprintf(stderr,"\n\n\n\n -> x: %d y: %d\n",(x>>i)&0xF,(y>>i)&0xF);
//  	}
//  	exit(0);
//  }

/// @brief get current time
/// @return time as char*
char *get_time()
{
	time_t current_time;
	struct tm *local_time;
	current_time = time(NULL);
	local_time = localtime(&current_time);
	return (asctime(local_time));
}

