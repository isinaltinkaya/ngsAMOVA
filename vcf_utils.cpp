
#ifndef __VCF_UTILS__
#define __VCF_UTILS__

#include "shared.h"

#include <cstddef>
#include <htslib/vcf.h>


//from angsd analysisFunction.cpp
extern const int bcf_allele_charToInt[256]={
	0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//15
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//31
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//47
	0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//63
	4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//79
	4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//95
	4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//111
	4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//127
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//143
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//159
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//175
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//191
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//207
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//223
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//239
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4//255
};



int bcf_alleles_get_gtidx(int a1, int a2){
	return bcf_alleles2gt(a1,a2);
}

int bcf_alleles_get_gtidx(char a1, char a2){
	return bcf_alleles2gt(bcf_allele_charToInt[(unsigned char)a1],bcf_allele_charToInt[(unsigned char)a2]);
}

int bcf_alleles_get_gtidx(unsigned char a1, unsigned char a2){
	return bcf_alleles2gt(bcf_allele_charToInt[a1],bcf_allele_charToInt[a2]);
}

#endif
