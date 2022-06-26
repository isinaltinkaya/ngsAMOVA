#ifndef __VCF_UTILS__
#define __VCF_UTILS__

#include "shared.h"

#include <cstddef>

int bcf_alleles_get_gtidx(int a1, int a2);

int bcf_alleles_get_gtidx(char a1, char a2);

extern const char bcf_allele_charToInt[256];


#endif
