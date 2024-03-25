#ifndef __EM__
#define __EM__

#include <pthread.h>

#include "argStruct.h"
#include "mathUtils.h"

#define EM_TERM_REASON_TOLE 1
#define EM_TERM_REASON_ITER 2

void jgtmat_get_run_em_optim(jgtmat_t* jgtm, paramStruct* pars, vcfData* vcfd, bblocks_t* bblocks);
void jgtmat_get_run_em_optim_bootstrap_reps(jgtmat_t* jgtm, paramStruct* pars, vcfData* vcfd, bblocks_t* bblocks);

#endif  // __EM__