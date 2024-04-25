/**
 * @file    em.h
 * @brief   header file for em.cpp
 * @details contains functions for performing expectation-maximization (EM) optimization
 */
#ifndef __EM_H__
#define __EM_H__

/* INCLUDES ----------------------------------------------------------------- */
/* END-OF-INCLUDES ---------------------------------------------------------- */

/* FORWARD-DECLARATIONS ----------------------------------------------------- */
typedef struct jgtmat_t jgtmat_t;
typedef struct paramStruct paramStruct;
typedef struct gldata_t gldata_t;
typedef struct bblocks_t bblocks_t;
/* END-OF-FORWARD-DECLARATIONS ---------------------------------------------- */

/* MACROS ------------------------------------------------------------------- */
// #define EM_TERM_REASON_THRES_TOLE 0
// #define EM_TERM_REASON_THRES_MAXITER 1
// #define EM_TERM_REASON_THRES_SHARED_NSITES 2

typedef enum {
    EM_TERM_REASON_THRES_TOLE = 0,
    EM_TERM_REASON_THRES_MAXITER,
    EM_TERM_REASON_THRES_SHARED_NSITES,
    EM_TERM_REASON_COUNT 
} em_term_reason;

extern const char* em_term_reason_strs[EM_TERM_REASON_COUNT];
const char* get_em_term_reason_str(em_term_reason reason);


/* END-OF-MACROS ------------------------------------------------------------ */

/* TYPEDEF-STRUCTS ---------------------------------------------------------- */
/* END-OF-TYPEDEF-STRUCTS --------------------------------------------------- */

/* FUNCTION-DECLARATIONS ----------------------------------------------------- */
void jgtmat_get_run_em_optim(jgtmat_t* jgtmat, paramStruct* pars, gldata_t* gldata, bblocks_t* bblocks);
/* END-OF-FUNCTION-DECLARATIONS ---------------------------------------------- */

#endif  // __EM_H__