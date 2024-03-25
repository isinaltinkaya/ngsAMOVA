#ifndef __PARAM_STRUCT__
#define __PARAM_STRUCT__

#include "shared.h"


extern const int acgt_charToInt[256];

/// ----------------------------------------------------------------------- ///

typedef struct strArray strArray;
typedef struct size_tArray size_tArray;
typedef struct ibdStruct ibdStruct;
typedef struct dmat_t dmat_t;
typedef struct jgtmat_t jgtmat_t;
typedef struct paramStruct paramStruct;
typedef struct alleles_t alleles_t;
typedef struct metadataStruct metadataStruct;

/// ----------------------------------------------------------------------- ///

void test_alleles_t();


/// @brief alleles_t - allele pair structure
/// @details 
///   
///   nContigs: cnames->size
///   nSites: pos->size
///   number of 64bit blocks: nSites/16 + 1
///   nSites in the last block64: nSites % 16
///
struct alleles_t {

    /// @var  d - per-site alleles data
    /// @size (nSites/16) + 1
    /// @details
    ///     - contains ordered a1a2 for each site
    ///       e.g. major/minor, or ancestral/derived, or reference/alternate alleles
    ///     - each uint64_t contains 16 packs of 4bit alleles information
    ///     - alleles are encoded in 4bits (ordered pairs of 2bit ACGTs)
    ///
    ///       base | 2-bit encoding
    ///       A    | 00
    ///       C    | 01
    ///       G    | 10
    ///       T    | 11
    ///
    ///       e.g. major: C, minor: T = [ major2bits, minor2bits ] = [ 0 1 , 1 1 ]
    ///
    uint64_t* d;

    /// @var  pos - 0-indexed position
    /// @size nSites
    size_tArray* pos;

    /// @var  cposidx - indices of contig start positions in pos array
    /// @size nContigs
    /// @details
    /// for each contig in cnames, the index of the first position belonging to that contig in pos array
    /// e.g. cnames={"chr1","chr2"}
    /// pos={0,42,64,100,20,21,23}  // 0-based positions
    ///      0, 1, 2,  3, 4, 5, 6   // indices of pos array
    ///      |chr1      ,|chr2 
    /// chr1 pos={0,42,64,100}
    /// chr2 pos={20,21,23}
    /// cposidx={0,4}  // 0-based index of first position of each contig in pos array
    size_tArray* cposidx;

    /// @var  cnames - contig names
    /// @size nContigs
    strArray* cnames;

};





/// ----------------------------------------------------------------------- ///

/*
 * @typedef paramStruct - parameter structure
 *
 * - get alleles at a given position and contig
 *    char a1, a2;
 *    size_t querypos = 1;
 *    char querychr[1024] = "chr23";
 *    alleles_get(pars->majorminor, querychr, querypos, &a1, &a2);
 *    DEVPRINT("Query result:\n%s\t%ld\t%c\t%c\n", querychr, querypos, a1, a2);
 *
 * - get alleles, position and contig name at a given position index
 *    size_t queryposidx = 4;
 *    char* retcontig = (char*)malloc(1024 * sizeof(char));
 *    size_t retpos;
 *    alleles_get(pars->majorminor, queryposidx, &retcontig, &retpos, &a1, &a2);
 *    DEVPRINT("Query result:\n%s\t%ld\t%c\t%c\n", retcontig, retpos, a1, a2);
 *    FREE(retcontig);
 *
 */
struct paramStruct {


    dmat_t* dm;
    jgtmat_t* jgtm;

    strArray* names;


    metadataStruct* metadata;

    // number of sites non skipped for all individuals
    // nSites may not be !=totSites if minInd is set
    // or if a site is missing for all inds
    size_t nSites;    // number of sites not-skipped for all individuals
    size_t totSites;  // total number of sites processed

    int nSites_arrays_size;
    int nContigs;

    // a1: major/ref/ancestral allele
    // a2: minor/alt/derived allele
    // a1a2[site] = {base_a1,base_a2}
    // where base_* is 0 for A, 1 for C, 2 for G, 3 for T, 4 for BASE_UNOBSERVED


    alleles_t* majorminor; // major/minor alleles for each site
    alleles_t* ancder;     // ancestral/derived alleles for each site

    // at each new vcf record, alleles_posidx contains last position we used
    // for search, start from alleles_posidx instead of alleles_posidx+1 to also handle the first ever contig
    int64_t alleles_posidx;  // current position in the alleles_t->pos
    int64_t alleles_contigidx;  // current contig in the alleles_t->cnames

    // for each site, 
    // a1a2[0] is the index of the allele1/major/ancestral allele in vcfd->rec->d.allele
    // a1a2[1] is the index of the allele2/minor/derived allele in vcfd->rec->d.allele
    int a1a2[2];


    /// @var formulaTokens
     // array of tokens in the formula, in hierarchical order
     // e.g. formula: "Individual ~ Region/Population/Subpopulation"
     // 		formulaTokens = {"Region","Population","Subpopulation","Individual"}
     //                   lvl =   1      , 2          , 3             , 0
     //                lvlidx =   0      , 1          , 2             , 3
    strArray* formulaTokens;

    // ------------
    ibdStruct* ibd;
    // ------------

    char* DATETIME;


};


/// @brief paramStruct_init initialize the paramStruct
/// @param args arguments argStruct
/// @return pointer to paramStruct
paramStruct* paramStruct_init(argStruct* args);
void paramStruct_destroy(paramStruct* p);

void check_consistency_args_pars(paramStruct* pars);


#endif  // __PARAM_STRUCT__
