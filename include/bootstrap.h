/**
 * @file    bootstrap.h
 * @brief   header file for bootstrap.cpp
 * @details contains block bootstrapping related data structures and functions
 */
#ifndef __BOOTSTRAP_H__
#define __BOOTSTRAP_H__

/* INCLUDES ----------------------------------------------------------------- */
#include "dataStructs.h"
#include "shared.h"
#include <htslib/hts.h>  // for hts_idx_t
/* END-OF-INCLUDES ---------------------------------------------------------- */


/* FORWARD-DECLARATIONS ----------------------------------------------------- */
typedef struct bblocks_t bblocks_t;
typedef struct strArray strArray;
typedef struct vcfData vcfData;
typedef struct paramStruct paramStruct;
typedef struct outfile_t outfile_t;
/* END-OF-FORWARD-DECLARATIONS ---------------------------------------------- */

/* MACROS ------------------------------------------------------------------- */
/* END-OF-MACROS ------------------------------------------------------------ */

/* TYPEDEF-STRUCTS ---------------------------------------------------------- */

/// @brief bblocks_t - structure for storing block bootstrapping blocks
/// @note positions: 0-based, [start, end)
struct bblocks_t {

    /// @var n_blocks - number of blocks
    size_t n_blocks;

    /// @var n_contigs - number of contigs
    size_t n_contigs;

    /// @var n_ind - number of individuals
    size_t n_ind;

    /// \def rblocks[nReps][n_blocks] = size_t index of the sampled block
    size_t** rblocks;

    /// @var nsites_per_block[n_blocks]
    /// nsites_per_block[i] == number of sites in block i
    uint64_t* nsites_per_block;

    strArray* contig_names;

    /// @var block_start_pos[n_blocks]
    /// block_start_pos[i] == position (vcfd->rec->pos) of the start of block i
    /// @note 0-based, inclusive [start
    size_t* block_start_pos;

    /// @var block_end_pos[n_blocks]
    /// block_end_pos[i] == position (vcfd->rec->pos) of the end of block i
    /// @note 0-based, exclusive end)
    size_t* block_end_pos;

    /// @var block_contig[n_blocks]
    /// block_contig[i] == index of the contig (in contig_names) that block i belongs to
    size_t* block_contig;

    /// @var block_start_siteidx[n_blocks]
    /// block_start_siteidx[i] == index of the first site that belongs to block i
    /// @note site index refers to the set of sites the program did not skip and are included in the analysis
    /// if all sites from all contigs are used, this set is of size for(c in contigs) size+=contignsites[c]
    /// the max siteidx is ((vcfd->nSites)-1), the min is 0, and it is global through all contigs
    /// vcfd->nSites==pars->nSites
    /// this information is collected during sites reading
    /// @note 0-based siteidx inclusive [start
    /// @note block_end = (wb == nBlocks - 1) ? max_nsites : block_start_siteidx[wb + 1];
    size_t* block_start_siteidx;

};
/* END-OF-TYPEDEF-STRUCTS --------------------------------------------------- */

/* FUNCTION-DECLARATIONS ----------------------------------------------------- */

/// @brief bblocks_sample_with_replacement - perform the block bootstrapping sampling with replacement
/// @details
/// sample the blocks with replacement for all individuals as to preserve the correlation structure
void bblocks_sample_with_replacement(bblocks_t* bblocks);

void bblocks_print_bootstrap_samples(bblocks_t* bblocks, outfile_t* outfile);
void bblocks_print_blocks_tab(bblocks_t* bblocks, outfile_t* outfile);

bblocks_t* bblocks_init(void);

void bblocks_destroy(bblocks_t* bblocks);

void bblocks_get(bblocks_t* bblocks, vcfData* vcfd, paramStruct* pars);

/* END-OF-FUNCTION-DECLARATIONS ---------------------------------------------- */


#endif  // __BOOTSTRAP_H__