#ifndef __IBD__
#define __IBD__

#include "shared.h"

#include "dataStructs.h"
#include "io.h"
#include "paramStruct.h"
#include "vcfReader.h"

struct paramStruct;

typedef struct ibds_t {

	/// @var size - current allocated n_sites size 
	/// @details at runtime, worst case scenario: total number of sites in the longest contig
	/// @note at the end_posidx of read_sites, size >= pars->nSites
	size_t size;

	/// @var step_size - step size for reallocating/expanding the arrays that were allocated using size
	/// @note step_size is the added size when reallocating
	//TODO RM?
	size_t step_size;

//TODO RM?
	size_t n_pairs;

	/// @var trim - indicator for trimming segments based on cumulative ibd scores
	/// @details true iff args->ibd_ibdtrim != 0.0
	bool trim;

	/// @var contig_name - current contig name
	/// @note at contig change, this is set to the latest contig name so not the same as vcfd->get_contig_name()
	kstring_t ks_contig_name;

	/// @var pos0 - 0-indexed positions (rec->pos) for each non-skipped site in the current contig
	hts_pos_t* pos0;

	/// @var pairs_ibd_scores[n_pairs][n_sites_in_contig] - array of pairwise per-site IBD scores for a contig
	/// @details lifetime: reset on contig change
	double** pairs_ibd_scores;

	outfile_t* out_ibd_segments;
	outfile_t* out_persite_ibd_scores;
	outfile_t* out_persite_smoothed_ibd_scores;

} ibds_t;

ibds_t* ibds_alloc(const size_t n_pairs, const size_t init_n_sites, const size_t init_step_size);
ibds_t* ibds_realloc(ibds_t* ibds);
void ibds_destroy(ibds_t* ibds);

size_t ibds_estimate_max_mem_use(const size_t n_pairs, const size_t max_percontig_n_sites);

void ibds_print(ibds_t* ibds);
void ibds_reset_percontig(ibds_t* ibds);

void identify_ibd_segments(ibds_t* ibds, paramStruct* pars);
void identify_ibd_segments_ibdseq(ibds_t* ibds, paramStruct* pars);



// --- ibdseq method

double get_ibd_likelihood(int dose1, int dose2, double e, double pfa2);
double get_null_likelihood(int dose1, int dose2, double fB);

#endif // __IBD__