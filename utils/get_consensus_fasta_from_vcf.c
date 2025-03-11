#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>

#define ASSERT_EXPAND(x) x
#define ASSERT(expr)                                                        \
    do {                                                                    \
        if (!(ASSERT_EXPAND(expr))){                                        \
            fprintf(stderr, "\n\n*******\n[ERROR](%s:%d) %s\n*******\n", \
                    __FILE__, __LINE__, #expr);               \
            exit(1);                                                        \
        }                                                                   \
    } while(0)

#define ERROR(...)  \
    do{ \
        fprintf(stderr, "\n\n*******\n[ERROR](%s:%d)\n\t", __FILE__, __LINE__); \
        fprintf(stderr, __VA_ARGS__); \
        fprintf(stderr, "\n*******\n"); \
        exit(1); \
    } while(0)



#define INITIAL_BUFFER 1024

typedef struct {
	char* seq;
	int len;
	int cap;
} Sequence;

int generate_consensus(const char* vcf_filename, const char* output_filename){
	htsFile* vcf = hts_open(vcf_filename, "r");
	if (!vcf) {
		fprintf(stderr, "Error opening VCF file %s\n", vcf_filename);
		return 1;
	}
	bcf_hdr_t* hdr = bcf_hdr_read(vcf);
	if (!hdr) {
		fprintf(stderr, "Error reading header from VCF file\n");
		hts_close(vcf);
		return 1;
	}
	int n_samples = bcf_hdr_nsamples(hdr);


	Sequence* consensus = malloc(n_samples * sizeof(Sequence));
	ASSERT(consensus!=NULL);
	for (int i = 0; i < n_samples; i++) {
		consensus[i].cap = INITIAL_BUFFER;
		consensus[i].len = 0;
		consensus[i].seq = malloc(consensus[i].cap * sizeof(char));
		ASSERT(consensus[i].seq!=NULL);
	}

	bcf1_t* rec = bcf_init();
	ASSERT(rec!=NULL);

	int ngt_arr = 0;
	int32_t* gt_arr = NULL;
	while (bcf_read(vcf, hdr, rec) == 0) {
		bcf_unpack(rec, BCF_UN_ALL);
		int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
		ASSERT(ngt>0);

		int ploidy = ngt / n_samples;
		ASSERT(ploidy==2);

		int32_t* sample_gts=NULL;

		for (int i = 0; i < n_samples; i++) {

			sample_gts = gt_arr + i * ploidy;
			int32_t allele1,allele2;
			char allele_char;
			if (bcf_gt_is_missing(sample_gts[0]) || bcf_gt_is_missing(sample_gts[1])) {
				allele_char = 'N';
			} else {
				allele1 = bcf_gt_allele(sample_gts[0]);
				allele2 = bcf_gt_allele(sample_gts[1]);
				int chosen_allele;
				if (allele1 == allele2) {
					chosen_allele = allele1;
				} else {
					chosen_allele = (rand() % 2) ? allele2 : allele1;
				}
				ASSERT(strlen(rec->d.allele[chosen_allele])==1);
				allele_char=rec->d.allele[chosen_allele][0];
			}
			if (consensus[i].len + 1 >= consensus[i].cap) {
				consensus[i].cap *= 2;
				consensus[i].seq = realloc(consensus[i].seq, consensus[i].cap * sizeof(char));
				ASSERT(consensus[i].seq!=NULL);
			}
			ASSERT(consensus[i].len>=0);
			ASSERT(consensus[i].len<consensus[i].cap);
			consensus[i].seq[consensus[i].len] = allele_char;
			consensus[i].len++;
		}
	}
	for (int i = 0; i < n_samples; i++) {
		if (consensus[i].len + 1 >= consensus[i].cap) {
			consensus[i].cap = consensus[i].len + 2;
			consensus[i].seq = realloc(consensus[i].seq, consensus[i].cap * sizeof(char));
			ASSERT(consensus[i].seq!=NULL);
		}
		consensus[i].seq[consensus[i].len] = '\0';
	}

	FILE* out = fopen(output_filename, "w");
	for (int i = 0; i < n_samples; i++) {
		fprintf(out, ">%s\n%s\n", hdr->samples[i], consensus[i].seq);
	}
	fclose(out);

	free(gt_arr);
	gt_arr=NULL;
	bcf_destroy(rec);
	bcf_hdr_destroy(hdr);
	hts_close(vcf);
	for (int i = 0; i < n_samples; i++) {
		free(consensus[i].seq);
		consensus[i].seq=NULL;
	}
	free(consensus);
	consensus=NULL;
	return 0;

}

int main(int argc, char* argv[]) {
	if (argc < 4) {
		fprintf(stderr, "Usage: %s <seed> <vcf_file> <output_file>\n", argv[0]);
		return 1;
	}
	unsigned int seed = (unsigned int)atoi(argv[1]);
	srand(seed);
	return generate_consensus(argv[2], argv[3]);
}
