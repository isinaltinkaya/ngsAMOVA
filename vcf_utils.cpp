#include "vcf_utils.h"





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

int VCF::read_GL10_to_GL3(bcf_hdr_t *hdr, bcf1_t *bcf, double **lngl, paramStruct *pars, argStruct *args, size_t nSites, int nInd){

	// fprintf(stderr,"\n\n\t-> Printing at site %d\n\n",nSites);
	VCF::get_data<float> lgl;

	lgl.n = bcf_get_format_float(hdr,bcf,"GL",&lgl.data,&lgl.size_e);

	if(lgl.n<0){
		fprintf(stderr,"\n[ERROR](File reading)\tVCF tag GL does not exist; will exit!\n\n");
		exit(1);
	}

	lngl[nSites]=(double*)malloc(nInd*3*sizeof(double));

	if(bcf_is_snp(bcf)){
		char a1=bcf_allele_charToInt[(unsigned char) bcf->d.allele[0][0]];
		char a2=bcf_allele_charToInt[(unsigned char) bcf->d.allele[1][0]];


		int n_missing_ind=0;
		for(int indi=0; indi<nInd; indi++){

			for(int ix=0;ix<3;ix++){
				lngl[nSites][(3*indi)+ix]=NEG_INF;
			}

			//TODO only checking the first for now
			//what is the expectation in real cases?
			//should we skip sites where at least one is set to missing?
			//
			if(isnan(lgl.data[(10*indi)+0])){
				if(args->minInd==0){
					return 1;
				}else{
					n_missing_ind++;
				}
			}else{
				lngl[nSites][(3*indi)+0]=(double) log2ln(lgl.data[(10*indi)+bcf_alleles_get_gtidx(a1,a1)]);
				lngl[nSites][(3*indi)+1]=(double) log2ln(lgl.data[(10*indi)+bcf_alleles_get_gtidx(a1,a2)]);
				lngl[nSites][(3*indi)+2]=(double) log2ln(lgl.data[(10*indi)+bcf_alleles_get_gtidx(a2,a2)]);
			}
		}

		//skip site if missing for all individuals
		if (nInd==n_missing_ind){
			return 1;
		}

		// //if there are only 2 individuals, skip site regardless of onlyShared val
		if (nInd==2){
			if(n_missing_ind>0){
				return 1;
			}
		}

		if(args->minInd>0){
			//skip site if minInd is defined and #non-missing inds=<nInd
			if( (nInd - n_missing_ind) < args->minInd ){
				// fprintf(stderr,"\n\nMinimum number of individuals -minInd is set to %d, but nInd-n_missing_ind==n_nonmissing_ind is %d at site %d\n\n",args->minInd,nInd-n_missing_ind,site);
				return 1;
			}else{
				return 0;
			}
		}else{
			return 0;
		}

		// if (set_lngl(lngl, lgl.data,args, nInd,a1,a2, nSites)==1){
		// free(lngl[nSites]);
		// lngl[nSites]=NULL;
		// // fprintf(stderr,"\nset_lngl returns 1 at site (idx: %lu, pos: %lu 1based: %lu)\n\n",nSites,bcf->pos,bcf->pos+1);
		// return 1;
		// }
	}else{
		//TODO check
		fprintf(stderr,"\n\nHERE BCF_IS_SNP==0!!!\n\n");
		free(lngl[nSites]);
		lngl[nSites]=NULL;
		return 1;
	}

	lgl.destroy();
	return 0;
}

int VCF::GT_to_i2i_SFS(bcf_hdr_t *hdr, bcf1_t *bcf, int **sfs, paramStruct *pars, argStruct *args, size_t nSites, int nInd, int** LUT_indPair_idx){

	//TODO add check if missing gt return 1

	VCF::get_data<int32_t> gt;
	gt.n = bcf_get_genotypes(hdr,bcf,&gt.data,&gt.size_e);
	gt.ploidy=gt.n/nInd;

	if(gt.n<0){
		fprintf(stderr,"\n[ERROR](File reading)\tProblem with reading GT; will exit!\n\n");
		exit(1);
	}

	if (gt.ploidy!=2){
		fprintf(stderr,"ERROR:\n\nploidy: %d not supported\n",gt.ploidy);
		exit(1);
	}

	for(int i1=0;i1<nInd-1;i1++){
		for(int i2=i1+1;i2<nInd;i2++){

			int32_t *ptr1 = gt.data + i1*gt.ploidy;
			int32_t *ptr2 = gt.data + i2*gt.ploidy;
			int pair_idx=LUT_indPair_idx[i1][i2];

			// fprintf(stderr,"\n->i1: %d, i2: %d, pair: %d, nInd: %d\n\n", i1, i2, pair_idx,nInd);
			// fprintf(stderr,"%d\t%d\t%d\n", i1, i2, pair_idx);

			int gti1=0;
			int gti2=0;

			//using binary input genotypes from VCF GT tag
			//assume ploidy=2
			for (int i=0; i<2;i++){
				gti1 += bcf_gt_allele(ptr1[i]);
				gti2 += bcf_gt_allele(ptr2[i]);
			}
			// fprintf(stderr,"\n%d:%d %d:%d -> %d\n",i1,gti1,i2,gti2,get_3x3_idx[gti1][gti2]);
			sfs[pair_idx][get_3x3_idx[gti1][gti2]]++;

		}
	}

	gt.destroy();
	return 0;

}
