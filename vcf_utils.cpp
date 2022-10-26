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


// return 1: skip site for all individuals
int VCF::read_GL10_to_GL3(bcf_hdr_t *hdr, bcf1_t *bcf, double **lngl, paramStruct *pars, argStruct *args, size_t site_i, DATA::pairStruct** PAIRS){

	// fprintf(stderr,"\n\n\t-> Printing at site %d\n\n",nSites);
	VCF::get_data<float> lgl;

	lgl.n = bcf_get_format_float(hdr,bcf,"GL",&lgl.data,&lgl.size_e);

	if(lgl.n<0){
		fprintf(stderr,"\n[ERROR](File reading)\tVCF tag GL does not exist; will exit!\n\n");
		exit(1);
	}

	lngl[site_i]=(double*)malloc(pars->nInd*3*sizeof(double));


	int* cmbArr;
	// ASSERT(pars->n_ind_cmb!=0);
	cmbArr=(int*) malloc(pars->n_ind_cmb * sizeof(int));
	ASSERT(cmbArr);
	for(int i=0; i<pars->n_ind_cmb; i++){
		cmbArr[i]=0;
	}

	//TODO check why neccessary
	if(bcf_is_snp(bcf)){
		char a1=bcf_allele_charToInt[(unsigned char) bcf->d.allele[0][0]];
		char a2=bcf_allele_charToInt[(unsigned char) bcf->d.allele[1][0]];


		for(int indi1=0; indi1<pars->nInd; indi1++){

			for(int ix=0;ix<3;ix++){
				lngl[site_i][(3*indi1)+ix]=NEG_INF;
			}

			//TODO only checking the first for now
			//what is the expectation in real cases?
			//should we skip sites where at least one is set to missing?
			//
			if(isnan(lgl.data[(10*indi1)+0])){

				//if only use sites shared across all individuals; skip site when you first encounter nan
				if(args->minInd==0){ 
					return 1;
				}else{

					//if there are only 2 individuals any missing will skip the site 
					if (pars->nInd==2){
						return 1;
					}

					lgl.n_missing_ind++;

					if (pars->nInd==lgl.n_missing_ind){
						return 1;
					}
					

					if(args->minInd!=2){
						//skip site if minInd is defined and #non-missing inds=<nInd
						if( (pars->nInd - lgl.n_missing_ind) < args->minInd ){
							// fprintf(stderr,"\n\nMinimum number of individuals -minInd is set to %d, but nInd-n_missing_ind==n_nonmissing_ind is %d at site %d\n\n",args->minInd,pars->nInd-n_missing_ind,site);
							return 1;
						}
					}

				}

			}else{

				//if not first individual, check previous individuals pairing with current ind
				if(indi1!=0){
					for (int indi2=indi1-1; indi2>-1; indi2--){
						if(cmbArr[pars->LUT_indPair_idx[indi2][indi1]]==1){
							//both inds has data

							int pidx=pars->LUT_indPair_idx[indi2][indi1];

							//append site to sSites shared sites list
							PAIRS[pidx]->sSites[PAIRS[pidx]->snSites]=site_i;

							PAIRS[pidx]->snSites++;

							if(PAIRS[pidx]->snSites == PAIRS[pidx]->_sSites){
								PAIRS[pidx]->_sSites = 2 * PAIRS[pidx]->_sSites;
								PAIRS[pidx]->sSites=(int*) realloc(PAIRS[pidx]->sSites, PAIRS[pidx]->_sSites * sizeof(int));
							}

						
						}
					}

				}

				//if not last individual, check latter individuals pairing with current ind
				if(indi1!= pars->nInd-1){
					for (int indi2=indi1+1; indi2<pars->nInd; indi2++){
						cmbArr[pars->LUT_indPair_idx[indi1][indi2]]++;
					}
				}

				lngl[site_i][(3*indi1)+0]=(double) LOG2LN(lgl.data[(10*indi1)+bcf_alleles_get_gtidx(a1,a1)]);
				lngl[site_i][(3*indi1)+1]=(double) LOG2LN(lgl.data[(10*indi1)+bcf_alleles_get_gtidx(a1,a2)]);
				lngl[site_i][(3*indi1)+2]=(double) LOG2LN(lgl.data[(10*indi1)+bcf_alleles_get_gtidx(a2,a2)]);
			}
		}


	}else{
		//TODO check
		fprintf(stderr,"\n\nHERE BCF_IS_SNP==0!!!\n\n");
		exit(1);
		// free(lngl[nSites]);
		// lngl[nSites]=NULL;
		// return 1;
	}

	free(cmbArr);
	return 0;
}



//minInd is already checked before this
//
int VCF::GL_to_GT_1_SFS(bcf_hdr_t *hdr, bcf1_t *bcf, int **sfs, paramStruct *pars, argStruct *args){

	// fprintf(stderr,"\n\n\t-> Printing at site %d\n\n",nSites);
	VCF::get_data<float> lgl;

	lgl.n = bcf_get_format_float(hdr,bcf,"GL",&lgl.data,&lgl.size_e);

	if(lgl.n<0){
		fprintf(stderr,"\n[ERROR](File reading)\tVCF tag GL does not exist; will exit!\n\n");
		exit(1);
	}


	VCF::get_data<int32_t> new_gt;

	// hts_expand(int32_t, pars->nInd*2, new_gt.n, new_gt.data);
	// new_gt.ploidy=new_gt.n/pars->nInd;
	new_gt.ploidy=2;

	if (new_gt.ploidy!=2){
		fprintf(stderr,"ERROR:\n\nploidy: %d not supported\n",new_gt.ploidy);
		exit(1);
	}



	if(bcf_is_snp(bcf)){


		for(int indi=0; indi<pars->nInd; indi++){


			//pointer to genotype likelihoods
			float *pgl = lgl.data + indi*new_gt.ploidy;

			if(bcf_float_is_missing(pgl[0])){
				new_gt.data[2*indi]=new_gt.data[2*indi+1]=bcf_gt_missing;
				continue;
			}
			if(isnan(lgl.data[(10*indi)+0])){
				new_gt.data[2*indi]=new_gt.data[2*indi+1]=bcf_gt_missing;
				continue;
			}


			float max_like=NEG_INF;
			int max_idx=-1;
			//TODO
			//max is 0 in simulated data input due to rescaling
			//ignore multiple 0
			//
			for(int i=0; i<10; i++){
				if(lgl.data[(10*indi)+i] >= max_like){
					max_like=lgl.data[(10*indi)+i];
					max_idx=i;
				}
			}

			int a,b;
			bcf_gt2alleles(max_idx, &a, &b);
			new_gt.data[2*indi] = bcf_gt_unphased(a);
			new_gt.data[2*indi+1] = bcf_gt_unphased(b);

		}
	}

	for(int i1=0;i1<pars->nInd-1;i1++){
		for(int i2=i1+1;i2<pars->nInd;i2++){

			int32_t *ptr1 = new_gt.data + i1*new_gt.ploidy;
			int32_t *ptr2 = new_gt.data + i2*new_gt.ploidy;
			int pair_idx=pars->LUT_indPair_idx[i1][i2];

			int gti1=0;
			int gti2=0;

			//using binary input genotypes from VCF GT tag
			//assume ploidy=2
			for (int i=0; i<2;i++){
				gti1 += bcf_gt_allele(ptr1[i]);
				gti2 += bcf_gt_allele(ptr2[i]);
			}
			sfs[pair_idx][get_3x3_idx[gti1][gti2]]++;

			//last field is for snSites
			sfs[pair_idx][9]++;

		}
	}

	return 0;
}


//if doAMOVA==3; both gl and gt; if gl returns 1 to skip all sites, this never runs
//if gl returns 0; this will check again for shared sites using DP tag as an indicator
// return 1: skip site for all individuals
//
// sfs[pair_idx][9] holds snSites

int VCF::GT_to_i2i_SFS(bcf_hdr_t *hdr, bcf1_t *bcf, int **sfs, paramStruct *pars, argStruct *args){
// fprintf(stderr,"\n\nHERE!!!\n\n");

	int nInd=pars->nInd;

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


	const char* TAG="DP";
	VCF::get_data<int32_t> dp;
	dp.n = bcf_get_format_int32(hdr,bcf,TAG,&dp.data,&dp.size_e);
	if(dp.n<0){
		fprintf(stderr,"\n[ERROR](File reading)\tVCF tag \"%s\" does not exist; will exit!\n\n",TAG);
		exit(1);
	}

	if(args->minInd>2){
		for(int indi=0; indi<nInd; indi++){
			if(dp.data[indi]==0){
				dp.n_missing_ind++;
			}
		}
		if ( (nInd - dp.n_missing_ind) < args->minInd ){
fprintf(stderr,"\n\nHERE!!!\n\n");
			return 1;
		}
	}


	for(int i1=0;i1<nInd-1;i1++){
		for(int i2=i1+1;i2<nInd;i2++){

			if(dp.data[i1]==0 || dp.data[i2]==0){
				//skip the pair
				continue;
			}

			int32_t *ptr1 = gt.data + i1*gt.ploidy;
			int32_t *ptr2 = gt.data + i2*gt.ploidy;
			int pair_idx=pars->LUT_indPair_idx[i1][i2];

			int gti1=0;
			int gti2=0;

			//using binary input genotypes from VCF GT tag
			//assume ploidy=2
			for (int i=0; i<2;i++){
				gti1 += bcf_gt_allele(ptr1[i]);
				gti2 += bcf_gt_allele(ptr2[i]);
			}
			sfs[pair_idx][get_3x3_idx[gti1][gti2]]++;

			//last field is for snSites
			sfs[pair_idx][9]++;

		}
	}

	return 0;
}




