#include "ibd.h"

#include  "vcfReader.h"


//////////////////////////
///
///
///// if (args->minInd ) does not work with doibd
///
//////////////////////////

const int dosePairIndex[3][3] = { {0,1,5}, {1,2,3}, {5,3,4} };

void readSites_doIbd(vcfData *vcfd, paramStruct *pars){

	const int nInd=pars->nInd;
	int skip_site = 0;

	int contig_i = -1;
	int site_i = 0;

	char prev_contig[1024];
	int nSites_contig_i = 0;

	int* site2pos=(int*)malloc(BUF_NSITES*sizeof(int));
	for (int i=0;i<BUF_NSITES;++i){
		site2pos[i]=-1;
	}

	double* lngl=NULL;
	if(args->doIbd==1){
		lngl=(double*) malloc(nInd*3*sizeof(double));
		for (int indi = 0; indi < nInd; indi++) {
			const int lngls_ind_start = indi * 3;
			for (int j = 0; j < 3; j++) {
				lngl[lngls_ind_start + j] = NEG_INF;
			}
		}
	}

	double* ibdScores_lut = NULL;
	double* hbdScores_lut = NULL;
	if(args->doIbd==2){
		ibdScores_lut = (double*) malloc(6 * sizeof(double));
		hbdScores_lut = (double*) malloc(3 * sizeof(double));

	}


	int *a1 = new int;
	int *a2 = new int;


	while (1 == (vcfd->records_next())) {


		ASSERT(site_i<BUF_NSITES-1);//TODO


		if(NULL != lngl){

			for (int indi = 0; indi < nInd; indi++) {
				const int lngls_ind_start = indi * 3;
				for (int j = 0; j < 3; j++) {
					lngl[lngls_ind_start + j] = NEG_INF;
				}
			}
		}

		// if at the very beginning
		if (-1 == contig_i) {
			++contig_i;
			strncpy(prev_contig, vcfd->get_contig_name(), sizeof(prev_contig));
			if (prev_contig[sizeof(prev_contig) - 1] != '\0') {
				ERROR("Contig name is too long");
			}
			site_i = 0;
		}


		// ** contig change **
		if (0 != strcmp(vcfd->get_contig_name(), prev_contig)) {
			if (contig_i >= vcfd->nContigs) {
				ERROR("Number of contigs in the VCF file header is larger than the number of contigs in the VCF file");
			}

			// if we skipped all sites in the previous contig, reuse contig_i index for the new contig
			if (0 == nSites_contig_i) {
			} else {
				++contig_i;
				nSites_contig_i = 0;
			}

			strncpy(prev_contig, vcfd->get_contig_name(), sizeof(prev_contig));
			if (prev_contig[sizeof(prev_contig) - 1] != '\0') {
				ERROR("Contig name is too long");
			}
			site_i = 0;
		}


		ASSERT(1 == bcf_is_snp(vcfd->rec));
		// if (1 != bcf_is_snp(vcfd->rec)) {
			// skip_site=1;
			// IO::vprint(2, "Skipping site at %s:%ld. Reason: VCF record is not a SNP.\n", vcfd->get_contig_name(), vcfd->rec->pos + 1);
		// }

		ASSERT(0== read_site_with_alleles(contig_i, site_i, vcfd, pars, a1, a2));
		// int ret = read_site_with_alleles(contig_i, site_i, vcfd, pars, a1, a2);
		// if (0 == ret) {
			// skip_site = 0;
		// } else if (1 == ret) {
			// // 1    if 'contig:site' in VCF cannot be found in alleles file
			// skip_site = 1;
			// break;
		// } else if (2 == ret) {
			// // 2    if an a1 is not a single character
			// skip_site = 1;
			// break;
		// } else {
			// NEVER;
		// }


		get_data<float> afData;
		int n = bcf_get_info_float(vcfd->hdr, vcfd->rec, "AF", &afData.data, &afData.n);
		if(n<=0){
			ERROR("AF tag not found.");
		}
		double maf = afData.data[0];
		if(maf>0.5){
			NEVER;
		}

		// ------------------------------------------------
		// FILTER: --min-af <DOUBLE>
		if (args->min_af>0){
			if(maf<args->min_af){
				skip_site=1;
				NEVER;
				break;
			}
		}
		// ------------------------------------------------

		skip_site=0;


		double score=0.0;


		/// -------------------
		if(args->doIbd==1){

			glData lgls(vcfd);
			bool ind_is_missing[nInd]={false};

			do {

				//TODO assume only 3 gls for now
				// int a1a1 = bcf_alleles_get_gtidx(*a1, *a1);D
				// int a1a2 = bcf_alleles_get_gtidx(*a1, *a2);
				// int a2a2 = bcf_alleles_get_gtidx(*a2, *a2);
				int a1a1=0;
				int a1a2=1;
				int a2a2=2;


				for (int indi = 0; indi < nInd; indi++) {
					if (1 == lgls.ind_data_isMissing(indi)) {
						ind_is_missing[indi]=true;
					}else{
						const int lngls_ind_start = indi * vcfd->nGT;

						if((lgls.ind_ptr(indi)[a1a1]==lgls.ind_ptr(indi)[a1a2]) && (lgls.ind_ptr(indi)[a1a1]==lgls.ind_ptr(indi)[a2a2]) ) {

							ind_is_missing[indi]=true;
						}else{

							lngl[lngls_ind_start + 0] = (double)LOG2LN(lgls.ind_ptr(indi)[a1a1]);
							lngl[lngls_ind_start + 1] = (double)LOG2LN(lgls.ind_ptr(indi)[a1a2]);
							lngl[lngls_ind_start + 2] = (double)LOG2LN(lgls.ind_ptr(indi)[a2a2]);


							double mmax=lngl[lngls_ind_start+0];
							for(int ii=1;ii<3;ii++){
								if(lngl[lngls_ind_start+ii]>mmax){
									mmax=lngl[lngls_ind_start+ii];
								}
							}
							for(int ii=0;(!std::isinf(mmax))&&ii<3;ii++){
								lngl[lngls_ind_start+ii]-=mmax;
								if(std::isnan(lngl[lngls_ind_start+ii])){
									fprintf(stderr,"mmax: %f\n",mmax);
									exit(0);
								}
							}

						}
					}
				}

									  //
				bool isCor=false;//TODO
				score=0.0;

				for (int i1 = 0; i1 < nInd; ++i1) {
					for (int i2 = i1+1; i2 < nInd; ++i2) {

						if(ind_is_missing[i1]==true||ind_is_missing[i2]==true){
							///if any ind in the pair is missing define LOD score of the pair at the locus as 0
							score=0.0;
						}else{
							if (isCor==false){
								score=pars->ibd->ibdScore_GL(i1,i2,vcfd, pars, maf,site_i,lngl);
							}else{
								NEVER;
								// score=0.0;
							}
						}

						int pair_idx=nCk_idx(nInd,i1,i2);
						pars->ibd->pairScores[pair_idx][site_i]=score;
						// VALIDATION GL1
						// fprintf(stdout, "VGL1,%d,%d,%d,%f\n",i1,i2,site_i,score);
						// fprintf(stdout, "VGL1,%d,%d,%d,%f\n",i1,i2,vcfd->rec->pos+1,score);


					}
				}

			} while (0);




		}else if (args->doIbd==2){




				// precalculate the dose pair ibd scores
				//
				// getIbdScores
				int nIndices = 6;
				int N_DOSES=3;
				int index=-1;
				for (int dose1=0; dose1<N_DOSES; ++dose1) {
					for (int dose2=dose1; dose2<N_DOSES; ++dose2) {
						ibdScores_lut[dosePairIndex[dose1][dose2]]=pars->ibd->ibdScore(dose1, dose2, maf);
						ibdScores_lut[dosePairIndex[dose2][dose1]]=pars->ibd->ibdScore(dose1, dose2, maf);
					}
				}

				// getHbdScores
				for (int dose=0; dose<N_DOSES; ++dose) {
					hbdScores_lut[dose] = pars->ibd->hbdScore(dose, maf);
				}


			do {

				vcfd->site_gts_get(*a1, *a2);
				if (!vcfd->gts->pass_minInd_threshold(nInd)) {
					skip_site = 1;
					break;
				}


				for (int ind = 0; ind < nInd; ++ind) {

					int ind_nder = vcfd->gts->get_n_derived_alleles_ind(ind);
					if (-2 == ind_nder) {
						// individual's allelic state cannot be found in d.allele
						NEVER;
					}
					if (-1 == ind_nder) {
						continue;
					}
				}



				// ------------------------------------------------
				// FILTER: --ibdseq-minalleles <INT>
				// 
				bool has_ac_tag=false;
				// First, if exists try using AC tag for this filter
				// AC : allele count in genotypes, for each ALT allele, in the same order as listed
				get_data<int32_t> acData;
				int nac = bcf_get_info_int32(vcfd->hdr, vcfd->rec, "AC", &acData.data, &acData.n);
				if(nac<=0){
					has_ac_tag=false;
					// ERROR("AC tag not found.");
				}else{
					has_ac_tag=true;
					if(acData.data[0] < args->ibdseq_minalleles){
						IO::vprint(2, "Skipping site at %s:%ld. Reason: --ibdseq-minalleles filter is not met.");
						skip_site=1;
						break;
					}
				}
				if(false==has_ac_tag){
					int sum_nder=0;
					for (int ind = 0; ind < nInd; ++ind) {
						int ind_nder = vcfd->gts->get_n_derived_alleles_ind(ind);
						if (ind_nder!=-1){
							sum_nder+=ind_nder;
						}
					}
					if(sum_nder < args->ibdseq_minalleles){
						IO::vprint(2, "Skipping site at %s:%ld. Reason: --ibdseq-minalleles filter is not met.");
						skip_site=1;
						break;
					}
				}
				// ------------------------------------------------

				score=0.0;
				bool isCor=false;//TODO

				for (int i1 = 0; i1 < nInd; ++i1) {
					for (int i2 = i1+1; i2 < nInd; ++i2) {

						if(i1==i2){
							NEVER;
//
							// int dose=vcfd->gts->get_n_derived_alleles_ind(i1);
//
							// ASSERT(dose>=0);
//
							// if (dose>=0 && (isCor==false || dose==1) ) {
								// score=hbdScores_lut[dose];
							// }else {
								// score=0.0;
							// }
//
							// pars->ibd->selfScores[i1][site_i]=score;
							// VALIDATION1
							// fprintf(stdout, "V1,%d,%d,%d,%d,%d,%f\n",dose,dose,i1,i2,site_i,score);

						}else{

							int dose1=vcfd->gts->get_n_derived_alleles_ind(i1);
							int dose2=vcfd->gts->get_n_derived_alleles_ind(i2);


							if(dose1==-1||dose2==-1){
								score=0.0;
							}else{


							// TODO HERE we dont have this in gl method
							bool notIbs = (dose1==0 && dose2==2) || (dose1==2 && dose2==0);
							if ( (notIbs || isCor==false) && dose1>=0 && dose2>=0 ) {
								score=ibdScores_lut[dosePairIndex[dose1][dose2]];
							}else {
								score=0.0;
							}

							int pair_idx=nCk_idx(nInd,i1,i2);
							pars->ibd->pairScores[pair_idx][site_i]=score;
							// VALIDATION1
							// VALIDATION GT1
							// fprintf(stdout, "VGT1,%d,%d,%d,%f\n",i1,i2,site_i,score);
							// fprintf(stdout, "VGT1,%d,%d,%d,%f\n",i1,i2,vcfd->rec->pos+1,score);
							}

						}


					}
				}

			} while (0);

		}else{
			NEVER;
		}
		/// -------------------

		site2pos[site_i]=vcfd->rec->pos+1;


		if (skip_site == 1) {
			IO::vprint(1, "Skipping site at %s:%ld", vcfd->get_contig_name(), vcfd->rec->pos + 1);
			pars->totSites++;
			site_i++;

			continue;
		}


		++site_i;
		pars->nSites++;
		pars->totSites++;
		nSites_contig_i++;


	}

	delete a1;
	delete a2;

	get_ibd_segments(vcfd->get_contig_name(), pars,site2pos);


	BEGIN_LOGSECTION;
	LOG("Finished reading sites.");
	LOG("Read %d (out of %d) contigs from the VCF file.", contig_i + 1, vcfd->nContigs);
	LOG("Number of contigs retained: %d", contig_i + 1);
	LOG("Number of contigs skipped: %d", vcfd->nContigs - contig_i - 1);
	LOG("Total number of sites processed: %lu", pars->totSites);
	LOG("Total number of sites skipped for all individual pairs: %lu", pars->totSites - pars->nSites);
	LOG("Total number of sites retained: %lu", pars->nSites);
	END_LOGSECTION;

	FREE(site2pos);
	if(NULL!=ibdScores_lut){
		FREE(ibdScores_lut);
	}


	if(NULL!=hbdScores_lut){
		FREE(hbdScores_lut);
	}

	if(NULL!=lngl){
		FREE(lngl);
	}

}

int get_ibd_segments(const char* contigName, paramStruct *pars, int* site2pos){
	
	double maxSum=0.0;
	int start=0;
	int end=0;

	double thisSum=0.0;
	double thisScore=0.0;
	int pair_idx=-1;
	const int nInd=pars->nInd;

	for (int i1 = 0; i1 < nInd; ++i1) {
		for (int i2 = i1+1; i2 < nInd; ++i2) {

			if(i1==i2){
				NEVER;
				pair_idx=-1;
			}else{
				pair_idx=nCk_idx(nInd,i1,i2);
			}

			thisSum=0.0;
			maxSum=0.0;
			start=0;
			end=0;

			//TODO work with multiple contigs

			for (int site=0;site<pars->nSites;site++){

				thisScore=0.0;

				if(pair_idx==-1){
					NEVER;
					// thisScore=pars->ibd->selfScores[i1][site];
				}else{
					thisScore=pars->ibd->pairScores[pair_idx][site];
				}

				thisSum+=thisScore;

				// VALIDATION
				// if(i1==0&&i2==1){
				// fprintf(stdout,"V4,%d,%d,%ld\n",i1,i2,site,thisScore,thisSum,maxSum);
					// fprintf(stdout,"%d,%d,%ld,%f\n",i1,i2,site,thisScore);
				// fprintf(stdout,"V3,%d,%d,%ld,%f,%f\n",i1,i2,site,maxSum,thisSum);
				// fprintf(stdout,"V4,%d,%d,%ld,%f,%f,%f\n",i1,i2,site,thisScore,thisSum,maxSum);
				// }

				if(thisSum > maxSum) {
					maxSum = thisSum;
					end = site;
				}else if (thisSum <= 0.0) {
					if (maxSum >= args->ibdseq_ibdlod) {
						start = trimmedStart(pars,i1, i2, start, end);
						end = trimmedEnd(pars,i1, i2, start, end);

						if (end > start) {

							if(-1==pair_idx){ //hbd
								NEVER;
								// fprintf(hbdOut,"%d\t%d\t%s\t%d\t%d\t%f\n",i1,i2,contigName,site2pos[start],site2pos[end], maxSum);
							}else{ //ibd
#if 1
								fprintf(stdout,"%d\t%d\t%s\t%d\t%d\t%f\n",i1,i2,contigName,site2pos[start],site2pos[end], maxSum);
#endif
							}
						}

					}
					start = site + 1;
					end = start;
					thisSum = 0.0;
					maxSum = 0.0;
				}


			}// sites loop

			if (maxSum >= args->ibdseq_ibdlod) {
				start = trimmedStart(pars,i1, i2, start, end);
				end = trimmedEnd(pars,i1, i2, start, end);
				if (end > start) {

					if(-1==pair_idx){ //hbd
									  // fprintf(hbdOut,"%d\t%d\t%s\t%d\t%d\t%f\n",i1,i2,contigName,site2pos[start],site2pos[end], maxSum);
						NEVER;
					}else{ //ibd
#if 1
						fprintf(stdout,"%d\t%d\t%s\t%d\t%d\t%f\n",i1,i2,contigName,site2pos[start],site2pos[end], maxSum);
#endif
					}
				}
			}


		}// i2 loop
	}// i1 loop

return 0;
}




int trimmedStart(paramStruct* pars, const int i1, const int i2, const int start, const int end) {
	if (args->ibdseq_ibdtrim <= 0.0) {
		return(start);
	}
	double sum = 0.0;
	int index = start;

	int pair_idx=-1;
	if(i1==i2){
		NEVER;
		pair_idx=-1;
	}else{
		pair_idx=nCk_idx(pars->nInd,i1,i2);
	}

	while (index <= end && sum < args->ibdseq_ibdtrim) {
		// sum += score(s1, s2, index, vcfData.isCorrelated(index));
		if(pair_idx==-1){
			NEVER;
			// sum+=pars->ibd->selfScores[i1][index];
		}else{
			sum+=pars->ibd->pairScores[pair_idx][index];
		}
		++index;
	}
	return(index - 1);
}

int trimmedEnd(paramStruct* pars, const int i1, const int i2, const int start, const int end) {
	if (args->ibdseq_ibdtrim <= 0.0) {
		return(end);
	}
	double sum = 0.0;
	int index = end;

	int pair_idx=-1;
	if(i1==i2){
		NEVER;
		pair_idx=-1;
	}else{
		pair_idx=nCk_idx(pars->nInd,i1,i2);
	}

	while (index >= start && sum < args->ibdseq_ibdtrim) {
		// sum += score(s1, s2, index, vcfData.isCorrelated(index));
		if(pair_idx==-1){
			NEVER;
			// sum+=pars->ibd->selfScores[i1][index];
		}else{
			sum+=pars->ibd->pairScores[pair_idx][index];
		}
		--index;
	}
	return(index + 1);
}


double* ibdStruct::errorArray(double e){
	double* errArray = (double*) malloc(5 * sizeof(double));
	double x = 1.0 - e;
	// errArray[j] = (errArray[j-1] * e) / x;// unroll
	errArray[0] = x*x*x*x;
	errArray[1] = (e) * (x*x*x);
	errArray[2] = (e*e) * (x*x);
	errArray[3] = (e*e*e) * (x);
	errArray[4] = (e*e*e*e);
	return (errArray);
}

double ibdStruct::errorRate(double fB){
	double e = args->ibdseq_errorprop * fB;
	return(MIN(e,args->ibdseq_errormax));
	// return(e); //TEST
}

double ibdStruct::ibdLike(int dose1, int dose2, double err, double pB){

	double ret=0.0;
	double pA = 1.0 - pB;
	
	double* e= NULL;

	if(err==args->ibdseq_errormax) {
		e=(double*) malloc (5 * sizeof (double));
		for(int i=0;i<5;++i){
			e[i]=this->maxErrorArray[i];
		}
	}else{
		e=this->errorArray(err);
	}

	switch (dose1 + dose2){
		case 0: 
			ret=(e[0]*(pA*pA*pA) + 2*e[1]*pA*pA*pB + e[2]*pA*pB);
			break;

		case 1: 
			ret=(2*(e[0]*pA*pA*pB + e[1]*(pA*pB + 2*(pA*pA*pA)) + 3*e[2]*pA*pB));
			break;

		case 2:
			if (dose1==dose2) {
				ret=((e[0]+4*e[1] + 2*e[2])*pA*pB + 4*e[2]*((pA*pA*pA) + (pB*pB*pB)));
			}
			else {
				ret=(2*( (e[1]+e[2]+e[3])*pA*pB + e[2]*((pA*pA*pA) + (pB*pB*pB)) ));
			}
			break;

		case 3: 
			ret=(2*(e[0]*pA*pB*pB + 2*e[1]*(pB*pB*pB) + (e[1] + 3*e[2] + e[3])*pA*pB + 2*e[3]*(pA*pA*pA) ));
			break;

		case 4: 
			ret=(( e[0]*(pB*pB*pB) + 2*e[1]*pA*pB*pB + e[2]*pA*pB + 2*e[3]*pA*pA*pB + e[4]*(pA*pA*pA) ));
			break;

		default: NEVER;

	}
	FREE(e);
	return(ret);

}

// IBD = 1 
double ibdStruct::ibdLike_GL(const int idx, vcfData* vcfd, paramStruct* pars, double fm){

	double fM=1.0-fm;
	//
	// double M[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	// M[0]=fM*fM*fM; // A
	// M[1]=fM*fM * fm; // B
	// // M[2]=0.0 // C
	// M[3]= M[1]; // D
	// M[4]=fM*fm; // E
	// M[5]=fM*fm*fm; // F
	// // M[6]=0.0; // G
	// M[7]=M[5]; // H
	// M[8]=fm*fm*fm; // I

	switch (idx){

		case 0:
			return(fM*fM*fM);

		case 1:
			return(fM*fM * fm);

		case 2:
			return  (0.0);

		case 3:
			return(fM*fM * fm);

		case 4:
			return (fM*fm);

		case 5:
			return (fM*fm*fm);

		case 6:
			return (0.0);

		case 7:
			return (fM*fm*fm);

		case 8:
			return (fm*fm*fm);


		default: NEVER;

	}



}


double ibdStruct::nullLike_GL(const int idx, vcfData* vcfd, paramStruct* pars,double fm){

	double fM=1.0-fm;

	// double M[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	// M[0]=fM*fM*fM; // A
	// M[1]=fM*fM * fm; // B
	// // M[2]=0.0 // C
	// M[3]= M[1]; // D
	// M[4]=fM*fm; // E
	// M[5]=fM*fm*fm; // F
	// // M[6]=0.0; // G
	// M[7]=M[5]; // H
	// M[8]=fm*fm*fm; // I

	switch (idx){

		case 0:
			return(fM*fM*fM*fM);

		case 1:
			return(2*fM*fM*fM*fm);

		case 2:
			return  (fM*fM*fm*fm);

		case 3:
			return  (2*fM*fM*fM*fm);

		case 4:
			return (4*fM*fM*fm*fm);

		case 5:
			return (2*fM*fm*fm*fm);

		case 6:
			return  (fM*fM*fm*fm);

		case 7:
			return (2*fM*fm*fm*fm);

		case 8:
			return (fm*fm*fm*fm);

		default: NEVER;

	}

}



double ibdStruct::ibdScore_GL(const int i1, const int i2, vcfData* vcfd, paramStruct* pars, const double fa, const size_t site_i, double* lngl){


	//TODO missing gls check

	double r=0.0;

	double i1gl1=lngl[(i1*vcfd->nGT)+0];
	double i1gl2=lngl[(i1*vcfd->nGT)+1];
	double i1gl3=lngl[(i1*vcfd->nGT)+2];
	if((i1gl1==i1gl2)&&(i1gl1==i1gl3)){//TODO delme
		NEVER;
	}

	double i2gl1=lngl[(i2*vcfd->nGT)+0];
	double i2gl2=lngl[(i2*vcfd->nGT)+1];
	double i2gl3=lngl[(i2*vcfd->nGT)+2];
	if((i2gl1==i2gl2)&&(i2gl1==i2gl3)){//TODO delme
		NEVER;
	}

	double GLS[9] = {0.0};
	GLS[0]=exp(i1gl1+i2gl1);
	GLS[1]=exp(i1gl2+i2gl1);
	GLS[2]=exp(i1gl3+i2gl1);
	GLS[3]=exp(i1gl1+i2gl2);
	GLS[4]=exp(i1gl2+i2gl2);
	GLS[5]=exp(i1gl3+i2gl2);
	GLS[6]=exp(i1gl1+i2gl3);
	GLS[7]=exp(i1gl2+i2gl3);
	GLS[8]=exp(i1gl3+i2gl3);


	double m=0.0;
	double n=0.0;
	for(int i=0;i<9;++i){
		if(GLS[i]!=0.0){

// fprintf(stderr,"%ld i1:%d i2:%d i1gl1:%f i1gl2:%f i1gl3:%f i2gl1:%f i2gl2:%f i2gl3:%f\n",site_i,i1,i2,i1gl1,i1gl2,i1gl3,i2gl1,i2gl1,i1gl3);

			m+=GLS[i]*ibdLike_GL(i,vcfd,pars,fa);
			n+=GLS[i]*nullLike_GL(i,vcfd,pars,fa);
		}
	}



	ASSERT(0.0!=n);
	ASSERT(0.0!=m);

	r=(double) (m/n);

	ASSERT(r>0.0);

	r=(double) log10(r);
	return r;
}

double ibdStruct::nullLike(int dose1, int dose2, double fB){

	double fA = 1.0 - fB;

	switch (dose1 + dose2) {
		case 0: 
			return (fA*fA*fA*fA);

		case 1: 
			return (4*(fA*fA*fA) * fB);

		case 2:
			if (dose1==dose2) {
				return ((2*fA*fB)*(2*fA*fB));
			}
			else {
				return ((fA*fB)*(fA*fB))+((fA*fB)*(fA*fB));
			}

		case 3: 
			return (4*fA*(fB*fB*fB));

		case 4: 
			return (fB*fB*fB*fB);

		default: 
			NEVER;
	}

}

double estimate_true_maf(double fB, double errorRate){
	// return (fB);//TEST
	return (( fB - errorRate ) / (1.0 - 2 * errorRate));
}


double ibdStruct::ibdScore(const int dose1, const int dose2, const double fB){
	if(dose1<0 || dose2<0){
		NEVER;
		// unknown dose
		return(0.0);
	}
	ASSERT(fB<=0.5);

	double e = this->errorRate(fB);
	double pB = estimate_true_maf(fB, e);
	double r = ibdLike(dose1, dose2, e, pB) / nullLike(dose1, dose2, fB);
	return log10(r);
}


double ibdStruct::hbdScore(int dose, double fB){
	ASSERT(fB > 0.0 && fB < 1.0);
	ASSERT(dose<3 && dose>-1);
	if(dose<0){
		NEVER;
		// return 0.0;
	}

	double e = this->errorRate(fB);
	double fA = 1.0 - fB;
	double pB = estimate_true_maf(fB, e);
	double pA = 1.0 - pB;
	switch (dose) {
		case 0: return log10((pA + e*e*pB)/(fA*fA));
		case 1: return log10( e*(1-e)/(fA*fB) );
		case 2: return log10( (e*e*pA + pB)/(fB*fB) );
		default: NEVER;
	}

}


// assuming  N_DOSES=3
// index[0][0]=0;
// index[0][1]=1;
// index[0][2]=5;
// index[1][0]=1;
// index[1][1]=2;
// index[1][2]=3;
// index[2][0]=5;
// index[2][1]=3;
// index[2][2]=4;





// void ibdStruct::add_doseData(const int site_i,const int ind_i, const int dose){
//
// this->doseData[site_i][ind_i]=dose;
// }
//

ibdStruct::ibdStruct(vcfData* vcfd, paramStruct* pars){


	this->nInd=vcfd->nInd;
	if(args->doIbd==2){
		this->maxErrorArray = errorArray(args->ibdseq_errormax);
	}

	// ----------
	// // nIndices=6 from ibdseq
	// ibdScores[nIndices][nMarkers]
	this->ibdScores = (double*) malloc(6 * sizeof(double));

	this->hbdScores = (double*) malloc(3 * sizeof(double));


	this->pairScores = (double**) malloc(pars->nIndCmb * sizeof(double*));
	for (size_t i=0;i<pars->nIndCmb;++i){
		this->pairScores[i]=(double*)malloc(BUF_NSITES*sizeof(double));
	}

	// this->selfScores = (double**) malloc(pars->nInd * sizeof(double*));
	// for (size_t i=0;i<pars->nInd;++i){
		// this->selfScores[i]=(double*)malloc(BUF_NSITES*sizeof(double));
	// }

}


ibdStruct::~ibdStruct(){
	FREE(ibdScores);
	FREE(hbdScores);
	if(NULL!=maxErrorArray){
		FREE(maxErrorArray);
	}

}



