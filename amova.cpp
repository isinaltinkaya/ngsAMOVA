#include "amova.h"

double calculate_SumOfSquares_Total(DATA::distanceMatrixStruct *dMS)
{

	double ssd_TOTAL = 0.0;

	for (int px = 0; px < dMS->nIndCmb; px++)
	{
		// the distance is already stored in squared form
		ssd_TOTAL += dMS->M[px];
	}

	ssd_TOTAL = ssd_TOTAL / (double)dMS->nInd;
	return ssd_TOTAL;
}

/// @brief calculate_SumOfSquares_Within - calculate the sum of squares within a hierarchical level
/// @param lvl hierarchical level
/// @param aS  pointer to amovaStruct
/// @param dMS pointer to distanceMatrixStruct
/// @param mS  pointer to metadataStruct
/// @param sS  pointer to sampleStruct
/// @param out_amova_ff    pointer to output file
/// @param LUT_inds2idx  lookup table for index pair to index in distance matrix
/// @return (double) sum of squares within a hierarchical level
double calculate_SumOfSquares_Within(int lvl, AMOVA::amovaStruct *aS, DATA::distanceMatrixStruct *dMS, DATA::metadataStruct *mS, DATA::sampleStruct *sS, FILE *out_amova_ff, int **LUT_inds2idx)
{

	ASSERT(lvl >= 0 && lvl < aS->nAmovaLevels);
	// if at the highest level, then calculate the total sum of squares
	if(lvl == aS->nAmovaLevels-2){
		return calculate_SumOfSquares_Total(dMS);
	}else{

		double sum = 0.0;
		double ssd = 0.0;
		int n=0;
		// loop over different stratas at this hierarchical level
		// e.g. if hierarchical level is region, then loop over all regions {region1, region2, region3 ...}
		for (int sti=0; sti < mS->hierArr[lvl]->nStrata; sti++)
		{

			// sum of squared distances within this strata
			sum = 0.0;
			n=mS->hierArr[lvl]->nIndPerStrata[sti];
			for (int i1 = 0; i1 < dMS->nInd - 1; i1++)
			{
				for (int i2 = i1 + 1; i2 < dMS->nInd; i2++)
				{
					// include only pairs that are in the same strata at this hierarchical level
					if(mS->pairInStrataAtLevel(i1, i2, lvl, sti) == 1){
						
						// the distance is already stored in squared form, so no need to square it again
						sum += dMS->M[LUT_inds2idx[i1][i2]]/n;
					}
				}
			}
			ssd += sum;
		}
		return ssd;
	}
}


/// @brief doAMOVA - perform AMOVA analysis
/// @param dMS  pointer to distanceMatrixStruct
/// @param mS   pointer to metadataStruct
/// @param sS   pointer to sampleStruct
/// @param out_amova_ff   pointer to output file
/// @param LUT_inds2idx lookup table for index pair to index in distance matrix
/// @return 
int AMOVA::doAMOVA(DATA::distanceMatrixStruct *dMS, DATA::metadataStruct *mS, DATA::sampleStruct *sS, FILE *out_amova_ff, int **LUT_inds2idx)
{

	// ----------------------------------------------- //

	// if(mS->nLevels == 1){
	// }else if(mS->nLevels==2){
	// }else{
	// 	fprintf(stderr, "[ERROR] NOT IMPLEMENTED YET");
	// 	ASSERT(0 == 1);
	// }
	amovaStruct *aS = new amovaStruct(mS);


	// ----------------------------------------------- //
	// DEGREES OF FREEDOM 
	// ----------------------------------------------- //
	// calculate the degrees of freedom for each level
	

	// highest level df is nStrata - 1
	// e.g. Individual ~ Region / Population / Subpopulation
	// if there are 3 regions, then nStrata is number of unique regions -1 = 3-1 = 2
	aS->df[0] = mS->hierArr[0]->nStrata - 1;
	// aS->df[0] = mS->sumUniqStrataAtEachLevel(1) - 1;

	// total df
	aS->df[aS->nAmovaLevels-1] = dMS->nInd - 1;

	if(mS->nLevels == 1){

		// df for lowest amova level
		// e.g. Individual ~ Region / Population / Subpopulation
		// nInd - sum(for each Population, nStrata (num_Subpopulation) in Population)
		aS->df[1] = dMS->nInd - mS->hierArr[0]->nStrata;

	}else if(mS->nLevels==2){

		aS->df[1] = mS->sumUniqStrataAtEachLevel(0) - mS->sumUniqStrataAtEachLevel(1);
		aS->df[2] = dMS->nInd - mS->sumUniqStrataAtEachLevel(0);

	}else{
		fprintf(stderr, "[ERROR] NOT IMPLEMENTED YET");
		ASSERT(0 == 1);
	}

	

	// ----------------------------------------------- //
	// SUM OF SQUARES
	// ----------------------------------------------- //
	// calculate the sum of squares within for each level
	//

	for(size_t i=0; i<aS->_ncoef; i++){
		// within level sum of squares
		aS->ss[i] = calculate_SumOfSquares_Within(i, aS, dMS, mS, sS, out_amova_ff, LUT_inds2idx);
	}

	// ----------------------------------------------- //
	// SUM OF SQUARED DEVIATIONS
	// ----------------------------------------------- //

	if( mS->nLevels == 1){
	
		aS->ssd[0] = aS->ss[1] - aS->ss[0];
		aS->ssd[1] = aS->ss[0];
		aS->ssd[2] = aS->ss[1];

	}else if(mS->nLevels == 2){


		aS->ssd[3] = aS->ss[2];
		aS->ssd[1] = aS->ss[0] - aS->ss[1];
		aS->ssd[2] = aS->ss[1];
		aS->ssd[0] = aS->ssd[3] - (aS->ssd[0]+aS->ssd[1]+aS->ssd[2]);

	}else{
		fprintf(stderr, "[ERROR] NOT IMPLEMENTED YET");
		ASSERT(0 == 1);
	}

	// get msd from ssd and df
	for(size_t i=0; i < (size_t) aS->nAmovaLevels; i++){
		aS->msd[i] = aS->ssd[i] / aS->df[i];
	}

	if(mS->nLevels == 1){

		// ----------------------------------------------- //
		// VARIANCE COEFFICIENTS
		// ----------------------------------------------- //
		// calculate variance coefficients (n) for 1 level AMOVA

		double n_gi = 0.0;

		// n = [ N - \sum_{g \in G} ( N^2_{g}/N) ) ]  /   G - 1
		for (int sti=0; sti < mS->hierArr[0]->nStrata; sti++)
		{
			n_gi += SQUARE(mS->hierArr[0]->nIndPerStrata[sti]) ;
		}
		n_gi = n_gi / (double) dMS->nInd;
		aS->ncoef[0] = (double)((double)dMS->nInd - (double)n_gi) / (double)(mS->hierArr[0]->nStrata - 1);


		// ----------------------------------------------- //
		// VARIANCE COMPONENTS
		// ----------------------------------------------- //
		// calculate variance components (sigma squared)
		// for 1 level AMOVA

		aS->sigmasq[0] = (aS->msd[0]-aS->msd[1]) / aS->ncoef[0];
		aS->sigmasq[1] = aS->msd[1];

		aS->phi[0] = aS->sigmasq[0] / (aS->sigmasq[0] + aS->sigmasq[1]);

	}else if(mS->nLevels == 2){


		// ----------------------------------------------- //
		// VARIANCE COEFFICIENTS
		// ----------------------------------------------- //
		// calculate variance coefficients for 2 level AMOVA
		//
		// eq.9a-c in Excoffier1992
		// estimate the method of moments: n, n', and n''
		// ncoef = {n, n', n''}

		double n_moment1 = 0.0;
		double n_moment2 = 0.0;
		for (int sti=0; sti < mS->hierArr[1]->nStrata; sti++)
		{
			n_moment1 += mS->hierArr[1]->nIndPerStrata[sti];
			n_moment2 += SQUARE(mS->hierArr[1]->nIndPerStrata[sti]);
		}
		aS->ncoef[1] = (double)((double)n_moment1 - (double)n_moment2 / (double)n_moment1) / (double)(mS->hierArr[1]->nStrata - 1);
		double x=((SQUARE(mS->hierArr[0]->nIndPerStrata[0]))/dMS->nInd) +((SQUARE(mS->hierArr[0]->nIndPerStrata[1]))/dMS->nInd);

		double n_gi = 0.0;
		for (int sti=0; sti < mS->hierArr[0]->nStrata; sti++)
		{
			n_gi += SQUARE(mS->hierArr[0]->nIndPerStrata[sti]) ;
		}
		n_gi = n_gi / (double) dMS->nInd;

		double sum=0.0;
		for(int i=0; i< mS->hierArr[0]->nStrata; i++){
			sum+=mS->nMemberStrataAtStrataAtLevel(i,0)-1;
		}
		aS->ncoef[0] = (dMS->nInd - (n_moment1 - (n_moment2/n_moment1))) / sum;

		aS->ncoef[2] = (double) (dMS->nInd -  x) / (double) (mS->hierArr[0]->nStrata - 1);

		// ----------------------------------------------- //
		// VARIANCE COMPONENTS
		// ----------------------------------------------- //
		// calculate variance components (sigma squared)
		aS->sigmasq[2] = aS->msd[2];
		aS->sigmasq[1] = (aS->msd[1]-aS->sigmasq[2]) / aS->ncoef[0];
		aS->sigmasq[0] = (aS->msd[0] - aS->msd[2] - (aS->ncoef[1] * aS->sigmasq[1])) / aS->ncoef[2];


		// ----------------------------------------------- //
		// PHI STATISTIC
		// ----------------------------------------------- //
		// calculate phi statistics
		//
		// Individual,Region,Population,Total
		// 
		

		//lvl0 (i.e. reg) in TOTAL
		//Phi_CT
		aS->phi[0] = aS->sigmasq[0] / (aS->sigmasq[0] + aS->sigmasq[1] + aS->sigmasq[2]);

		//lvl1 in lvl0
		//Phi_SC
		aS->phi[1] = aS->sigmasq[1] / (aS->sigmasq[1] + aS->sigmasq[2]);

		//lvl1 (i.e. pop) in TOTAL
		//Phi_ST
		aS->phi[2] = (aS->sigmasq[0] + aS->sigmasq[1])/ (aS->sigmasq[0] + aS->sigmasq[1] + aS->sigmasq[2]);


	}else{
		// use approximation
	}

	aS->print_as_table(stderr,mS);
	aS->print_as_csv(out_amova_ff,mS);

	delete aS;
	return 0;
}

