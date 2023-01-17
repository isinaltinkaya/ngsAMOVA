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

double calculate_SumOfSquares(int lvl_a, DATA::distanceMatrixStruct *dMS, DATA::metadataStruct *mS, DATA::samplesStruct *sS, FILE *out_amova_ff, int **LUT_indPair_idx, const char *analysis_type)
{

	int ss_lvl = lvl_a - 1;

	// if at the highest level, then calculate the total sum of squares
	if(ss_lvl == mS->nLevels){
		return calculate_SumOfSquares_Total(dMS);
	}else{

		double sum = 0.0;
		double ssd = 0.0;

		// loop over different stratas at this hierarchical level
		// e.g. if hierarchical level is region, then loop over all regions {region1, region2, region3 ...}
		for (int sti=0; sti < mS->hierArr[ss_lvl]->nStrata; sti++)
		{

			// sum of squared distances within this strata
			sum = 0.0;
			for (int i1 = 0; i1 < dMS->nInd - 1; i1++)
			{
				for (int i2 = i1 + 1; i2 < dMS->nInd; i2++)
				{
					// include only pairs that are in the same strata at this hierarchical level
					if(mS->pairInStrataAtLevel(i1, i2, ss_lvl, sti) == 1){
						
						// the distance is already stored in squared form, so no need to square it again
						sum += dMS->M[LUT_indPair_idx[i1][i2]];
					}
				}
			}
			fprintf(stderr,"\n\n\t->nIndPerStrata[%i] = %i", sti, mS->hierArr[ss_lvl]->nIndPerStrata[sti]);
			ssd += sum / (double)mS->hierArr[ss_lvl]->nIndPerStrata[sti];
		}
		return ssd;

	}

}


int AMOVA::doAMOVA(DATA::distanceMatrixStruct *dMS, DATA::metadataStruct *mS, DATA::samplesStruct *sS, FILE *out_amova_ff, int **LUT_indPair_idx, const char *analysis_type)
{

	amovaResultStruct *aRS = new amovaResultStruct(mS);



	// ----------------------------------------------- //
	// DEGREES OF FREEDOM 
	// ----------------------------------------------- //
	// calculate the degrees of freedom for each level
	
	if(mS->nLevels == 1){
		// highest level df is nStrata - 1
		// e.g. Individual ~ Region / Population / Subpopulation
		// if there are 3 regions, then df = 3 - 1 = 2
		aRS->df[0] = mS->hierArr[0]->nStrata - 1;

		
		// df for lowest amova level
		// e.g. Individual ~ Region / Population / Subpopulation
		// nInd - sum(for each Population, nStrata (num_Subpopulation) in Population)

		aRS->df[mS->nLevels] = dMS->nInd - mS->hierArr[0]->nStrata;

		// last element in df array is total df
		aRS->df[aRS->_ssd-1] = dMS->nInd - 1;

	}else if(mS->nLevels==2){
		fprintf(stderr, "\ndf[0] = %i, df[1] = %i, df[2] = %i, df[3] = %i", aRS->df[0], aRS->df[1], aRS->df[2], aRS->df[3]);
		// last element: total df
		aRS->df[aRS->_ssd-1] = dMS->nInd - 1;

		aRS->df[0] = mS->sumUniqStrataAtEachLevel(1) - 1;

		aRS->df[1] = mS->sumUniqStrataAtEachLevel(0) - mS->sumUniqStrataAtEachLevel(1);
		
		aRS->df[mS->nLevels] = dMS->nInd - mS->sumUniqStrataAtEachLevel(0);

	}else{
		ASSERT(0==1);
	}

	


	int lvl_a = aRS->nAmovaLevels;
	while(lvl_a > 0){
		// within level sum of squares
		aRS->ss[lvl_a] = calculate_SumOfSquares(lvl_a, dMS, mS, sS, out_amova_ff, LUT_indPair_idx, analysis_type);
		lvl_a--;
	}

	if( mS->nLevels == 1){
		aRS->ssd[0] = aRS->ss[2] - aRS->ss[1];
		aRS->ssd[1] = aRS->ss[1];
		aRS->ssd[2] = aRS->ss[2];
	}else{
		
		aRS->ssd[0] = aRS->ss[aRS->_ssd-1] - aRS->ss[1];
		// aRS->ssd[1] = aRS->ss[1] - aRS->ss[0];
		// aRS->ssd[2] = aRS->ss[2] - aRS->ss[1];
		aRS->ssd[aRS->nAmovaLevels] = aRS->ss[aRS->nAmovaLevels];
	}



	// get msd from ssd and df
	for(size_t i=0; i < aRS->_ssd; i++){
		aRS->msd[i] = aRS->ssd[i] / aRS->df[i];
	}



	// variance coefficients



	if(mS->nLevels == 1){
		double n_gi = 0.0;

		// n = [ N - \sum_{g \in G} ( N^2_{g}/N) ) ]  /   G - 1
		for (int sti=0; sti < mS->hierArr[0]->nStrata; sti++)
		{
			n_gi += SQUARE(mS->hierArr[0]->nIndPerStrata[sti]) / (double) dMS->nInd;
		}
		aRS->ncoef[0] = (double)((double)dMS->nInd - (double)n_gi) / (double)(mS->hierArr[0]->nStrata - 1);

		aRS->sigmasq[0] = (aRS->msd[0]-aRS->msd[1]) / aRS->ncoef[0];
		aRS->sigmasq[1] = aRS->msd[1];


		aRS->phi[0] = aRS->sigmasq[0] / (aRS->sigmasq[0] + aRS->sigmasq[1]);
	}else if(mS->nLevels == 2){
		// eq.9a-c in Excoffier1992
		// ncoef = {n, n', n''}
		
		
		// estimate the method of moments: n, n', and n''

		// // n = [ N - \sum_{g \in G} ( N^2_{g}/N) ) ]  /   G - 1
		// double n_gi = 0.0;
		// for (int sti=0; sti < mS->hierArr[0]->nStrata; sti++)
		// {
		// 	n_gi += SQUARE(mS->hierArr[0]->nIndPerStrata[sti]) / (double) dMS->nInd;
		// }
		// aRS->ncoef[0] = (double)((double)dMS->nInd - (double)n_gi) / (double)(mS->hierArr[0]->nStrata - 1);

		// double n_moment1 = 0.0;
		// double n_moment2 = 0.0;
		// double n_moment3 = 0.0;
		// for (int sti=0; sti < mS->hierArr[1]->nStrata; sti++)
		// {
		// 	n_moment1 += mS->hierArr[1]->nIndPerStrata[sti];
		// 	n_moment2 += SQUARE(mS->hierArr[1]->nIndPerStrata[sti]);
		// }
		// aRS->ncoef[1] = (double)((double)n_moment1 - (double)n_moment2 / (double)n_moment1) / (double)(mS->hierArr[1]->nStrata - 1);

		// n'' = [ N - \sum_{g \in G} ( N^2_{g}/N) ) ]  /   G - 1


		// // n'' = [ N - \sum_{g \in G} ( N^2_{g}/N) ) ]  /   G - 1
		// double n_gi3 = 0.0;
		// for (int sti=0; sti < mS->hierArr[2]->nStrata; sti++)
		// {
		// 	n_gi3 += SQUARE(mS->hierArr[2]->nIndPerStrata[sti]) / (double) dMS->nInd;
		// }
		// aRS->ncoef[2] = (double)((double)dMS->nInd - (double)n_gi3) / (double)(mS->hierArr[2]->nStrata - 1);




	}else{
		// use approximation

	}

	aRS->print_as_table(stderr,mS);
	// aRS->print_as_csv(out_amova_ff, analysis_type);

	delete aRS;
	return 0;
}
