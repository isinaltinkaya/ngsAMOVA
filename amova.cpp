#include "amova.h"

double get_SSD_total(DATA::distanceMatrixStruct *dMS)
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

double get_SSD_at_level(int amova_lvl_i, DATA::distanceMatrixStruct *dMS, DATA::metadataStruct *mS, DATA::samplesStruct *sS, FILE *out_amova_ff, int **LUT_indPair_idx, const char *analysis_type)
{
	int lvl_i = amova_lvl_i-1;

	double n = 0.0;
	double sum = 0.0;
	double ssd = 0.0;
	int px;

	if(lvl_i >= mS->nLevels){
		fprintf(stderr, "[ERROR] lvl_i (%i) > mS->nLevels (%i)\n", lvl_i, mS->nLevels);
		exit(1);
	}

	for (int sti=0; sti < mS->hierArr[lvl_i]->nStrata; sti++)
	{

		sum = 0.0;
		for (int i1 = 0; i1 < dMS->nInd - 1; i1++)
		{
			for (int i2 = i1 + 1; i2 < dMS->nInd; i2++)
			{
				if(mS->pair_in_strata(sti, i1, i2, sS, lvl_i) == 1){
					
					px = LUT_indPair_idx[i1][i2];

					// the distance is already stored in squared form
					sum += dMS->M[px];
				}
			}
		}
		n = (double)mS->hierArr[lvl_i]->nIndPerStrata[sti];
		ssd += sum / n;
	}
	return ssd;
}


int AMOVA::doAMOVA(DATA::distanceMatrixStruct *dMS, DATA::metadataStruct *mS, DATA::samplesStruct *sS, FILE *out_amova_ff, int **LUT_indPair_idx, const char *analysis_type)
{


// give mS directly
// ars_get func 
// also get names of levels to auto use in the table

	amovaResultStruct *aRS = new amovaResultStruct(mS);

	// ssd=[0=AG,1=AIWG,2=TOTAL]


	// fprintf(stderr, "\n\nssd_TOTAL = %f\n\n", aRS->ssd[aRS->nAmovaLevels]);

	// fprintf(stderr,"\n\nAmong %s within %s\n\n", mS->levelNames[0], mS->levelNames[1]);



	// calculate the degrees of freedom for each level
	aRS->df[0] = mS->hierArr[mS->nLevels-1]->nStrata - 1;
	for (size_t i=1; i < aRS->_ssd-1; i++){
		for (int j = 0; j < mS->hierArr[i-1]->nStrata; j++){
			aRS->df[i] += mS->hierArr[i-1]->nIndPerStrata[j]-1;
		}
	}
	// last element in df array is total df
	aRS->df[aRS->_ssd-1] = dMS->nInd - 1;



	// ssd_TOTAL 
	aRS->ssd[aRS->nAmovaLevels] = get_SSD_total(dMS);

	int lvl_i = aRS->nAmovaLevels-1;
	while(lvl_i > 0){
		// fprintf(stderr,"\n\n\nlvl_i %i\n\n\n",lvl_i);
		aRS->ssd[lvl_i] = get_SSD_at_level(1, dMS, mS, sS, out_amova_ff, LUT_indPair_idx, analysis_type);
		lvl_i--;
	}

	aRS->ssd[0] = aRS->ssd[aRS->nAmovaLevels] - aRS->ssd[aRS->nAmovaLevels-1];


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


	// variance components
	// for (size_t i=0; i < aRS->nAmovaLevels-1; i++){
	// }


	aRS->print_as_table(stderr,mS);
	aRS->print_as_csv(out_amova_ff, analysis_type);

	delete aRS;
	return 0;
}
