#include "amova.h"

double calculate_SumOfSquares_Total(distanceMatrixStruct *dMS)
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
/// @param pars pointer to paramstruct 
/// @return (double) sum of squares within a hierarchical level
double calculate_SumOfSquares_Within(int lvl, AMOVA::amovaStruct *aS, distanceMatrixStruct *dMS, metadataStruct *mS, sampleStruct *sS, paramStruct *pars)
{

	ASSERT(lvl >= 0 && lvl < aS->nAmovaLevels);
	// if at the highest level, then calculate the total sum of squares
	if (lvl == aS->nAmovaLevels - 2)
	{
		return calculate_SumOfSquares_Total(dMS);
	}
	else
	{
		double sum = 0.0;
		double ssd = 0.0;
		int n = 0;

		for (int g=0; g < mS->hierArr[lvl]->nStrata; g++)
		{
			sum=0.0;
			n = mS->hierArr[lvl]->nIndPerStrata[g];

			for(int i1=0; i1 < dMS->nInd-1; i1++){
				for(int i2=i1+1; i2 < dMS->nInd; i2++){
					//TODO do this by bitshifting
					if(mS->pairToAssoc[pars->lut_indsToIdx[i1][i2]][lvl][g] == 1){
						sum += dMS->M[pars->lut_indsToIdx[i1][i2]] / n;
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
/// @param pars pointer to paramStruct
/// @return
AMOVA::amovaStruct *AMOVA::doAmova(distanceMatrixStruct *dMS, metadataStruct *mS, sampleStruct *sS, paramStruct *pars)
{

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
	aS->df[aS->nAmovaLevels - 1] = dMS->nInd - 1;

	if (mS->nLevels == 1)
	{

		// df for lowest amova level
		// e.g. Individual ~ Region / Population / Subpopulation
		// nInd - sum(for each Population, nStrata (num_Subpopulation) in Population)
		aS->df[1] = dMS->nInd - mS->hierArr[0]->nStrata;
	}
	else if (mS->nLevels == 2)
	{

		aS->df[1] = mS->sumUniqStrataAtEachLevel(0) - mS->sumUniqStrataAtEachLevel(1);
		aS->df[2] = dMS->nInd - mS->sumUniqStrataAtEachLevel(0);
	}
	else
	{
		fprintf(stderr, "[ERROR] NOT IMPLEMENTED YET");
		ASSERT(0 == 1);
	}

	// ----------------------------------------------- //
	// SUM OF SQUARES
	// ----------------------------------------------- //
	// calculate the sum of squares within for each level
	//

	for (size_t i = 0; i < aS->_ncoef; i++)
	{
		// within level sum of squares
		aS->ss[i] = calculate_SumOfSquares_Within(i, aS, dMS, mS, sS, pars);
	}

	// ----------------------------------------------- //
	// SUM OF SQUARED DEVIATIONS
	// ----------------------------------------------- //

	if (mS->nLevels == 1)
	{

		aS->ssd[0] = aS->ss[1] - aS->ss[0];
		aS->ssd[1] = aS->ss[0];
		aS->ssd[2] = aS->ss[1];
	}
	else if (mS->nLevels == 2)
	{

		aS->ssd[3] = aS->ss[2];
		aS->ssd[1] = aS->ss[0] - aS->ss[1];
		aS->ssd[2] = aS->ss[1];
		aS->ssd[0] = aS->ssd[3] - (aS->ssd[0] + aS->ssd[1] + aS->ssd[2]);
	}
	else
	{
		fprintf(stderr, "[ERROR] NOT IMPLEMENTED YET");
		ASSERT(0 == 1);
	}

	// get msd from ssd and df
	for (size_t i = 0; i < (size_t)aS->nAmovaLevels; i++)
	{
		aS->msd[i] = aS->ssd[i] / aS->df[i];
	}

	if (mS->nLevels == 1)
	{

		// ----------------------------------------------- //
		// VARIANCE COEFFICIENTS
		// ----------------------------------------------- //
		// calculate variance coefficients (n) for 1 level AMOVA

		double n_gi = 0.0;

		// n = [ N - \sum_{g \in G} ( N^2_{g}/N) ) ]  /   G - 1
		for (int sti = 0; sti < mS->hierArr[0]->nStrata; sti++)
		{
			n_gi += SQUARE(mS->hierArr[0]->nIndPerStrata[sti]);
		}
		n_gi = n_gi / (double)dMS->nInd;
		aS->ncoef[0] = (double)((double)dMS->nInd - (double)n_gi) / (double)(mS->hierArr[0]->nStrata - 1);

		// ----------------------------------------------- //
		// VARIANCE COMPONENTS
		// ----------------------------------------------- //
		// calculate variance components (sigma squared)
		// for 1 level AMOVA
		aS->sigmasq_total = aS->sigmasq[0] + aS->sigmasq[1];

		aS->sigmasq[0] = (aS->msd[0] - aS->msd[1]) / aS->ncoef[0];
		aS->sigmasq[1] = aS->msd[1];

		aS->phi[0] = aS->sigmasq[0] / (aS->sigmasq[0] + aS->sigmasq[1]);
	}
	else if (mS->nLevels == 2)
	{

		// ----------------------------------------------- //
		// VARIANCE COEFFICIENTS
		// ----------------------------------------------- //
		// calculate variance coefficients for 2 level AMOVA
		//
		// eq.9a-c in Excoffier1992
		// estimate the method of moments: n, n', and n''
		// ncoef[2] = {n, n', n''} == {n, n1, n2}

		// n =
		// 		(NUL) \frac{ \sum^G_{g=1} \sum^{I_g}_{i=1} N_{ig} - 
		// 		(NUR)		\sum^G_{g=1} ( \frac{ \sum^{I_g}_{i=1} N^2_{ig}}{\sum^{I_g}_{i=1} N_{ig}} ) }
		// 		(NL)	{ \sum^G_{g=1} (I_g - 1) }
		//
		// n = (NUL-NUR) / NL
		// 		NUL = n formula upper left
		// 		NUR = n formula upper right
		// 		NL = n formula lower
		// 
		// NUL=0;
		// for each region
		// 		for each population in region
		// 			NUL += number of individuals in population (nIndPerStrata)
		//
		// loop through the groups in the highest hierarchy level hierArr[0] (e.g. regions)
		double NUL=0.0;
		for (int g=0; g < mS->hierArr[0]->nStrata; g++)
		{
			for (int sg=0; sg < mS->hierArr[0]->nSubStrata[g]; sg++){
				int sg_idx = mS->hierArr[0]->subStrataIdx[g][sg];
				NUL += mS->hierArr[1]->nIndPerStrata[sg_idx];
			}
		}


		//
		// NUR=0;
		// for each region
		// 		Nig2_g=0;
		//		Nig_g=0;
		// 		for each population in region
		// 			Nig2_g += nIndPerStrata^2
		// 		for each population in region
		// 			Nig_g += nIndPerStrata
		// 		NUR += (Nig2_g / Nig_g)
		//
		double NUR=0.0;
		double Nig2_g=0.0;
		double Nig_g=0.0;

		for (int g=0; g < mS->hierArr[0]->nStrata; g++)
		{
			Nig2_g=0.0;
			Nig_g=0.0;
			for (int sg=0; sg < mS->hierArr[0]->nSubStrata[g]; sg++){
				int sg_idx = mS->hierArr[0]->subStrataIdx[g][sg];
				Nig2_g += SQUARE(mS->hierArr[1]->nIndPerStrata[sg_idx]);
				Nig_g += mS->hierArr[1]->nIndPerStrata[sg_idx];
			}
			NUR += Nig2_g / Nig_g;
		}

		// 
		// NL=0;
		// for each region
		// 		NL += (number of populations in region - 1 )
		//
		// 
		double NL=0.0;
		for (int g=0; g < mS->hierArr[0]->nStrata; g++)
		{
			NL += mS->hierArr[0]->nSubStrata[g]-1;
		}

		aS->ncoef[0] = (NUL-NUR) / NL;


		// n' =
		//  (NUR) \frac{ \sum^G_{g=1} (\frac{\sum^{I_g}_{j=1} N^2_{ig} }{\sum^{I_g}_{i=1} N{ig}  } ) -
		//  (N1UR) 	\frac{\sum^G_{g=1} \sum^{I_g}_{j=1} N^2_{ig} }
		//  (N1L)	{ G - 1 }
		//
		// n' = (NUR-N1UR)/N1L
		//
		// NUR => already calculated above
		// 
		// N1UR =0;
		// Nig2=0;
		// Nig=0;
		// for each region
		// 		for each population in region
		// 			Nig2 += nIndPerStrata^2
		//
		// for each region
		// 		for each population in region
		// 			Nig += nIndPerStrata
		//
		// N1UR = Nig2 / Nig
		// 				Nig => already calculated above as NUL
		//
		double N1UR=0.0;
		double Nig2=0.0;

		// g=group, sg=subgroup
		for (int g=0; g < mS->hierArr[0]->nStrata; g++)
		{
			for (int sg=0; sg < mS->hierArr[0]->nSubStrata[g]; sg++){

				int sg_idx = mS->hierArr[0]->subStrataIdx[g][sg];
				int n_sgi = mS->hierArr[1]->nIndPerStrata[sg_idx];
				Nig2 += SQUARE(n_sgi);
			}
		}

		N1UR = Nig2 / NUL;

		// N1L=0;
		// N1L = nRegions - 1;
		//
		double N1L=0.0;
		N1L = mS->hierArr[0]->nStrata - 1;

		aS->ncoef[1] = (NUR-N1UR)/N1L;

		// n'' =
		// 	(NUL) frac{\sum^G_{g=1} \sum^{I_g}_{i=1} N_{ig} - 
		//	(N2UR) 	\frac{\sum^G_{g=1} \left( \sum^{I_g}_{i=1} N_{ig} \right)^2 }
		// 	(N1L) {\sum^G_{g=1} \sum^{I_g}_{i=1} N_{ig} }}{ G - 1 } 
		//
		// n'' = (NUL-N2UR)/N1L
		//
		// NUL, N1L => already calculated above 
		//
		// N2UR=0.0;
		// Nig=Nig; => already calculated above
		// Ni2g=0.0;
		// for each region
		//		Ni2g=0.0
		// 		for each population in region
		// 			Ni2g += nIndPerStrata
		//		Ni2g = Ni2g ^ 2
		//
		// N2UR = Ni2g/Nig;
		//  
		double N2UR=0.0;
		double Ni2g=0.0;
		for (int g=0; g < mS->hierArr[0]->nStrata; g++)
		{
			double xx=0.0;
			for (int sg=0; sg < mS->hierArr[0]->nSubStrata[g]; sg++){
				int sg_idx = mS->hierArr[0]->subStrataIdx[g][sg];
				xx += mS->hierArr[1]->nIndPerStrata[sg_idx];
			}
			Ni2g += SQUARE(xx);
		}
		N2UR=Ni2g/NUL;


		aS->ncoef[2] = (NUL - N2UR) / N1L;



		// ----------------------------------------------- //
		// VARIANCE COMPONENTS
		// ----------------------------------------------- //
		// calculate variance components (sigma squared)
		aS->sigmasq[2] = aS->msd[2];
		aS->sigmasq[1] = (aS->msd[1] - aS->sigmasq[2]) / aS->ncoef[0];
		aS->sigmasq[0] = (aS->msd[0] - aS->msd[2] - (aS->ncoef[1] * aS->sigmasq[1])) / aS->ncoef[2];
		aS->sigmasq_total = aS->sigmasq[0] + aS->sigmasq[1] + aS->sigmasq[2];

		// ----------------------------------------------- //
		// PHI STATISTIC
		// ----------------------------------------------- //
		// calculate phi statistics
		//
		// Individual,Region,Population,Total
		//

		// lvl0 (i.e. reg) in TOTAL
		// Phi_CT
		ASSERT(aS->sigmasq_total>0);
		aS->phi[0] = aS->sigmasq[0] / (aS->sigmasq_total);

		// lvl1 in lvl0
		// Phi_SC
		aS->phi[1] = aS->sigmasq[1] / (aS->sigmasq[1] + aS->sigmasq[2]);

		// lvl1 (i.e. pop) in TOTAL
		// Phi_ST
		aS->phi[2] = (aS->sigmasq[0] + aS->sigmasq[1]) / (aS->sigmasq_total);
	}
	// else
	// {
	// 	// use approximation
	// }
	// aS->print_as_table(stderr,mS);

	fprintf(stderr, "[INFO]\t-> Finished running AMOVA\n");

	return aS;
}
