#include "amova.h"

double get_SSD_total(double *M_PWD, int n_ind_cmb, int nInd)
{

	double ssd_TOTAL = 0.0;

	for (int px = 0; px < n_ind_cmb; px++)
	{
		ssd_TOTAL += SQUARE(M_PWD[px]);
	}

	ssd_TOTAL /= (double)nInd;
	return ssd_TOTAL;
}

double get_SSD_level(double *M_PWD, int n_ind_cmb, int nInd, DATA::metadataStruct *MTD, DATA::samplesStruct *SAMPLES, FILE *out_amova_ff, int sqDist, int **LUT_indPair_idx, const char *analysis_type)
{

	double sum = 0.0;
	double d_sq = 0.0;
	double ssd_WG = 0.0;
	int px;

	// group pairs by strata
	for (int sti = 0; sti < MTD->nStrata; sti++)
	{
		// fprintf(stderr, "\n-> nStrata: %i at hierarchical level ..\n", MTD->nStrata);
		// fprintf(stderr, "\n-> Strata (%s,idx:%i) has %i individuals\n",
		// 		MTD->strataArr[sti].id,
		// 		sti,
		// 		MTD->strataArr[sti].nInds);

		
		sum = 0.0;
		for (int i1 = 0; i1 < nInd - 1; i1++)
		{
			for (int i2 = i1 + 1; i2 < nInd; i2++)
			{

				//TODO do it (get associations) once, store in struct->LUT
				//TODO add print to struct as func
				// check if pair belongs to strata
				if ((SAMPLES->sampleArr[i1] & (1 << sti)) && (SAMPLES->sampleArr[i2] & (1 << sti)))
				{

#if 0
					fprintf(stderr, "\n-> Pair %i,idx:(%i,%i)) belongs to strata (%s,idx:%i)\n",
							LUT_indPair_idx[i1][i2],
							i1,
							i2,
							MTD->strataArr[sti].id,
							sti);
#endif
					px = LUT_indPair_idx[i1][i2];

					d_sq = SQUARE(M_PWD[px]);
					sum += d_sq;
				}
			}
		}
		ssd_WG += sum / (double)MTD->strataArr[sti].nInds;
	}
	return ssd_WG;
}

int AMOVA::doAMOVA(double *M_PWD, int n_ind_cmb, int nInd, DATA::metadataStruct *MTD, DATA::samplesStruct *SAMPLES, FILE *out_amova_ff, int sqDist, int **LUT_indPair_idx, const char *analysis_type)
{

	int nLevels = 1;
	amovaResultStruct *amv = new amovaResultStruct(nLevels);
	

	double ssd_TOTAL = 0.0;
	double msd_TOTAL = 0.0;
	// double sum = 0.0;
	int df_TOTAL = 0;
	// double delta_sq = 0.0;



	// ssd=[0=AG,1=AIWG,2=TOTAL]

	ssd_TOTAL = get_SSD_total(M_PWD, n_ind_cmb, nInd);
	amv->ssd[2]=ssd_TOTAL;

	df_TOTAL = nInd - 1;
	amv->df[2]=df_TOTAL;

	msd_TOTAL = ssd_TOTAL / df_TOTAL;
	amv->msd[2] = msd_TOTAL;

	double ssd_AG = 0.0;
	double ssd_WG = 0.0;
	double msd_AG = 0.0;
	int df_AG = 0;

	df_AG = MTD->nStrata - 1;
	amv->df[0] = df_AG;

	double ssd_AIWG = 0.0;
	double msd_AIWG = 0.0;

	int df_AIWG = 0;

	df_AIWG = nInd - MTD->nStrata;
	amv->df[1] = df_AIWG;

	// double s = 0.0;
	// double d_sq = 0.0;
	// int px = 0;
	ssd_WG = get_SSD_level(M_PWD, n_ind_cmb, nInd, MTD, SAMPLES, out_amova_ff, sqDist, LUT_indPair_idx, analysis_type);


	// TODO only because we have one strata level. change this
	ssd_AIWG = ssd_WG;
	amv->ssd[1] = ssd_AIWG;
	msd_AIWG = ssd_AIWG / (double)df_AIWG;
	amv->msd[1] = msd_AIWG;


	ssd_AG = ssd_TOTAL - ssd_WG;
	amv->ssd[0] = ssd_AG;
	msd_AG = ssd_AG / (double)df_AG;
	amv->msd[0] = msd_AG;

	// n variance coefficient
	// n = [ N - \sum_{g \in G} ( N^2_{g}/N) ) ]  /   G - 1
	double n_gi = 0.0;

	for (int sti = 0; sti < MTD->nStrata; sti++)
	{
		n_gi += (double)SQUARE(MTD->strataArr[sti].nInds) / (double)nInd;
	}

	// TODO double and castings are probably not necessary here
	double coef_n = (double)((double)nInd - (double)n_gi) / (double)(MTD->nStrata - 1);

	double sigmasq_a = 0.0;
	double sigmasq_b = 0.0;
	double phi_a = 0.0;
	sigmasq_a = (double)(msd_AG - msd_AIWG) / (double)coef_n;
	sigmasq_b = msd_AIWG;

	phi_a = (double)sigmasq_a / (double)(sigmasq_a + sigmasq_b);


	amv->phi[0] = phi_a;
	amv->sigmasq[0] = sigmasq_a;
	amv->sigmasq[1] = sigmasq_b;

	// print results as amova table
	// amv->print_as_table(stderr);

	amv->print_as_csv(out_amova_ff, analysis_type);


	delete amv;

	return 0;



}

