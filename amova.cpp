#include "amova.h"

double get_SSD_total(DATA::distanceMatrixStruct *dMS)
{

	double ssd_TOTAL = 0.0;

	for (int px = 0; px < dMS->nIndCmb; px++)
	{
		ssd_TOTAL += SQUARE(dMS->M[px]);
	}

	ssd_TOTAL /= (double)dMS->nInd;
	return ssd_TOTAL;
}

double get_SSD_level(int lvl_i, DATA::distanceMatrixStruct *dMS, DATA::metadataStruct *MTD, DATA::samplesStruct *SAMPLES, FILE *out_amova_ff, int **LUT_indPair_idx, const char *analysis_type)
{

	double sum = 0.0;
	double d_sq = 0.0;
	double ssd_WG = 0.0;
	int px;

	for (int sti=0; sti < MTD->hierArr[lvl_i]->nStrata; sti++)
	{
		// fprintf(stderr, "\n-> Strata (%s,idx:%i) has %i individuals\n", MTD->hierArr[lvl_i]->strataNames[sti], sti, MTD->hierArr[lvl_i]->nIndPerStrata[sti]);
		// fprintf(stderr, "\n\n\n --- sti = %i ---\n\n\n", sti);
		// MTD->print(stderr);

		sum = 0.0;
		for (int i1 = 0; i1 < dMS->nInd - 1; i1++)
		{
			for (int i2 = i1 + 1; i2 < dMS->nInd; i2++)
			{
				if(MTD->cmp_assoc_at_lvl(lvl_i, i1, i2) == 0){
#if 0
					fprintf(stderr, "\n-> Pair %i,idx:(%i,%i)) belongs to strata (%s,idx:%i)\n",
							LUT_indPair_idx[i1][i2],
							i1,
							i2,
							MTD->strataArr[sti].id,
							sti);
#endif
					px = LUT_indPair_idx[i1][i2];

					d_sq = SQUARE(dMS->M[px]);
					sum += d_sq;
				}
			}
		}
		ssd_WG += sum / (double)MTD->hierArr[lvl_i]->nIndPerStrata[sti];
	}
	return ssd_WG;
}

int AMOVA::doAMOVA(DATA::distanceMatrixStruct *dMS, DATA::metadataStruct *MTD, DATA::samplesStruct *SAMPLES, FILE *out_amova_ff, int **LUT_indPair_idx, const char *analysis_type)
{



	int nLevels = 1;
	amovaResultStruct *amv = new amovaResultStruct(nLevels);

	double ssd_TOTAL = 0.0;
	double msd_TOTAL = 0.0;
	// double sum = 0.0;
	int df_TOTAL = 0;
	// double delta_sq = 0.0;

	// ssd=[0=AG,1=AIWG,2=TOTAL]

	ssd_TOTAL = get_SSD_total(dMS);
	amv->ssd[2] = ssd_TOTAL;

	df_TOTAL = dMS->nInd - 1;
	amv->df[2] = df_TOTAL;

	msd_TOTAL = ssd_TOTAL / df_TOTAL;
	amv->msd[2] = msd_TOTAL;

	double ssd_AG = 0.0;
	double ssd_WG = 0.0;
	double msd_AG = 0.0;
	int df_AG = 0;

	// df_AG = MTD->nLevels - 1;
	df_AG = MTD->hierArr[0]->nStrata - 1;

	amv->df[0] = df_AG;

	double ssd_AIWG = 0.0;
	double msd_AIWG = 0.0;

	int df_AIWG = 0;

	// df_AIWG = nInd - MTD->nLevels;
	df_AG = dMS->nInd - MTD->hierArr[0]->nStrata;
	amv->df[1] = df_AIWG;

	amv->print_as_table(stderr);


	
	// double s = 0.0;
	// double d_sq = 0.0;
	// int px = 0;
	int lvl_i=0;
	ssd_WG = get_SSD_level(lvl_i, dMS, MTD, SAMPLES, out_amova_ff, LUT_indPair_idx, analysis_type);
	fprintf(stderr, "\n\n\nssd_WG = %f\n\n\n", ssd_WG);

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

	for (int lvl_i = 0; lvl_i < MTD->nLevels; lvl_i++)
	{
		for (int sti = 0; sti < MTD->hierArr[lvl_i]->nStrata; sti++)
		{
			n_gi += (double)SQUARE(MTD->hierArr[lvl_i]->nIndPerStrata[sti]) / (double)dMS->nInd;
		}
	}

	// TODO double and castings are probably not necessary here
	double coef_n = (double)((double)dMS->nInd - (double)n_gi) / (double)(MTD->nLevels - 1);

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
