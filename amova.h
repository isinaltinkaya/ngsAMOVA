#ifndef __AMOVA__
#define __AMOVA__

#include "mathUtils.h"

namespace AMOVA
{

    /// @brief amovaStruct - stores the results of the AMOVA analysis
    /// @details
    ///
    /// levels in arrays df, ssd are as follows:
    ///     -> arr[0] = highest level
    ///         e.g. Individual ~ Region/Population/Subpopulation
    ///         ss[0] = SS within regions
    ///     -> arr[1] = lowest level
    ///     -> arr[2] = total
    ///
    /// @param df  - array of degrees of freedom
    ///               e.g. Individual ~ Region/Population
    ///            - df[0] = highest level df :: number of regions - 1
    ///            - df[2] = lowest level df
    ///            - df[3] = total df :: nInd - 1
    ///
    /// @param ss - array of sum of squares within levels
    ///               e.g. Individual ~ Region/Population
    ///           - ss[0] = highest level ss ::  SS within regions
    ///           - ss[1] = SS within populations
    ///
    /// @param ssd - array of sum of squared distances
    ///               e.g. Individual ~ Region/Population
    ///           - ssd[0] = highest level ssd :: SSD among regions
    ///               ssd[0] = SS_total - SS_within_regions
    ///           - ssd[1] = among populations within regions
    ///               ssd[1] = SS_within_regions - SS_within_populations
    ///           - ssd[2] = total ssd :: SS_total
    ///
    /// @param msd - array of mean squared distances
    /// @param ncoef - array of n coefficients
    /// @param sigmasq - array of variance components
    ///               e.g. Individual ~ Region/Population
    ///           - sigmasq[0] = variance component among regions
    ///           - sigmasq[1] = variance component among populations within regions
    ///           - sigmasq[2] = variance component among individuals within populations
    /// @param phi - array of phi statistics
    /// @param _ncoef - size of ncoef array :: nLevels+1
    /// @param _phi - size of phi array :: choose(nLevels+1, 2)
    /// @param nLevels - number of metadata levels
    /// @param nAmovaLevels - number of amova levels
    ///       equal to the number of lines in the amova table
    ///       - e.g. nLevels=2 == {region, population}
    ///           corresponds to 4 amova levels
    ///           {among regions within total,
    ///            among populations within regions,
    ///            among individuals within populations,
    ///            total}
    ///
    typedef struct amovaStruct
    {

        int *df = NULL;
        double *ss = NULL;
        double *ssd = NULL;
        double *msd = NULL;

        double *ncoef = NULL;
        double *sigmasq = NULL; // variance component
        double *phi = NULL;

        int nAmovaLevels = 0;
        int nLevels = 0;

        size_t _ncoef = 0;
        size_t _phi = 0;

        amovaStruct(metadataStruct *mS)
        {

            // storing levels in order of highest to lowest level + total
            // {level1, level2, level3, ..., total}
            // where level1 is the highest level
            // (e.g. region in Individual ~ Region/Population/Subpopulation)

            nLevels = mS->nLevels;
            // set number of amova levels
            // number of metadata levels + 2
            nAmovaLevels = nLevels + 2;
            _ncoef = nLevels + 1;
            _phi = nCk(_ncoef, 2);

            ssd = new double[nAmovaLevels];
            df = new int[nAmovaLevels];
            msd = new double[nAmovaLevels];
            for (size_t i = 0; i < (size_t)nAmovaLevels; i++)
            {
                ssd[i] = 0.0;
                df[i] = 0;
                msd[i] = 0.0;
            }

            ss = new double[_ncoef];
            ncoef = new double[_ncoef];
            sigmasq = new double[_ncoef];
            for (size_t i = 0; i < _ncoef; i++)
            {
                ncoef[i] = 0.0;
                sigmasq[i] = 0.0;
                ss[i] = 0.0;
            }

            phi = new double[_phi];
            for (size_t i = 0; i < _phi; i++)
            {
                phi[i] = 0.0;
            }
        }

        ~amovaStruct()
        {
            delete[] ssd;
            delete[] ss;
            delete[] df;
            delete[] msd;
            delete[] ncoef;
            delete[] phi;
            delete[] sigmasq;
        }

        void print_variables(FILE *fp)
        {
            fprintf(fp, "\nnAmovaLevels = %d\n", nAmovaLevels);
            fprintf(fp, "nLevels = %d\n", nLevels);
            fprintf(fp, "\n_ncoef = %zu\n", _ncoef);
            fprintf(fp, "_phi = %zu\n", _phi);
            for (size_t i = 0; i < (size_t)nAmovaLevels; i++)
            {
                fprintf(fp, "ssd[%zu] = %f\n", i, ssd[i]);
                fprintf(fp, "df[%zu] = %d\n", i, df[i]);
                fprintf(fp, "msd[%zu] = %f\n", i, msd[i]);
            }
            for (size_t i = 0; i < _ncoef; i++)
            {
                fprintf(fp, "ss[%zu] = %f\n", i, ss[i]);
                fprintf(fp, "ncoef[%zu] = %f\n", i, ncoef[i]);
                fprintf(fp, "sigmasq[%zu] = %f\n", i, sigmasq[i]);
            }
            for (size_t i = 0; i < _phi; i++)
            {
                fprintf(fp, "phi[%zu] = %f\n", i, phi[i]);
            }
        }

        // TODO print to a table file
        void print_as_table(FILE *fp, metadataStruct *mS)
        {
            fprintf(fp, "\n");
            fprintf(fp, "\n\n");
            fprintf(fp, "==========================================  AMOVA  ==========================================");
            fprintf(fp, "\n");
            fprintf(fp, "Source of variation\t\t\t\t\td.f.\tSSD\t\tMSD");
            fprintf(fp, "\n");
            fprintf(fp, "---------------------------------------------------------------------------------------------");
            fprintf(fp, "\n\n");
            fprintf(fp, "\n");

            // TODO print formula
            int x = 1;
            fprintf(fp, "Among %-15s\t\t\t\t\t%d\t%f\t%f", mS->levelNames[x], df[0], ssd[0], msd[0]);

            while (x < mS->nLevels + 1)
            {
                if (x == mS->nLevels)
                {
                    fprintf(fp, "\n");
                    fprintf(fp, "Among %s within %-25s\t%d\t%f\t%f", mS->levelNames[0], mS->levelNames[mS->nLevels], df[x], ssd[x], msd[x]);
                    fprintf(fp, "\n");
                }
                else
                {
                    fprintf(fp, "\n");
                    fprintf(fp, "Among %s within %-25s\t%d\t%f\t%f", mS->levelNames[x + 1], mS->levelNames[x], df[x], ssd[x], msd[x]);
                    fprintf(fp, "\n");
                }
                x++;
            }

            fprintf(fp, "\n");
            fprintf(fp, "Total\t\t\t\t\t\t\t%d\t%f\t%f", df[nAmovaLevels - 1], ssd[nAmovaLevels - 1], msd[nAmovaLevels - 1]);
            fprintf(fp, "\n\n\n");
            fprintf(fp, "Variance components:\n\n");
            for (size_t i = 0; i < _ncoef - 1; i++)
            {
                fprintf(fp, "\n\t%-20s", mS->levelNames[i + 1]);
                fprintf(fp, "\t%f", sigmasq[i]);
            }
            // Lowest level (i.e. Individual)
            fprintf(fp, "\n\t%-20s", mS->levelNames[0]);
            fprintf(fp, "\t%f", sigmasq[_ncoef - 1]);

            fprintf(fp, "\n\n\n");
            fprintf(fp, "\nVariance coefficients:\n\n\t");
            if (nLevels == 1)
            {
                fprintf(fp, "%f", ncoef[0]);
            }
            else
            {
                for (size_t i = 0; i < _ncoef; i++)
                {
                    fprintf(fp, "%f\t", ncoef[i]);
                }
            }
            fprintf(fp, "\n\n\n");
            fprintf(fp, "Phi-statistic:\n\n");
            for (size_t i = 0; i < _phi; i++)
            {
                fprintf(fp, "\t%f", phi[i]);
            }
            fprintf(fp, "\n\n");
            fprintf(fp, "=============================================================================================");
            fprintf(fp, "\n\n");
        }

        void print_as_csv(FILE *fp, metadataStruct *mS)
        {

            // TODO add percentage total?
            // header
            //  type,label,value
            //  SSD,Among_region,0.1234
            //  fprintf(fp, "type,label,value\n");
            fprintf(fp, "df,Total,%d\n", df[nAmovaLevels - 1]);
            fprintf(fp, "SSD,Total,%f\n", ssd[nAmovaLevels - 1]);
            fprintf(fp, "MSD,Total,%f\n", msd[nAmovaLevels - 1]);

            int x = 1;
            fprintf(fp, "df,Among_%s_within_%s,%d\n", mS->levelNames[x], "Total", df[0]);
            fprintf(fp, "SSD,Among_%s_within_%s,%f\n", mS->levelNames[x], "Total", ssd[0]);
            fprintf(fp, "MSD,Among_%s_within_%s,%f\n", mS->levelNames[x], "Total", msd[0]);
            while (x < mS->nLevels + 1)
            {
                if (x == mS->nLevels)
                {
                    fprintf(fp, "df,Among_%s_within_%s,%d\n", mS->levelNames[0], mS->levelNames[mS->nLevels], df[x]);
                    fprintf(fp, "SSD,Among_%s_within_%s,%f\n", mS->levelNames[0], mS->levelNames[mS->nLevels], ssd[x]);
                    fprintf(fp, "MSD,Among_%s_within_%s,%f\n", mS->levelNames[0], mS->levelNames[mS->nLevels], msd[x]);
                }
                else
                {
                    fprintf(fp, "df,Among_%s_within_%s,%d\n", mS->levelNames[x + 1], mS->levelNames[x], df[x]);
                    fprintf(fp, "SSD,Among_%s_within_%s,%f\n", mS->levelNames[x + 1], mS->levelNames[x], ssd[x]);
                    fprintf(fp, "MSD,Among_%s_within_%s,%f\n", mS->levelNames[x + 1], mS->levelNames[x], msd[x]);
                }
                x++;
            }

            if (nLevels == 1)
            {
                fprintf(fp, "Phi,%s_in_%s,%f\n", mS->levelNames[1], "Total", phi[0]);
                fprintf(fp, "Variance_coefficient,a,%f\n", ncoef[0]);
                fprintf(fp, "Variance_component,%s,%f\n", mS->levelNames[1], sigmasq[0]);
                fprintf(fp, "Variance_component,%s,%f\n", mS->levelNames[0], sigmasq[1]);
            }
            else if (nLevels == 2)
            {
                fprintf(fp, "Phi,%s_in_%s,%f\n", mS->levelNames[1], "Total", phi[0]);
                fprintf(fp, "Phi,%s_in_%s,%f\n", mS->levelNames[2], mS->levelNames[1], phi[1]);
                fprintf(fp, "Phi,%s_in_%s,%f\n", mS->levelNames[2], "Total", phi[2]);
                fprintf(fp, "Variance_coefficient,a,%f\n", ncoef[0]);
                fprintf(fp, "Variance_coefficient,b,%f\n", ncoef[1]);
                fprintf(fp, "Variance_coefficient,c,%f\n", ncoef[2]);
                fprintf(fp, "Variance_component,%s,%f\n", mS->levelNames[1], sigmasq[0]);
                fprintf(fp, "Variance_component,%s,%f\n", mS->levelNames[2], sigmasq[1]);
                fprintf(fp, "Variance_component,%s,%f\n", mS->levelNames[0], sigmasq[2]);
            }
            else
            {
                fprintf(stderr, "[ERROR]: nLevels > 2 not supported yet\n");
            }
        }

    } amovaStruct;

    int doAMOVA(distanceMatrixStruct *dMS,
                metadataStruct *MTD,
                sampleStruct *SAMPLES, FILE *out_amova_ff, int **lut_indsToIdx);

    amovaStruct *amovaStruct_doAmova(distanceMatrixStruct *dMS,
                                     metadataStruct *MTD,
                                     sampleStruct *SAMPLES, int **lut_indsToIdx);

}

#endif