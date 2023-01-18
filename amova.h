#include "shared.h"
#include "math_utils.h"




namespace AMOVA {

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
    typedef struct amovaStruct{

        int* df = NULL;
        double* ss = NULL;
        double* ssd = NULL;
        double* msd = NULL;

        double* ncoef = NULL;
        double* sigmasq = NULL; //variance component
        double* phi = NULL;

        
        int nAmovaLevels = 0;
        int nLevels = 0;

        size_t _ncoef = 0;
        size_t _phi = 0;

        amovaStruct(DATA::metadataStruct* mS){

            // storing levels in order of highest to lowest level + total
            // {level1, level2, level3, ..., total}
            // where level1 is the highest level 
            // (e.g. region in Individual ~ Region/Population/Subpopulation)

            nLevels = mS->nLevels;
            // set number of amova levels
            // number of metadata levels + 2
            nAmovaLevels = nLevels + 2;
            _ncoef=nLevels+1;
            _phi = nCk(_ncoef, 2);

            ssd = new double[nAmovaLevels];
            df = new int[nAmovaLevels];
            msd = new double[nAmovaLevels];
            for (size_t i=0; i<(size_t) nAmovaLevels; i++){
                ssd[i] = 0.0;
                df[i] = 0;
                msd[i] = 0.0;
            }


            ss = new double[_ncoef];
            ncoef = new double[_ncoef];
            sigmasq = new double[_ncoef];
            for (size_t i=0; i<_ncoef; i++){
                ncoef[i] = 0.0;
                sigmasq[i] = 0.0;
                ss[i] = 0.0;
            }

            phi = new double[_phi];
            for (size_t i=0; i< _phi; i++){
                phi[i] = 0.0;
            }

        }

        ~amovaStruct(){
            delete[] ssd;
            delete[] ss;
            delete[] df;
            delete[] msd;
            delete[] ncoef;
            delete[] phi;
            delete[] sigmasq;

        }

        void print_variables(FILE *fp){
            fprintf(fp, "\nnAmovaLevels = %d\n", nAmovaLevels);
            fprintf(fp, "nLevels = %d\n", nLevels);
            fprintf(fp, "\n_ncoef = %zu\n", _ncoef);
            fprintf(fp, "_phi = %zu\n", _phi);
            for (size_t i=0; i<(size_t)nAmovaLevels; i++){
                fprintf(fp, "ssd[%zu] = %f\n", i, ssd[i]);
                fprintf(fp, "df[%zu] = %d\n", i, df[i]);
                fprintf(fp, "msd[%zu] = %f\n", i, msd[i]);
            }
            for (size_t i=0; i<_ncoef; i++){
                fprintf(fp, "ss[%zu] = %f\n", i, ss[i]);
                fprintf(fp, "ncoef[%zu] = %f\n", i, ncoef[i]);
                fprintf(fp, "sigmasq[%zu] = %f\n", i, sigmasq[i]);
            }
            for(size_t i=0; i<_phi; i++){
                fprintf(fp, "phi[%zu] = %f\n", i, phi[i]);
            }
        }

        void print_as_table(FILE *fp, DATA::metadataStruct* mS){
            fprintf(fp,"\n");
            fprintf(fp,"\n\n");
            fprintf(fp,"==========================================  AMOVA  =========================================="); 
            fprintf(fp,"\n");
            fprintf(fp,"Source of variation\t\t\t\t\td.f.\tSSD\t\tMSD");
            fprintf(fp,"\n");
            fprintf(fp,"---------------------------------------------------------------------------------------------");
            fprintf(fp,"\n\n");
            fprintf(fp,"\n");

            int x=1;
            fprintf(fp, "Among %-15s\t\t\t\t\t%d\t%f\t%f", mS->levelNames[x], df[0], ssd[0], msd[0]);

            while(x<mS->nLevels+1){
                if(x==mS->nLevels){
                    fprintf(fp,"\n");
                    fprintf(fp, "Among %s within %-25s\t%d\t%f\t%f", mS->levelNames[0], mS->levelNames[mS->nLevels], df[x], ssd[x], msd[x]);
                    fprintf(fp,"\n");
                }else{
                    fprintf(fp,"\n");
                    fprintf(fp, "Among %s within %-25s\t%d\t%f\t%f", mS->levelNames[x+1], mS->levelNames[x], df[x], ssd[x], msd[x]);
                    fprintf(fp,"\n");
                }
                x++;
            }

            fprintf(fp,"\n");
            fprintf(fp,"Total\t\t\t\t\t\t\t%d\t%f\t%f", df[nAmovaLevels-1], ssd[nAmovaLevels-1], msd[nAmovaLevels-1]);
            fprintf(fp,"\n\n\n");
            fprintf(fp,"Variance components:");
            fprintf(fp,"\n\n");
            for(size_t i=0; i<_ncoef; i++){
                fprintf(fp,"\n");
                fprintf(fp,"sigma^2");
                fprintf(fp,"\t");
                // fprintf(fp,"%20s", mS->levelNames[i]);
                fprintf(fp,"\t");
                fprintf(fp,"%f", sigmasq[i]);
                fprintf(fp,"\n");
                fprintf(fp,"\n\t%f",ncoef[i]);
            }
            fprintf(fp,"\n\n");
            fprintf(fp,"Phi-statistic:");
            for(size_t i=0; i<_phi; i++){
                fprintf(fp,"\n");
                fprintf(fp,"\t%f", phi[i]);
            }
            fprintf(fp,"\n\n");
            fprintf(fp,"============================================================================================="); 
            fprintf(fp,"\n\n");

        }

        void print_as_csv(FILE *fp, const char *analysis_type){
            fprintf(fp, "%s,%i,%f,%f,%i,%f,%f,%i,%f,%f,%f,%f,%f,%f\n",
                analysis_type, 
                df[0], ssd[0], msd[0],
                df[1], ssd[1], msd[1],
                df[2], ssd[2], msd[2],
                ncoef[0], sigmasq[0], sigmasq[1], phi[0]);
        }

    
    }amovaStruct;

    int doAMOVA(DATA::distanceMatrixStruct *dMS,
        DATA::metadataStruct *MTD, 
        DATA::samplesStruct *SAMPLES, FILE *out_amova_ff, int **LUT_indPair_idx, const char *type);



}