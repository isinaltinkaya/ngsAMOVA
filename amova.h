#include "shared.h"
#include "math_utils.h"




namespace AMOVA {

    typedef struct amovaResultStruct{

        // initialize to lowest possible values, 1 level case
        int* df = NULL;
        double* ssd = NULL;
        double* msd = NULL;

        double* ncoef = NULL;

        double* sigmasq = NULL; //variance component
        double* phi = NULL;

        int nAmovaLevels = 0;


        size_t _ssd;
        size_t _ncoef;
        size_t _phi;

        amovaResultStruct(DATA::metadataStruct* mS){

            // set number of amova levels
            nAmovaLevels = mS->nLevels + 1;

            // number of SSD values
            // = number of amova levels + 1 (total)
            _ssd = nAmovaLevels + 1;
            ssd = new double[_ssd];
            df = new int[_ssd];
            msd = new double[_ssd];
            for (size_t i=0; i<_ssd; i++){
                ssd[i] = 0;
                df[i] = 0;
                msd[i] = 0;
            }



            _ncoef = nCk(nAmovaLevels,2);
            ncoef = new double[_ncoef];
            phi = new double[_ncoef];
            for (size_t i=0; i<_ncoef; i++){
                ncoef[i] = 0;
                phi[i] = 0;
            }

            sigmasq = new double[nAmovaLevels];
            for (size_t i=0; i<nAmovaLevels; i++){
                sigmasq[i] = 0;
            }

        }

        ~amovaResultStruct(){
            delete[] ssd;
            delete[] df;
            delete[] msd;
            delete[] ncoef;
            delete[] phi;
            delete[] sigmasq;

        }

        void print_variables(FILE *fp){
            fprintf(stderr, "\n_ssd = %zu\n", _ssd);
            fprintf(stderr, "_ncoef = %zu\n", _ncoef);
            fprintf(stderr, "nAmovaLevels = %d\n", nAmovaLevels);
            for (size_t i=0; i<nAmovaLevels; i++){
                fprintf(fp, "sigmasq[%zu] = %f\n", i, sigmasq[i]);
            }
            for (size_t i=0; i<_ssd; i++){
                fprintf(fp, "ssd[%zu] = %f\n", i, ssd[i]);
                fprintf(fp, "df[%zu] = %d\n", i, df[i]);
                fprintf(fp, "msd[%zu] = %f\n", i, msd[i]);
            }
            for (size_t i=0; i<_ncoef; i++){
                fprintf(fp, "ncoef[%zu] = %f\n", i, ncoef[i]);
                fprintf(fp, "phi[%zu] = %f\n", i, phi[i]);
            }
        }

        void print_as_table(FILE *fp, DATA::metadataStruct* mS){
            fprintf(fp,"\n");
            fprintf(fp,"\n");
            fprintf(fp,"\n");
            fprintf(fp,"==========================================  AMOVA  =========================================="); 
            fprintf(fp,"\n");
            fprintf(fp,"Source of variation\t\t\t\t\td.f.\tSSD\t\tMSD");
            fprintf(fp,"\n");
            fprintf(fp,"---------------------------------------------------------------------------------------------");
            fprintf(fp,"\n");
            fprintf(fp,"\n");
            int x=1;
            fprintf(fp,"\n");
            fprintf(fp, "Among %-15s", mS->levelNames[x]);
            fprintf(fp,"\t\t\t\t");
            fprintf(fp,"\t");
            fprintf(fp,"%d",df[x-1]);
            fprintf(fp,"\t");
            fprintf(fp,"%f",ssd[x-1]);
            fprintf(fp,"\t");
            fprintf(fp,"%f",msd[x-1]);
            fprintf(fp,"\n");

            while(x<mS->nLevels+1){
                if(x==mS->nLevels){
                    fprintf(fp,"\n");
                    fprintf(fp, "Among %s within %-25s\t%d\t%f\t%f", mS->levelNames[0], mS->levelNames[mS->nLevels], df[x], ssd[x], msd[x]);
                    // fprintf(fp,"\t\t\t");
                    fprintf(fp,"\n");
                }else{
                    fprintf(fp,"\n");
                    fprintf(fp, "Among %s within %-25s\t%d\t%f\t%f", mS->levelNames[x+1], mS->levelNames[x], df[x], ssd[x], msd[x]);
                    fprintf(fp,"\n");
                }
                x++;

            }

            fprintf(fp,"\n");
            fprintf(fp,"Total\t\t\t\t\t\t\t%d\t%f\t%f", df[nAmovaLevels], ssd[nAmovaLevels], msd[nAmovaLevels]);
            fprintf(fp,"\n\n\n");
            fprintf(fp,"Variance components:");
            fprintf(fp,"\n\n");
            int i=nAmovaLevels;
            while(i>0){
                fprintf(fp,"sigma^2");
                fprintf(fp,"\t");
                fprintf(fp,"%s", mS->levelNames[i-1]);
                fprintf(fp,"\t");
                fprintf(fp,"%f", sigmasq[nAmovaLevels-i]);
                fprintf(fp,"\n");
                i--;
            }



            //TODO add sigmasqb
            //todo add this nCk style sigma_which_which
            fprintf(fp,"\n");
            fprintf(fp,"Variance coefficients:");
            for(int i=0; i<_ncoef; i++){
                fprintf(fp,"\n");
                // fprintf(fp,"a");
                fprintf(fp,"\t");
                fprintf(fp,"%f",ncoef[i]);
            }
            fprintf(fp,"\n");
            fprintf(fp,"\n");
            // fprintf(fp,"a");
            // fprintf(fp,"\t");
            // fprintf(fp,"%f",ncoef[0]);
            fprintf(fp,"\n");
            fprintf(fp,"\n");
            fprintf(fp,"\n");
            fprintf(fp,"Phi-statistic:");
            //TODO x in y
            fprintf(fp,"\n");
            fprintf(fp,"a");
            fprintf(fp,"\t");
            fprintf(fp,"%f",phi[0]);
            fprintf(fp,"\n");
            fprintf(fp,"\n");
            fprintf(fp,"============================================================================================="); 
            fprintf(fp,"\n");
            fprintf(fp,"\n");
            fprintf(fp,"\n");

        }

        void print_as_csv(FILE *fp, const char *analysis_type){
            fprintf(fp, "%s,%i,%f,%f,%i,%f,%f,%i,%f,%f,%f,%f,%f,%f\n",
                analysis_type, 
                df[0], ssd[0], msd[0],
                df[1], ssd[1], msd[1],
                df[2], ssd[2], msd[2],
                ncoef[0], sigmasq[0], sigmasq[1], phi[0]);
        }

    
    }amovaResultStruct;

    int doAMOVA(DATA::distanceMatrixStruct *dMS,
        DATA::metadataStruct *MTD, 
        DATA::samplesStruct *SAMPLES, FILE *out_amova_ff, int **LUT_indPair_idx, const char *type);



}