#include "shared.h"
#include "math_utils.h"




namespace AMOVA {

    typedef struct amovaResultStruct{

        // initialize to lowest possible values, 1 level case
        int* df = NULL;
        double* ssd = NULL;
        double* msd = NULL;
        double* coef_n = NULL;
        double* sigmasq = NULL;
        double* phi = NULL;



        int nLevels=1;
        size_t _df;
        size_t _ssd;
        size_t _msd;
        size_t _coef_n;
        size_t _sigmasq;
        size_t _phi;

        amovaResultStruct(int nLevels_){

            _df = 3;
            _ssd = 3;
            _msd = 3;
            _coef_n = 1;
            _sigmasq = 2;
            _phi = 1;

            nLevels = nLevels_;
            df = (int*) malloc( _df * sizeof(int));
            ssd = (double*) malloc( _ssd * sizeof(double));
            msd = (double*) malloc( _msd * sizeof(double));
            coef_n = (double*) malloc( _coef_n * sizeof(double));
            sigmasq = (double*) malloc( _sigmasq * sizeof(double));
            phi = (double*) malloc( _phi * sizeof(double));

            for (size_t i=0; i<_df; i++) df[i] = 0;
            for (size_t i=0; i<_ssd; i++) ssd[i] = 0;
            for (size_t i=0; i<_msd; i++) msd[i] = 0;
            for (size_t i=0; i<_coef_n; i++) coef_n[i] = 0;
            for (size_t i=0; i<_sigmasq; i++) sigmasq[i] = 0;
            for (size_t i=0; i<_phi; i++) phi[i] = 0;

        }

        ~amovaResultStruct(){
            free(df);
            free(ssd);
            free(msd);
            free(coef_n);
            free(sigmasq);
            free(phi);
        }

        void print_as_table(FILE *fp){

            fprintf(fp,"\n");
            fprintf(fp,"\n");
            fprintf(fp,"\n");
            fprintf(fp,"==========================================  AMOVA  =========================================="); 
            fprintf(fp,"\n");
            fprintf(fp,"Source of variation\t\t\td.f.\tSSD\t\tMSD");
            fprintf(fp,"\n");
            fprintf(fp,"---------------------------------------------------------------------------------------------");
            fprintf(fp,"\n");
            fprintf(fp,"\n");
            fprintf(fp,"Among groups");
            fprintf(fp,"\t\t\t\t");
            fprintf(fp,"%d",df[0]);
            fprintf(fp,"\t");
            fprintf(fp,"%f",ssd[0]);
            fprintf(fp,"\t");
            fprintf(fp,"%f",msd[0]);
            fprintf(fp,"\n");
            fprintf(fp,"Among individuals within groups");
            fprintf(fp,"\t\t");
            fprintf(fp,"%d",df[1]);
            fprintf(fp,"\t");
            fprintf(fp,"%f",ssd[1]);
            fprintf(fp,"\t");
            fprintf(fp,"%f",msd[1]);
            fprintf(fp,"\n");
            fprintf(fp,"\n");
            fprintf(fp,"\n");
            fprintf(fp,"Total");
            fprintf(fp,"\t\t\t\t\t");
            fprintf(fp,"%d",df[2]);
            fprintf(fp,"\t");
            fprintf(fp,"%f",ssd[2]);
            fprintf(fp,"\t");
            fprintf(fp,"%f",msd[2]);
            fprintf(fp,"\n");
            fprintf(fp,"\n");
            fprintf(fp,"\n");
            fprintf(fp,"\n");
            fprintf(fp,"Variance coefficients:");
            fprintf(fp,"\n");
            fprintf(fp,"a");
            fprintf(fp,"\t");
            fprintf(fp,"%f", sigmasq[0]);
            fprintf(fp,"\n");
            fprintf(fp,"\n");
            fprintf(fp,"\n");
            fprintf(fp,"Phi-statistic:");
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
                coef_n[0], sigmasq[0], sigmasq[1], phi[0]);
        }

    
    }amovaResultStruct;

    int doAMOVA(DATA::distanceMatrixStruct *dMS,
        DATA::metadataStruct *MTD, 
        DATA::samplesStruct *SAMPLES, FILE *out_amova_ff, int **LUT_indPair_idx, const char *type);



}