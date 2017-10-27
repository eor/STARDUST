#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_integration.h>

#include "table_settings.h"
#include "constants.h"

#include "prototype.h"
#include "allvars.h"
#include "log.h"

/* This program generates the table required for Compton heating 

 REFERENCES:
 [1] Dijkstra, M., Haiman, Z. 2004,ApJ, 613, 646-654 
 [2] Fukugita, M & Kawasaki, M. 1994, MNRAS 269, 563
 [3] Madau P., Efstathiou G., 1999,  ApJ, 517, L9-L12
 [4] Zaroubi, S & Silk, J. 2005, MNRAS (TBS)
 [5] Zhao D.H., Mo H.J., Jing Y.P., Borner G. 2003, MNRAS,339,12.
*/ 
#ifdef POWERLAWSRC
#define Ec1 2e3
#define Ec2 1e4           /** cut-off Ref (Zaroubi & Silk, 2005) **/
#endif

/***************************************************************
 * global variables (scope - this file)
 ***************************************************************/
double comp_term;

int table_compton (void){
  
    printf(" Creating tables: COMPTON\n");
    if(DEBUG)log_debug("Starting table ctreation (COMPTON)");
    
    
    double l1,l2,l3,l4;

    double Cell_Max_Log1,Cell_Max_Log2;
    double previous = 100.; /* Some relatively large number */

    double absError = 1e-3;     /* integration errors */
    double relError = 1e-3;
    
    l1  = log(HIONIZEeV);
    l2  = log(He1IONIZEeV);
    l3  = log(He2IONIZEeV);
    l4  = log(myConfig.sourceEHigh); 
    
    size_t intResolution = 100;

    Cell_Max_Log1 = UPLIM1+TABLERES;
    Cell_Max_Log2 = UPLIM2+TABLERES;

    /* Setting up the GSL integration environment */
    gsl_integration_workspace *wtable = gsl_integration_workspace_alloc (intResolution);

    double result1,result2, error;

    gsl_function Func1;
    gsl_function Func2;


    comp_term = THOMSON_ELEC_CROSS / MASSELECTRONeV ; /* [cm^2 eV^-1] */
    
    FILE *fpComp1, *fpComp2;
    char *tmpPath1 = NULL;
    char *tmpPath2 = NULL;
        
    const char *tableFile_comp1_p1 = "table_comp1_p1.dat";
    const char *tableFile_comp1_p2 = "table_comp1_p2.dat";
    const char *tableFile_comp1_p3 = "table_comp1_p3.dat";
    
    const char *tableFile_comp2_p1 = "table_comp2_p1.dat";
    const char *tableFile_comp2_p2 = "table_comp2_p2.dat";
    const char *tableFile_comp2_p3 = "table_comp2_p3.dat";
  
    
    /***************************************************************
     * Compton 1 & 2 tables for first intervall (p1)
     ***************************************************************/     

    tmpPath1 = utils_concat_path( myConfig.pathOutDir, (char*) tableFile_comp1_p1 );
    tmpPath2 = utils_concat_path( myConfig.pathOutDir, (char*) tableFile_comp2_p1 );
    fpComp1  = fopen( tmpPath1 , "w");       
    fpComp2  = fopen( tmpPath2 , "w");
    
    
    Func1.function = &compton_heating1_p1;
    Func1.params = 0;

    Func2.function = &compton_heating2_p1;
    Func2.params = 0;

    for(nHX1=LOWLIM;nHX1<=Cell_Max_Log1;nHX1+=TABLERES){
        
        if(previous>1.e-100){
            
            /* Integration */
            gsl_integration_qags( &Func1, l1, l2, absError, relError, intResolution, wtable, &result1, &error );
            gsl_integration_qags( &Func2, l1, l2, absError, relError, intResolution, wtable, &result2, &error );
            
            previous = result1;

            fprintf(fpComp1,"%e \t %e \n",nHX1,previous);
            fprintf(fpComp2,"%e \t %e \n",nHX1,result2);
        }
        else{
            fprintf(fpComp1,"%e \t %e\n",nHX1,1.e-100);
            fprintf(fpComp2,"%e \t %e\n",nHX1,1.e-100); 
        }
    }
    
    fclose(fpComp1);
    fclose(fpComp2);

    if(DEBUG)log_debug(" Computed Compton tables for H in first energy interval (l1-l2)");

    /***************************************************************
     * Compton 1 & 2 tables for second intervall (p2)
     ***************************************************************/     
 
    tmpPath1 = utils_concat_path( myConfig.pathOutDir, (char*) tableFile_comp1_p2 );
    tmpPath2 = utils_concat_path( myConfig.pathOutDir, (char*) tableFile_comp2_p2 );     
    fpComp1  = fopen( tmpPath1 , "w");       
    fpComp2  = fopen( tmpPath2 , "w");
    
    Func1.function = &compton_heating1_p2;
    Func1.params = 0;

    Func2.function = &compton_heating2_p2;
    Func2.params = 0;

    for(nHX1=LOWLIM;nHX1<=Cell_Max_Log1;nHX1+=TABLERES){  
        
        nHeX2 = LOWLIM;
        
        /* Integration */
        gsl_integration_qags (&Func1, l2, l3, absError, relError, intResolution, wtable, &result1, &error);
        gsl_integration_qags (&Func2, l2, l3, absError, relError, intResolution, wtable, &result2, &error);

        if(result1 > 1.e-100)
            previous = 100.;
        else 
            previous = 1.e-200;
        
        for(nHeX2=LOWLIM;nHeX2<=Cell_Max_Log1;nHeX2+=TABLERES){
                
            if(previous>1.e-100){
                previous = result1;

                fprintf(fpComp1,"%e\t%e\t %e\n",nHX1,nHeX2,previous);
                fprintf(fpComp2,"%e\t%e\t %e\n",nHX1,nHeX2,result2);
            }
            else{
                fprintf(fpComp1,"%e\t%e\t %e\n",nHX1,nHeX2,1.e-100);
                fprintf(fpComp2,"%e\t%e\t %e\n",nHX1,nHeX2,1.e-100);
            }
        }
    }

    fclose(fpComp1);
    fclose(fpComp2);
    
    if(DEBUG)log_debug(" Computed Compton tables for H, He in second energy interval (l2-l3)");
    
    /***************************************************************
     * Compton 1 & 2 tables for third intervall (p3)
     ***************************************************************/     
    
    tmpPath1 = utils_concat_path( myConfig.pathOutDir, (char*) tableFile_comp1_p3 );
    tmpPath2 = utils_concat_path( myConfig.pathOutDir, (char*) tableFile_comp2_p3 );
    fpComp1  = fopen( tmpPath1 , "w");       
    fpComp2  = fopen( tmpPath2 , "w"); 
    
    Func1.function = &compton_heating1_p3;
    Func1.params = 0;

    Func2.function = &compton_heating2_p3;
    Func2.params = 0;

    for(nHX13=LOWLIM;nHX13<=Cell_Max_Log2;nHX13+=TABLERES){ 
        
        nHeX2 = LOWLIM;
        
        /* Integration */
        gsl_integration_qags (&Func1, l3, l4, absError, relError, intResolution, wtable, &result1, &error);
        gsl_integration_qags (&Func2, l3, l4, absError, relError, intResolution, wtable, &result2, &error);

        if(result1 > 1.e-100)
            previous = 100.;
        else 
            previous = 1.e-200;

        for(nHeX2=LOWLIM;nHeX2<=Cell_Max_Log1;nHeX2+=TABLERES){
            
            if(previous>1.e-100){
                previous = result1;

                fprintf(fpComp1,"%e\t%e\t%e\n",nHX13,nHeX2,previous);
                fprintf(fpComp2,"%e\t%e\t%e\n",nHX13,nHeX2,result2);
            }
            else{
                fprintf(fpComp1,"%e\t%e\t%e\n",nHX13,nHeX2,1.e-100);
                fprintf(fpComp2,"%e\t%e\t%e\n",nHX13,nHeX2,1.e-100);
            }
        }
    }

    fclose(fpComp1);
    fclose(fpComp2);

    if(DEBUG)log_debug(" Computed Compton tables for H, He in third energy interval (l3-l4)");

    /* freeing the GSL integration environment */
    gsl_integration_workspace_free (wtable);
    
    if(DEBUG)log_debug("Leaving table creation (COMPTON)");    
    
    return 1;
}


/***************************************************************
 * NOTE: Integration with the functions below is done using 
 * the logarithmic integration trick ... 
 * hence the extra "energy term = egamm" in the equations.
 ***************************************************************/

double compton_heating1_p1(double egamm,void * params){

    double nHx1;    
    double decay_term,F_e1h1;
    
    
    egamm   = exp(egamm);
    nHx1    = exp(nHX1);
    F_e1h1  = cross_sec_e1h1(egamm); //---eq(14) ref ([3])
    decay_term = exp(-(nHx1*F_e1h1)*drKpc); 

    return egamm*(egamm*comp_term*sed_get_value(egamm)*decay_term) ;

  }


double compton_heating1_p2(double egamm,void * params)
{
 egamm = exp(egamm);
 double nHx1,nHex2;
 nHx1 = exp(nHX1);
 nHex2 = exp(nHeX2);

 double decay_term,F_e1h1,F_ehe1;

 F_e1h1=cross_sec_e1h1(egamm);
 F_ehe1=cross_sec_ehe1(egamm);
   //---eq(14) ref ([3])
 decay_term = exp(-(nHx1*F_e1h1+nHex2*F_ehe1)*drKpc);  

 return  egamm*(egamm*comp_term*sed_get_value(egamm)*decay_term);

  }

double compton_heating1_p3(double egamm,void * params){
    
    egamm = exp(egamm);
    double nHx13,nHex2;
    nHx13 = exp(nHX13);
    nHex2 = exp(nHeX2);
    double decay_term,F_e1h1,F_ehe1;
    F_e1h1=cross_sec_e1h1(egamm);
    F_ehe1=cross_sec_ehe1(egamm);
    //---eq(14) ref ([3])
    decay_term = exp(-(nHx13*F_e1h1+nHex2*F_ehe1)*drKpc);
    /* nHx13 = (xH + (k/12.5)*xHeII) got from the program */

    return  egamm*(egamm*comp_term*sed_get_value(egamm)*decay_term) ;

}

double compton_heating1_p4(double egamm,void * params){
    
    egamm = exp(egamm);
    double nHx13,nHex2;
    nHx13 = exp(nHX13);
    nHex2 = exp(nHeX2);
    double decay_term,F_e1h1,F_ehe1;
    F_e1h1=cross_sec_e1h1(egamm);
    F_ehe1=cross_sec_ehe1(egamm);
    //*exp(-egamm/Ec2);

    decay_term = exp(-(nHx13*F_e1h1+nHex2*F_ehe1)*drKpc);
    /* nHx13 = (xH + (k/12.5)*xHeII) got from the program */

    return  egamm*(egamm*comp_term*sed_get_value(egamm)*decay_term) ;

}

/*-----------Second part of the integral ----------------------*/

double compton_heating2_p1(double egamm,void * params){
    
    egamm = exp(egamm);
    double nHx1;
    nHx1 = exp(nHX1);

    double decay_term,F_e1h1;

    F_e1h1=cross_sec_e1h1(egamm);
    //---eq(14) ref ([3])
    decay_term = exp(-(nHx1*F_e1h1)*drKpc); 

    return  egamm*(4*k_BeV*comp_term*sed_get_value(egamm)*decay_term) ;

}


double compton_heating2_p2(double egamm,void * params){
    
    egamm = exp(egamm);
    double nHx1,nHex2;
    nHx1 = exp(nHX1);
    nHex2 = exp(nHeX2);

    double decay_term,F_e1h1,F_ehe1;

    F_e1h1=cross_sec_e1h1(egamm);
    F_ehe1=cross_sec_ehe1(egamm);
    
    //---eq(14) ref ([3])
    
    decay_term = exp(-(nHx1*F_e1h1+nHex2*F_ehe1)*drKpc); 

    return  egamm*(4*k_BeV*comp_term*sed_get_value(egamm)*decay_term) ;

}

double compton_heating2_p3(double egamm,void * params){
    
    egamm = exp(egamm);
    double nHx13,nHex2;
    nHx13 = exp(nHX13);
    nHex2 = exp(nHeX2);

    double decay_term,F_e1h1,F_ehe1;

    F_e1h1=cross_sec_e1h1(egamm);
    F_ehe1=cross_sec_ehe1(egamm);
    //---eq(14) ref ([3])
    decay_term = exp(-(nHx13*F_e1h1+nHex2*F_ehe1)*drKpc);
    /* nHx13 = (xH + (k/12.5)*xHeII) got from the program
    k=63.983; (B13) / B(16) Ref.... Fukugita */

    return egamm*(4*k_BeV*comp_term*sed_get_value(egamm)*decay_term) ;
}

double compton_heating2_p4(double egamm,void * params){
    egamm = exp(egamm);
    double nHx13,nHex2;
    nHx13 = exp(nHX13);
    nHex2 = exp(nHeX2);

    double decay_term,F_e1h1,F_ehe1;

    F_e1h1=cross_sec_e1h1(egamm);
    F_ehe1=cross_sec_ehe1(egamm);

    //*exp(-egamm/Ec2);

    decay_term = exp(-(nHx13*F_e1h1+nHex2*F_ehe1)*drKpc);
    /* nHx13 = (xH + (k/12.5)*xHeII) got from the program
    k=63.983; (B13) / B(16) Ref.... Fukugita */

    return egamm*(4*k_BeV*comp_term*sed_get_value(egamm)*decay_term) ;

}
