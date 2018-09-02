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

#ifdef POWERLAWSRC
#define Ec1 2e3   /**** cut-off Ref (Sazonov et al., 2004) ****/
#define Ec2 1e4      
#endif



// TODO: add adaptive relError handling to the integration functions, just like we did for the SED integration


int table_ion (void){       

    printf(" Creating tables: IONIZATION\n");
    if(DEBUG)log_debug("Starting table creation (IONIZATION)");

    /***************************************************************
     * local variables 
     ***************************************************************/      
    FILE *fpTau;
    char *tmpPath = NULL;
        
    const char *tableFile_e1h1_p1 = "table_ion_e1h1_p1.dat";
    const char *tableFile_e1h1_p2 = "table_ion_e1h1_p2.dat";
    const char *tableFile_e1h1_p3 = "table_ion_e1h1_p3.dat";
    
    const char *tableFile_ehe1_p2 = "table_ion_ehe1_p2.dat";
    const char *tableFile_ehe1_p3 = "table_ion_ehe1_p3.dat";
    
    const char *tableFile_ehe2_p3 = "table_ion_ehe2_p3.dat";
    
    double l1, l2, l3, l4;
    double Max_Log1, Max_Log2;
    double previous = 100.;     /* some relatively large number */
    
    double absError = 1e-3;     /* integration errors */
    double relError = 1e-2;
    
    size_t intResolution = 1000;

    l1 = log(HIONIZEeV);
    l2 = log(He1IONIZEeV);
    l3 = log(He2IONIZEeV);
    l4 = log(myConfig.sourceEHigh);    
    
    Max_Log1 = UPLIM1+TABLERES;
    Max_Log2 = UPLIM2+TABLERES;

    /* Setting up the GSL integration environment */
    gsl_integration_workspace       *wtable  = gsl_integration_workspace_alloc (intResolution);
    gsl_integration_cquad_workspace *wtable2 =  gsl_integration_cquad_workspace_alloc(intResolution);
    
    double result, error;
    gsl_function myFunction; 
  
    /***************************************************************
     * E1 Hydrogen
     ***************************************************************/
    
    tmpPath = utils_concat_path( myConfig.pathOutDir, (char*) tableFile_e1h1_p1 );    
    fpTau   = fopen( tmpPath , "w");     
    
    myFunction.function = &tau_e1h1_p1_log;
    myFunction.params = 0;

//     size_t nevals =1000;
    size_t nevals = intResolution;
    
    for( nHX1=LOWLIM; nHX1<=Max_Log1; nHX1+=TABLERES){
        
        if(previous>1.e-100){
            /* Integration */
//             gsl_integration_qags(&myFunction, l1, l2, absError, relError, 1000, wtable, &result, &error);
            gsl_integration_cquad(&myFunction, l1, l2, absError, relError, wtable2, &result, &error, &intResolution); 
            
            //result = Integrate_using_Patterson_adaptive(l1, l2, absError, relError, &tau_e1h1_p1_log, myFunction.params);
            previous = result;
            fprintf(fpTau, "%e \t %e \n", nHX1, previous);
        }else{            
            fprintf(fpTau, "%e \t %e\n", nHX1, 1.e-100);            
        }      
    }
    
    fclose(fpTau);
  
    if(DEBUG)log_debug(" Computed ionization (tau) table for HI in first energy interval (l1-l2)");

    tmpPath = utils_concat_path( myConfig.pathOutDir, (char*) tableFile_e1h1_p2 );
    fpTau  = fopen( tmpPath , "w");     
    
    
    myFunction.function = &tau_e1h1_p2_log;
    myFunction.params = 0;
     
    for(nHX1=LOWLIM;nHX1<=Max_Log1;nHX1+=TABLERES){ 
        
        nHeX2 = LOWLIM;
//         gsl_integration_qags(&myFunction, l2, l3, absError, relError, intResolution, wtable, &result, &error);
        gsl_integration_cquad(&myFunction, l2, l3, absError, relError, wtable2, &result, &error, &intResolution); 
        
        if(result > 1.e-100)
            previous = 100.;
        else 
            previous = 1.e-200;        
        
       for(nHeX2=LOWLIM;nHeX2<=Max_Log1;nHeX2+=TABLERES){
           
           if(previous>1.e-100){
                previous = result;
                fprintf(fpTau,"%e \t %e \t %e \n",nHX1,nHeX2,previous);
            }      
            else
                fprintf(fpTau,"%e \t %e \t %e \n",nHX1,nHeX2,1.e-100);
           
        }        
    }

    fclose(fpTau);
         
    if(DEBUG)log_debug(" Computed ionization (tau) table for HI in second energy interval (l2-l3)");

    tmpPath = utils_concat_path( myConfig.pathOutDir, (char*) tableFile_e1h1_p3 );
    fpTau  = fopen( tmpPath , "w");    
    
    myFunction.function = &tau_e1h1_p3_log;
    myFunction.params = 0;

    previous=100; 
    
    for(nHX13=LOWLIM;nHX13<=Max_Log2;nHX13+=TABLERES){
        
        nHeX2 = LOWLIM;
         gsl_integration_qags(&myFunction, l3, l4, absError, relError, intResolution, wtable, &result, &error);
        //gsl_integration_cquad(&myFunction, l3, l4, absError, relError, wtable2, &result, &error, &intResolution); 
        
        if(result>1.e-100)
            previous = 100.; // re-initialize to arbritary high value
        else
            previous = 1.e-200;        
        
        for(nHeX2=LOWLIM;nHeX2<=Max_Log1;nHeX2+=TABLERES){
            
            if(previous>1.e-100){
                previous = result;
                fprintf(fpTau,"%e \t %e \t %e \n",nHX13,nHeX2,previous);  
            }
            else
                fprintf(fpTau,"%e \t %e \t %e \n",nHX13,nHeX2,1.e-100);
        }
    }
    
    fclose(fpTau);
       
    if(DEBUG)log_debug(" Computed ionization (tau) table for HI in third energy interval (l3-l4)");  
	
     /***************************************************************
      * E Helium I
      ***************************************************************/    
  
    tmpPath = utils_concat_path( myConfig.pathOutDir, (char*) tableFile_ehe1_p2 );
    fpTau  = fopen( tmpPath , "w");   
    
    myFunction.function = &tau_ehe1_p2_log;
    myFunction.params = 0;

    previous=100;  

    for(nHX1=LOWLIM;nHX1<=Max_Log1;nHX1+=TABLERES){ 
        
        nHeX2 = LOWLIM;
//         gsl_integration_qags (&myFunction, l2, l3, absError, relError, intResolution, wtable, &result, &error);
        gsl_integration_cquad(&myFunction, l2, l3, absError, relError, wtable2, &result, &error, &intResolution); 
        

        if(result > 1.e-100)
            previous = 100.;
        else 
            previous = 1.e-200;
      
        for(nHeX2=LOWLIM;nHeX2<=Max_Log1;nHeX2+=TABLERES){
            
            if(previous>1.e-100){
                previous = result;
                fprintf(fpTau,"%e \t %e \t %e \n",nHX1,nHeX2,previous);
            }      
            else
                fprintf(fpTau,"%e \t %e \t %e \n",nHX1,nHeX2,1.e-100);
        }
    }
    
    fclose(fpTau);

    if(DEBUG)log_debug(" Computed ionization (tau) tables for H, He in second energy interval (l2-l3)");
    
    tmpPath = utils_concat_path( myConfig.pathOutDir, (char*) tableFile_ehe1_p3 );
    fpTau  = fopen( tmpPath , "w");   
    
    myFunction.function = &tau_ehe1_p3_log;
    myFunction.params = 0;

    previous=100; 

    for(nHX13=LOWLIM;nHX13<=Max_Log2;nHX13+=TABLERES){ 
        
        nHeX2 = LOWLIM;
        gsl_integration_qags(&myFunction, l3, l4, absError, relError, intResolution, wtable, &result, &error);

        if(result > 1.e-100)
            previous = 100.; // re-initialize to arbritary high value
        else
            previous = 1.e-200;
        
        for(nHeX2=LOWLIM;nHeX2<=Max_Log1;nHeX2+=TABLERES){
            
            
            if(previous>1.e-100){                
                previous = result;
                
                fprintf(fpTau,"%e \t %e \t %e \n", nHX13, nHeX2, previous);  
            }
            else
                fprintf(fpTau,"%e \t %e \t %e \n", nHX13, nHeX2, 1.e-100);
        }
    }
    
    fclose(fpTau);
  
    if(DEBUG)log_debug(" Computed ionization (tau) table for H, He in third energy interval (l3-l4)");

    /***************************************************************
     * E Helium II
     ***************************************************************/

    tmpPath = utils_concat_path( myConfig.pathOutDir, (char*) tableFile_ehe2_p3 );
    fpTau  = fopen( tmpPath , "w");   
    
    myFunction.function = &tau_ehe2_p3_log;
    myFunction.params   = 0;
    previous            = 100; 

    for(nHX13=LOWLIM;nHX13<=Max_Log2;nHX13+=TABLERES){ 
        
        nHeX2 = LOWLIM;
      
        gsl_integration_qags (&myFunction, l3, l4, absError, relError, intResolution,  wtable, &result, &error);
 
        if(result > 1.e-100)
            previous = 100.; // re-initialize to arbritary high value
        else
            previous = 1.e-200;
        
        for(nHeX2=LOWLIM;nHeX2<=Max_Log1;nHeX2+=TABLERES){
            
            if(previous > 1.e-100){
                previous = result;
                fprintf(fpTau,"%e \t %e \t %e \n", nHX13, nHeX2, previous);  
            }
            else
                fprintf(fpTau,"%e \t %e \t %e \n", nHX13, nHeX2, 1.e-100);            
        }
    }
    
    fclose(fpTau);  

    if(DEBUG)log_debug(" Computed ionization (tau) table for H, He in third energy interval (l3-l4)");
  
    /* freeing the GSL integration environment */
    gsl_integration_workspace_free (wtable);
    gsl_integration_cquad_workspace_free(wtable2);
  
    if(DEBUG)log_debug("Leaving table creation (IONIZATION)");  
    
    return 1; 
}


/***************************************************************
 * NOTE: Integration with the functions below is done using 
 * the logarithmic integration trick ... 
 * hence the extra "energy term = egamm" in the equations.
 ***************************************************************/



/*******************************E1H1*************************************/
double tau_e1h1_p1_log(double egamm, void * params){    
    
    double nHx1;
    double decay_term,F_e1h1;
    
    egamm   = exp(egamm);
    nHx1    = exp(nHX1);
    F_e1h1  = cross_sec_e1h1(egamm);
    
    decay_term = exp(-(nHx1*F_e1h1)*drKpc);

    return egamm*F_e1h1*sed_get_value(egamm)*decay_term/egamm;       /* returning function to be integrated */

}

double tau_e1h1_p2_log(double egamm,void * params){
    
    double decay_term,F_e1h1,F_ehe1;  
    double nHx1,nHex2;
    
    egamm   = exp(egamm);
    nHx1    = exp(nHX1);
    nHex2   = exp(nHeX2);

    F_e1h1  = cross_sec_e1h1(egamm);
    F_ehe1  = cross_sec_ehe1(egamm);
    decay_term = exp(-(nHx1*F_e1h1+nHex2*F_ehe1)*drKpc);

    return egamm*F_e1h1*sed_get_value(egamm)*decay_term/egamm; 
}

double tau_e1h1_p3_log(double egamm,void * params){
        
    double nHx13,nHex2;
    double decay_term, F_e1h1,F_ehe1;
    
    egamm = exp(egamm);
    nHx13 = exp(nHX13);
    nHex2 = exp(nHeX2);
    
    F_e1h1=cross_sec_e1h1(egamm);
    F_ehe1=cross_sec_ehe1(egamm);
    decay_term = exp(-(nHx13*F_e1h1+nHex2*F_ehe1)*drKpc);

    return egamm*F_e1h1*sed_get_value(egamm)*decay_term/egamm; 
}

double tau_e1h1_p4_log(double egamm,void * params){
    
    
    double nHx13,nHex2;
    double decay_term, F_e1h1,F_ehe1;
    
    egamm = exp(egamm);
    nHx13 = exp(nHX13);
    nHex2 = exp(nHeX2);
    
    F_e1h1=cross_sec_e1h1(egamm);
    F_ehe1=cross_sec_ehe1(egamm);

    decay_term = exp(-(nHx13*F_e1h1+nHex2*F_ehe1)*drKpc);
    /* nHx13 = (xH + (k/12.5)*xHeII) got from the program */

    return egamm*F_e1h1*sed_get_value(egamm)*decay_term/egamm;
}
/*******************************EHe1*********************************/


double tau_ehe1_p2_log(double egamm,void * params)
{
 egamm = exp(egamm);
 double nHx1,nHex2;
 nHx1 = exp(nHX1);
 nHex2 = exp(nHeX2);
 double decay_term, F_e1h1,F_ehe1;
 F_e1h1=cross_sec_e1h1(egamm);
 F_ehe1=cross_sec_ehe1(egamm);

  decay_term = exp(-(nHx1*F_e1h1+nHex2*F_ehe1)*drKpc);

 return egamm*F_ehe1*sed_get_value(egamm)*decay_term/egamm;
}
 

double tau_ehe1_p3_log(double egamm,void * params)
{ 
 egamm = exp(egamm);
 double nHx13,nHex2;
 nHx13 = exp(nHX13);
 nHex2 = exp(nHeX2);
 double decay_term, F_e1h1,F_ehe1;
 F_e1h1=cross_sec_e1h1(egamm);
 F_ehe1=cross_sec_ehe1(egamm);
 decay_term = exp(-(nHx13*F_e1h1+nHex2*F_ehe1)*drKpc);

 return egamm*F_ehe1*sed_get_value(egamm)*decay_term/egamm;    
}


double tau_ehe1_p4_log(double egamm,void * params){
    
    egamm = exp(egamm);
    double nHx13,nHex2;
    nHx13 = exp(nHX13);
    nHex2 = exp(nHeX2);
    double decay_term,F_e1h1,F_ehe1;
    F_e1h1=cross_sec_e1h1(egamm);
    F_ehe1=cross_sec_ehe1(egamm);

    decay_term = exp(-(nHx13*F_e1h1+nHex2*F_ehe1)*drKpc);

    return egamm*F_ehe1*sed_get_value(egamm)*decay_term/egamm; 
}
/***********************************EHe2*******************************/

double tau_ehe2_p3_log(double egamm,void * params){

    double nHx13,nHex2;
    double decay_term, F_e1h1,F_ehe1;
    float  k=63.983; // (B13) / B(16) Ref.... Fukugita
    
    egamm   = exp(egamm);
    nHx13   = exp(nHX13);
    nHex2   = exp(nHeX2);
    
    F_e1h1  = cross_sec_e1h1(egamm);
    F_ehe1  = cross_sec_ehe1(egamm);
    
    decay_term = exp(-(nHx13*F_e1h1+nHex2*F_ehe1)*drKpc);  

    return egamm*(k*F_ehe1*sed_get_value(egamm)*decay_term/egamm);
}

double tau_ehe2_p4_log(double egamm,void * params){

    double  decay_term, F_e1h1,F_ehe1;
    float   k = 63.983; // (B13) / B(16) Ref.... Fukugita    
    double  nHx13,nHex2;
    
    egamm   = exp(egamm);
    nHx13   = exp(nHX13);
    nHex2   = exp(nHeX2);

    F_e1h1  = cross_sec_e1h1(egamm);
    F_ehe1  = cross_sec_ehe1(egamm);

    decay_term = exp(-(nHx13*F_e1h1+nHex2*F_ehe1)*drKpc);  

    return egamm*(k*F_ehe1*sed_get_value(egamm)*decay_term/egamm);

}
/*------------------------------------------------------------------*/
