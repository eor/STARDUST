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
/* 

 REFERENCES:
 [1] Dijkstra, M., Haiman, Z. 2004,ApJ, 613, 646-654 
 [2] Fukugita, M & Kawasaki, M. 1994, MNRAS 269, 563
 [3] Shull, J.M., van Steenberg, M. 1985, ApJ, 298, 268-274
 [4] Zaroubi, S & Silk, J. 2005, MNRAS (TBS)
 [5] Zhao D.H., Mo H.J., Jing Y.P., Borner G. 2003, MNRAS,339,12.
*/ 
#ifdef POWERLAWSRC
#define Ec1 2e3      /*** cut-off Ref (Sazonov et al., 2004) ***/
#define Ec2 1e4      /*** cut-off Ref (Zaroubi & Silk, 2005) ***/
#endif

/***************************************************************
 * global variables (scope - this file)
 ***************************************************************/
double estart;


int table_temp(void){
    
    printf(" Creating tables: TEMPERATURE\n");
    if(DEBUG)log_debug("Starting table creation (TEMPERATURE)");    
    
    /***************************************************************
     * local variable declaration
     ***************************************************************/

    FILE *fpTemp;    
    char *tmpPath = NULL;
        
    const char *tableFile_e1h1_p1 = "table_temp_e1h1_p1.dat";
    const char *tableFile_e1h1_p2 = "table_temp_e1h1_p2.dat";
    const char *tableFile_e1h1_p3 = "table_temp_e1h1_p3.dat";
    
    const char *tableFile_ehe1_p2 = "table_temp_ehe1_p2.dat";
    const char *tableFile_ehe1_p3 = "table_temp_ehe1_p3.dat";
    
    const char *tableFile_ehe2_p3 = "table_temp_ehe2_p3.dat";
    
    double l1,l2,l3,l4;
    double Max_Log1,Max_Log2;
    double previous = 100.; /* Some relatively large number */
    
    double absError = 1e-3;     /* integration errors */
    double relError = 1e-3;
    
    l1 = log(HIONIZEeV);
    l2 = log(He1IONIZEeV);
    l3 = log(He2IONIZEeV);
    l4 = log(myConfig.sourceEHigh);

    Max_Log1 = UPLIM1+TABLERES;
    Max_Log2 = UPLIM2+TABLERES;
    
    size_t intResolution = 100;


    /* Setting up the GSL integration environment */

    gsl_integration_workspace *wtable = gsl_integration_workspace_alloc (intResolution);

    double result, error;

    gsl_function Func;
    
    /***************************************************************
     * E1 HI
     ***************************************************************/
    
    tmpPath = utils_concat_path( myConfig.pathOutDir, (char*) tableFile_e1h1_p1 );    
    fpTemp  = fopen( tmpPath , "w");     
    
    Func.function = &temp_e1h1_p1_log;
    Func.params = 0;

    estart  = HIONIZEeV;
       
    for(nHX1=LOWLIM;nHX1<=Max_Log1;nHX1+=TABLERES){
        
        if(previous>1.e-100){
            
            /* Integration */
            gsl_integration_qags (&Func, l1, l2, absError, relError, intResolution, wtable, &result, &error);

            previous = result;
            fprintf(fpTemp,"%e \t %e \n",nHX1,previous);            
        }
        else
            fprintf(fpTemp,"%e \t %e\n",nHX1,1.e-100);
        
    }
    
    fclose(fpTemp);

    if(DEBUG)log_debug(" Computed temperature tables for H in first energy interval (l1-l2)");        

    tmpPath = utils_concat_path( myConfig.pathOutDir, (char*) tableFile_e1h1_p2 );    
    fpTemp  = fopen( tmpPath , "w");     
    
    Func.function = &temp_e1h1_p2_log;
    Func.params = 0;
     
    for(nHX1=LOWLIM;nHX1<=Max_Log1;nHX1+=TABLERES){ 
        
        nHeX2 = LOWLIM;
        gsl_integration_qags (&Func, l2, l3, absError, relError, intResolution, wtable, &result, &error);
    
        if(result > 1.e-100)
            previous = 100.;
        else 
            previous = 1.e-200;
   
        for(nHeX2=LOWLIM;nHeX2<=Max_Log1;nHeX2+=TABLERES){

            if(previous>1.e-100){                
                previous = result;
                fprintf(fpTemp,"%e \t %e \t %e \n",nHX1,nHeX2,previous);
            }      
            else
                fprintf(fpTemp,"%e \t %e \t %e \n",nHX1,nHeX2,1.e-100);
            
        }        
    }
    
    fclose(fpTemp);
    
    if(DEBUG)log_debug(" Computed temperature tables for H in SECOND energy interval (l2-l3)");         

    tmpPath = utils_concat_path( myConfig.pathOutDir, (char*) tableFile_e1h1_p3 );    
    fpTemp  = fopen( tmpPath , "w");     
    
    Func.function = &temp_e1h1_p3_log;
    Func.params = 0;
       
    for(nHX13=LOWLIM;nHX13<=Max_Log2;nHX13+=TABLERES){ 
        
        nHeX2 = LOWLIM;
        
        gsl_integration_qags (&Func, l3, l4, absError, relError, intResolution, wtable, &result, &error);
    
        if(result>1.e-100)
            previous = 100.;    // re-initialize to arbritary high value            
        else
            previous = 1.e-200;
        
        for(nHeX2=LOWLIM;nHeX2<=Max_Log1;nHeX2+=TABLERES){
            
            if(previous>1.e-100){
                previous = result;
                fprintf(fpTemp,"%e \t %e \t %e \n",nHX13,nHeX2,previous);  
            }
            else
                fprintf(fpTemp,"%e \t %e \t %e \n",nHX13,nHeX2,1.e-100);
        }
        
    }
    
    fclose(fpTemp);
    
    if(DEBUG)log_debug(" Computed temperature tables for H in third energy interval (l3-l4)");         
       
    /***************************************************************
     * E Helium I
     ***************************************************************/    
           
    estart  = He1IONIZEeV;       
 
    tmpPath = utils_concat_path( myConfig.pathOutDir, (char*) tableFile_ehe1_p2 );    
    fpTemp  = fopen( tmpPath , "w");     
       
    Func.function = &temp_ehe1_p2_log;
    Func.params = 0;
       
    for(nHX1=LOWLIM;nHX1<=Max_Log1;nHX1+=TABLERES){ 
        
        nHeX2 = LOWLIM;
        gsl_integration_qags (&Func, l2, l3, absError, relError, intResolution,  wtable, &result, &error);

        if(result >1.e-100)
            previous = 100.; 
        else
            previous = 1.e-200;
        
        for(nHeX2=LOWLIM;nHeX2<=Max_Log1;nHeX2+=TABLERES){
            
            if(previous>1.e-100){
                previous = result;
                fprintf(fpTemp,"%e \t %e \t %e \n",nHX1,nHeX2,previous);
            }
            else
                fprintf(fpTemp,"%e \t %e \t %e \n",nHX1,nHeX2,1.e-100);
        }        
    }
    fclose(fpTemp);
          
    if(DEBUG)log_debug(" Computed temperature tables for HeI in second energy interval (l2-l3)");      

    tmpPath = utils_concat_path( myConfig.pathOutDir, (char*) tableFile_ehe1_p3 );
    fpTemp  = fopen( tmpPath , "w");     
    
    Func.function = &temp_ehe1_p3_log;
    Func.params = 0;

    previous = 100.; /* re-initialize to arbritary high value  */
       
    for(nHX13=LOWLIM;nHX13<=Max_Log2;nHX13+=TABLERES){  
        
        nHeX2 = LOWLIM;
        
        gsl_integration_qags (&Func, l3, l4, absError, relError, intResolution, wtable, &result, &error);

        if(result >1.e-100)
            previous = 100.; 
        else
            previous = 1.e-200;
        
        for(nHeX2=LOWLIM;nHeX2<=Max_Log1;nHeX2+=TABLERES){
            
            if(previous>1.e-100){
                previous = result;
                fprintf(fpTemp,"%e \t %e \t %e \n",nHX13,nHeX2,previous);
            }
            else
                fprintf(fpTemp,"%e \t %e \t %e \n",nHX13,nHeX2,1.e-100);
        }
    }
    
    fclose(fpTemp);
       
    if(DEBUG)log_debug(" Computed temperature tables for HeI in third energy interval (l3-l4)");             

    /***************************************************************
     * E Helium II
     ***************************************************************/    
       
    estart  = He2IONIZEeV;
       
    tmpPath = utils_concat_path( myConfig.pathOutDir, (char*) tableFile_ehe2_p3 );
    fpTemp  = fopen( tmpPath , "w");     

    Func.function = &temp_ehe2_p3_log;
    Func.params = 0;
       
    for(nHX13=LOWLIM;nHX13<=Max_Log2;nHX13+=TABLERES){   
        
        nHeX2 = LOWLIM;
        gsl_integration_qags (&Func, l3, l4, absError, relError, intResolution, wtable, &result, &error);
 
        if(result > 1.e-100)
            previous = 100.; 
        else
            previous = 1.e-200;
   
        for(nHeX2=LOWLIM;nHeX2<=Max_Log1;nHeX2+=TABLERES){
            
            if(previous>1.e-100){
                previous = result;
                fprintf(fpTemp,"%e \t %e \t %e \n",nHX13,nHeX2,previous);
            } 
            else
                fprintf(fpTemp,"%e \t %e \t %e \n",nHX13,nHeX2,1.e-100);
        }
    }
    
    fclose(fpTemp);
       
    if(DEBUG)log_debug(" Computed Temperature tables for HeII in third energy interval (l3-l4)");      
    
       
    /* freeing the GSL integration environment */ 
    gsl_integration_workspace_free (wtable);
 

    if(DEBUG)log_debug("Leaving table creation (TEMPERATURE)");
    
    return 1;
}

/*******************************E1H1*************************************/
double temp_e1h1_p1_log(double egamm,void * params)
{
    egamm = exp(egamm);
    double nHx1;
    nHx1 = exp(nHX1);

    double decay_term, F_e1h1;
    F_e1h1=cross_sec_e1h1(egamm);

    decay_term = exp(-(nHx1*F_e1h1)*drKpc);

    return egamm * (egamm-estart) * F_e1h1 * sed_get_value(egamm) * decay_term / egamm;
    /* returning Function to be integrated*/
 
}

/* the multiplication by egamm in the "function to be integrated is
   because the integration is done in the LOG and we need to transform */

double temp_e1h1_p2_log(double egamm,void * params)
{
 egamm = exp(egamm);
 double nHx1,nHex2;
 nHx1 = exp(nHX1);
 nHex2 = exp(nHeX2);
 double decay_term,F_e1h1,F_ehe1;
 F_e1h1=cross_sec_e1h1(egamm);
 F_ehe1=cross_sec_ehe1(egamm);
 
 decay_term = exp(-(nHx1*F_e1h1+nHex2*F_ehe1)*drKpc);

 return egamm*(egamm-estart)*F_e1h1*sed_get_value(egamm)*decay_term/egamm; 
 }

double temp_e1h1_p3_log(double egamm,void * params)
{
 egamm = exp(egamm);
 double nHx13,nHex2;
 nHx13 = exp(nHX13);
 nHex2 = exp(nHeX2);
 double decay_term, F_e1h1,F_ehe1;
 /* float k=63.983; // (B13) / B(16) Ref.... Fukugita */
 F_e1h1=cross_sec_e1h1(egamm);
 F_ehe1=cross_sec_ehe1(egamm);
 
 decay_term = exp(-(nHx13*F_e1h1+nHex2*F_ehe1)*drKpc);
 /* nHx13 = (xH + (k/12.5)*xHeII) got from the program */

 return egamm*(egamm-estart)*F_e1h1*sed_get_value(egamm)*decay_term/egamm;
}

double temp_e1h1_p4_log(double egamm,void * params)
{
 egamm = exp(egamm);
 double nHx13,nHex2;
 nHx13 = exp(nHX13);
 nHex2 = exp(nHeX2);
 double decay_term, F_e1h1,F_ehe1;
 /* float k=63.983; // (B13) / B(16) Ref.... Fukugita */
 F_e1h1=cross_sec_e1h1(egamm);
 F_ehe1=cross_sec_ehe1(egamm);
  //*exp(-egamm/Ec2);

 decay_term = exp(-(nHx13*F_e1h1+nHex2*F_ehe1)*drKpc);

return egamm*(egamm-estart)*F_e1h1*sed_get_value(egamm)*decay_term/egamm;

}
/*******************************EHe1*********************************/


double temp_ehe1_p2_log(double egamm,void * params)
{
 egamm = exp(egamm);
 double nHx1,nHex2;
 nHx1 = exp(nHX1);
 nHex2 = exp(nHeX2);
 double decay_term, F_e1h1,F_ehe1;
 F_e1h1=cross_sec_e1h1(egamm);
 F_ehe1=cross_sec_ehe1(egamm);
 
 decay_term = exp(-(nHx1*F_e1h1+nHex2*F_ehe1)*drKpc);
 
 return egamm*(egamm-estart)*F_ehe1*sed_get_value(egamm)*decay_term/egamm;
}


double temp_ehe1_p3_log(double egamm,void * params)
{ 
 egamm = exp(egamm);
 double nHx13,nHex2;
 nHx13 = exp(nHX13);
 nHex2 = exp(nHeX2);
 double decay_term, F_e1h1, F_ehe1;
 F_e1h1=cross_sec_e1h1(egamm);
 F_ehe1=cross_sec_ehe1(egamm);
 
 decay_term = exp(-(nHx13*F_e1h1+nHex2*F_ehe1)*drKpc);

 return egamm*(egamm-estart)*F_ehe1*sed_get_value(egamm)*decay_term/egamm;
}

double temp_ehe1_p4_log(double egamm,void * params)
{ 
 egamm = exp(egamm);
 double nHx13,nHex2;
 nHx13 = exp(nHX13);
 nHex2 = exp(nHeX2);
 double decay_term, F_e1h1,F_ehe1;
 F_e1h1=cross_sec_e1h1(egamm);
 F_ehe1=cross_sec_ehe1(egamm);
  //*exp(-egamm/Ec2);

 decay_term = exp(-(nHx13*F_e1h1+nHex2*F_ehe1)*drKpc);

 return egamm*(egamm-estart)*F_ehe1*sed_get_value(egamm)*decay_term/egamm;     
}

/***********************************EHe2*******************************/

double temp_ehe2_p3_log(double egamm,void * params)
{
 egamm = exp(egamm);
 double nHx13,nHex2;
 nHx13 = exp(nHX13);
 nHex2 = exp(nHeX2);
 double decay_term,F_e1h1,F_ehe1;
 float k=63.983; // (B13) / B(16) Ref.... Fukugita
 F_e1h1=cross_sec_e1h1(egamm);
 F_ehe1=cross_sec_ehe1(egamm);
 
 decay_term = exp(-(nHx13*F_e1h1+nHex2*F_ehe1)*drKpc);  

 return egamm*(egamm-estart)*(k*F_e1h1)*sed_get_value(egamm)*decay_term/egamm;
}

double temp_ehe2_p4_log(double egamm,void * params)
{
 egamm = exp(egamm);
 double nHx13,nHex2;
 nHx13 = exp(nHX13);
 nHex2 = exp(nHeX2);
 double decay_term, F_e1h1,F_ehe1;
 float k=63.983; // (B13) / B(16) Ref.... Fukugita
 F_e1h1=cross_sec_e1h1(egamm);
 F_ehe1=cross_sec_ehe1(egamm);
                     //*exp(-egamm/Ec2);

 decay_term = exp(-(nHx13*F_e1h1+nHex2*F_ehe1)*drKpc);  

 return egamm*(egamm-estart)*(k*F_e1h1)*sed_get_value(egamm)*decay_term/egamm; 
 /* because: F_ehe2=k*F_e1h1 */
 
}
/*------------------------------------------------------------------*/
