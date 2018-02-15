#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <dirent.h>
#include <libconfig.h>


#include "allvars.h"
#include "prototype.h"
#include "config_defaults.h"
#include "constants.h"

extern struct global_config myConfig;
struct config_t cfg;

void config_load_from_file(char *fileName){

    /* 1. Set default parameter values */
    
    
    /* paths */
    char const *pathOutDirDefault   = CONFIG_DEFAULT_PATH_OUTDIR; 
    char const *pathSEDDefault      = CONFIG_DEFAULT_PATH_SED;
    char const *pathDensityDefault  = CONFIG_DEFAULT_PATH_DENSITY;      // for future use
    char const *pathIDDefault       = CONFIG_DEFAULT_PATH_ID;
    
    /* simulation */
    double  sourceELowDefault       = CONFIG_DEFAULT_SOURCE_ELOW;
    double  sourceEHighDefault      = CONFIG_DEFAULT_SOURCE_EHIGH;     
    double  sourceLifetimeDefault   = CONFIG_DEFAULT_SOURCE_LIFETIME;        
    
    double  haloMass                = CONFIG_DEFAULT_HALOMASS;      

    double  redshiftLowDefault      = CONFIG_DEFAULT_REDSHIFT_LOW;
    double  redshiftHighDefault     = CONFIG_DEFAULT_REDSHIFT_HIGH;
    double  redshiftStrideDefault   = CONFIG_DEFAULT_REDSHIFT_STRIDE;      
    
    /* cosmology  */    
    double  cosmoOmegaMDefault      = CONFIG_DEFAULT_COSMO_OMEGA_M;     //  Total matter density 
    double  cosmoOmegaLDefault      = CONFIG_DEFAULT_COSMO_OMEGA_L;     //  Dark energy / cosmological constant
    double  cosmoOmegaBDefault      = CONFIG_DEFAULT_COSMO_OMEGA_B;     //  Baryon density
    double  cosmoH0Default          = CONFIG_DEFAULT_COSMO_H0;          //  Hubble parameter [km/sec/Mpc] at z = 0
    double  cosmoH100Default        = CONFIG_DEFAULT_COSMO_H100;        //  Hubble parameter / 100
    double  cosmoSigma8Default      = CONFIG_DEFAULT_COSMO_SIGMA8;      //  power spectrum normalization    
    double  cosmoTauThomDefault     = CONFIG_DEFAULT_COSMO_TAU_THOM;    //  Optical depth of Thomson scattering
    double  cosmoTCMBDefault        = CONFIG_DEFAULT_COSMO_TCMB0;       //  CMB temperature at z = 0
    
    /* general setting */
    int     settingsDebugDefault    = CONFIG_DEFAULT_SETTINGS_DEBUG;       
    
    double  settingsRMaxDefault     = CONFIG_DEFAULT_SETTINGS_RMAX;     // kpc    
    double  settingsRStartDefault   = CONFIG_DEFAULT_SETTINGS_RSTART;    
    
    double  settingsDeltaRDefault   = CONFIG_DEFAULT_SETTINGS_DELTA_R;  
    double  settingsDeltaTDefault   = CONFIG_DEFAULT_SETTINGS_DELTA_T;  //  in Myr
    
    double  settingsWriteTDefault   = CONFIG_DEFAULT_SETTINGS_WRITE_T;  // write interval in Myr
    
    
    /* 2. Read config */
  
    /* Initialize the configuration */
    config_init(&cfg); 
  
    /* Load the file */
    printf(" Loading config file \"%s\" ... ", fileName);
    if (!config_read_file(&cfg, fileName)){
        
        printf("\n Error: Line: %d -> %s\n Exiting.\n", config_error_line(&cfg), config_error_text(&cfg));
        config_destroy(&cfg);  
        exit(1);
    }
    else{
    
        printf("ok\n");
        printf(" Loading parameters\n");
        printf(" Using the following values:\n\n");
    
        
       /* path parameters */       
        
       /* pathOutDir */
        const char *pathOutDirTmp;
        if (config_lookup_string(&cfg, "paths.pathOutDir", &pathOutDirTmp)){
            
            /* check if path actually exists */
            DIR* checkDir = opendir( pathOutDirTmp );
            
            if(checkDir){
                closedir(checkDir); 
                strcpy(myConfig.pathOutDir, pathOutDirTmp);
                myConfig.pathOutDir[strlen(pathOutDirTmp)] = '\0'; 
            }
            else if(ENOENT == errno){
                printf("\tWarning: pathOutDir '%s' does not exist.\n", pathOutDirTmp);
                printf("\tFalling back to default.\n\n");
                strcpy(myConfig.pathOutDir, pathOutDirDefault);
            }else{
                printf("\tWarning: Could not access pathOutDir '%s'. Check permissions!\n", pathOutDirTmp);
                printf("\tFalling back to default.\n\n");
                strcpy(myConfig.pathOutDir, pathOutDirDefault);
            }
            

            printf("\tpathOutDir \t = %s\n", myConfig.pathOutDir);
            
            
        }else{
            strcpy(myConfig.pathOutDir, pathOutDirDefault);
            printf("\tpathOutDir \t = %s  (not set, using default value)\n", myConfig.pathOutDir);      
        }
        
       /* pathSED */
        const char *pathSEDTmp;
        if (config_lookup_string(&cfg, "paths.pathSED", &pathSEDTmp)){
            strcpy(myConfig.pathSED, pathSEDTmp);
            myConfig.pathSED[strlen(pathSEDTmp)] = '\0'; 
            printf("\tpathSED \t = %s\n", myConfig.pathSED);
        }else{
            strcpy(myConfig.pathSED, pathSEDDefault);
            printf("\tpathSED \t = %s  (not set, using default value)\n", myConfig.pathSED);      
        }

        /* pathDensity */
        const char *pathDensityTmp;
        if (config_lookup_string(&cfg, "paths.pathDensity", &pathDensityTmp)){
            strcpy(myConfig.pathDensity, pathDensityTmp);
            myConfig.pathDensity[strlen(pathDensityTmp)] = '\0'; 
            printf("\tpathDensity \t = %s\n", myConfig.pathDensity);
        }else{
            strcpy(myConfig.pathDensity, pathDensityDefault);
            printf("\tpathDensity \t = %s  (not set, using default value)\n", myConfig.pathDensity);      
        }       

        /* pathID */
        const char *pathIDTmp;
        if (config_lookup_string(&cfg, "paths.pathID", &pathIDTmp)){
            strcpy(myConfig.pathID, pathIDTmp);
            myConfig.pathID[strlen(pathIDTmp)] = '\0'; 
            printf("\tpathID \t\t = %s\n", myConfig.pathID);
        }else{
            strcpy(myConfig.pathID, pathIDDefault);
            printf("\tpathID \t\t = %s  (not set, using default value)\n", myConfig.pathID);      
        }        
       
       
       
       
        printf("\n");       
        
        /* Simulation parameters */
       
        /* sourceELow */
        if (config_lookup_float(&cfg, "simulation.sourceELow", &myConfig.sourceELow)){
            printf("\tsourceELow \t = %e\n", myConfig.sourceELow);
        }else{
            myConfig.sourceELow = sourceELowDefault;
            printf("\tsourceELow \t = %e (not set, using default value)\n",myConfig.sourceELow);
        }       
       
       /* sourceEHigh */  
        double tmpSourceEHigh = 0.0;
        if (config_lookup_float(&cfg, "simulation.sourceEHigh",&tmpSourceEHigh)){
            
            if (tmpSourceEHigh > He2IONIZEeV){
                myConfig.sourceEHigh = tmpSourceEHigh;
                printf("\tsourceEHigh \t = %e\n", myConfig.sourceEHigh); 
            }else{
                myConfig.sourceEHigh = sourceEHighDefault;
                printf("\tsourceEHigh \t = %e (not set or too small, using default value)\n",myConfig.sourceEHigh);
            
            }
        }   
        
        /* sourceLifetime */   
        if (config_lookup_float(&cfg, "simulation.sourceLifetime", &myConfig.sourceLifetime)){
            printf("\tsourceLifetime \t = %f\n", myConfig.sourceLifetime);
        }else{
            myConfig.sourceLifetime = sourceLifetimeDefault;
            printf("\tsourceLifetime \t = %f (not set, using default value)\n",myConfig.sourceLifetime);
        } 
       
        /* sourceEfficiency */ 
//         if (config_lookup_float(&cfg, "simulation.sourceEfficiency", &myConfig.sourceEfficiency)){
//             printf("\tsourceEfficiency = %f\n", myConfig.sourceEfficiency);
//         }else{
//             myConfig.sourceEfficiency = sourceEfficiencyDefault;
//             printf("\tsourceEfficiency = %f (not set, using default value)\n",myConfig.sourceEfficiency);
//         }       
        
        /* haloMass */      
        if (config_lookup_float(&cfg, "simulation.haloMass", &myConfig.haloMass)){
            printf("\thaloMass \t = %f\n", myConfig.haloMass);
        }else{
            myConfig.haloMass = haloMass;
            printf("\thaloMass \t = %f (not set, using default value)\n",myConfig.haloMass);
        }   
        
        /* haloMassHigh */     
//         if (config_lookup_float(&cfg, "simulation.haloMassHigh", &myConfig.haloMassHigh)){
//             printf("\thaloMassHigh \t = %f\n", myConfig.haloMassHigh);
//         }else{
//             myConfig.haloMassHigh = haloMassHighDefault;
//             printf("\thaloMassHigh \t = %f (not set, using default value)\n",myConfig.haloMassHigh);
//         }        
        
        /* haloMassStride */   
//         if (config_lookup_float(&cfg, "simulation.haloMassStride", &myConfig.haloMassStride)){
//             printf("\thaloMassStride \t = %f\n", myConfig.haloMassStride);
//         }else{
//             myConfig.haloMassStride = haloMassStrideDefault;
//             printf("\thaloMassStride \t = %f (not set, using default value)\n",myConfig.haloMassStride);
//         }
        
        /* redshiftLow */      
        if (config_lookup_float(&cfg, "simulation.redshiftLow", &myConfig.redshiftLow)){
            printf("\tredshiftLow \t = %f\n", myConfig.redshiftLow);
        }else{
            myConfig.redshiftLow = redshiftLowDefault;
            printf("\tredshiftLow \t = %f (not set, using default value)\n",myConfig.redshiftLow);
        }
        
        /* redshiftHigh */     
        if (config_lookup_float(&cfg, "simulation.redshiftHigh", &myConfig.redshiftHigh)){
            printf("\tredshiftHigh \t = %f\n", myConfig.redshiftHigh);
        }else{
            myConfig.redshiftHigh = redshiftHighDefault;
            printf("\tredshiftHigh \t = %f (not set, using default value)\n",myConfig.redshiftHigh);
        }
        
        /* redshiftStride */   
        if (config_lookup_float(&cfg, "simulation.redshiftStride", &myConfig.redshiftStride)){
            printf("\tredshiftStride \t = %f\n", myConfig.redshiftStride);
        }else{
            myConfig.redshiftStride = redshiftStrideDefault;
            printf("\tredshiftStride \t = %f (not set, using default value)\n",myConfig.redshiftStride);
        }       
        
        
        printf("\n");       
        
        /* Cosmological parameters */
       
        /* cosmoOmegaM */    
        if (config_lookup_float(&cfg, "cosmology.cosmoOmegaM", &myConfig.cosmoOmegaM)){
            printf("\tcosmoOmegaM \t = %f\n", myConfig.cosmoOmegaM);
        }else{
            myConfig.cosmoOmegaM = cosmoOmegaMDefault;
            printf("\tcosmoOmegaM \t = %f (not set, using default value)\n",myConfig.cosmoOmegaM);
        } 
        /* cosmoOmegaL */   
        if (config_lookup_float(&cfg, "cosmology.cosmoOmegaL", &myConfig.cosmoOmegaL)){
            printf("\tcosmoOmegaL \t = %f\n", myConfig.cosmoOmegaL);
        }else{
            myConfig.cosmoOmegaL = cosmoOmegaLDefault;
            printf("\tcosmoOmegaL \t = %f (not set, using default value)\n",myConfig.cosmoOmegaL);
        }         
        /* cosmoOmegaB */ 
        if (config_lookup_float(&cfg, "cosmology.cosmoOmegaB", &myConfig.cosmoOmegaB)){
            printf("\tcosmoOmegaB \t = %f\n", myConfig.cosmoOmegaB);
        }else{
            myConfig.cosmoOmegaB = cosmoOmegaBDefault;
            printf("\tcosmoOmegaB \t = %f (not set, using default value)\n",myConfig.cosmoOmegaB);
        }         
        /* cosmoH0 */    
        if (config_lookup_float(&cfg, "cosmology.cosmoH0", &myConfig.cosmoH0)){
            printf("\tcosmoH0 \t = %f\n", myConfig.cosmoH0);
        }else{
            myConfig.cosmoH0 = cosmoH0Default;
            printf("\tcosmoH0 \t = %f (not set, using default value)\n",myConfig.cosmoH0);
        }       
        
        /* cosmoH100 */    
        if (config_lookup_float(&cfg, "cosmology.cosmoH100", &myConfig.cosmoH100)){
            printf("\tcosmoH100 \t = %f\n", myConfig.cosmoH100);
        }else{
            myConfig.cosmoH100 = cosmoH100Default;
            printf("\tcosmoH100 \t = %f (not set, using default value)\n",myConfig.cosmoH100);
        }
        
        /* cosmoSigma8 */   
        if (config_lookup_float(&cfg, "cosmology.cosmoSigma8", &myConfig.cosmoSigma8)){
            printf("\tcosmoSigma8 \t = %f\n", myConfig.cosmoSigma8);
        }else{
            myConfig.cosmoSigma8 = cosmoSigma8Default;
            printf("\tcosmoSigma8 \t = %f (not set, using default value)\n",myConfig.cosmoSigma8);
        }
        
        /* cosmoTauThom */   
        if (config_lookup_float(&cfg, "cosmology.cosmoTauThom", &myConfig.cosmoTauThom)){
            printf("\tcosmoTauThom \t = %f\n", myConfig.cosmoTauThom);
        }else{
            myConfig.cosmoTauThom = cosmoTauThomDefault;
            printf("\tcosmoTauThom \t = %f (not set, using default value)\n",myConfig.cosmoTauThom);
        }
        
        /* cosmoTCMB */   
        if (config_lookup_float(&cfg, "cosmology.cosmoTCMB", &myConfig.cosmoTCMB)){
            printf("\tcosmoTCMB \t = %f\n", myConfig.cosmoTCMB);
        }else{
            myConfig.cosmoTCMB = cosmoTCMBDefault;
            printf("\tcosmoTCMB \t = %f (not set, using default value)\n",myConfig.cosmoTCMB);
        }
        
        
        printf("\n");       
        
        /* General parameters */
        
        /* settingsDebug */         
        int settingsDebugTmp = 0; 
        //long settingsDebugTmp = 0; // needed to get rid of old libconfig warnings (int / long int conversion)
        if (config_lookup_int(&cfg, "settings.settingsDebug",&settingsDebugTmp)){
                        
            
            if(settingsDebugTmp!=0 && settingsDebugTmp!=1){
                myConfig.settingsDebug = settingsDebugDefault;
                printf("\tsettingsDebug \t = %d (entered value is out of range, falling back on default value)\n",myConfig.settingsDebug);
                
            }else{
                //myConfig.settingsDebug = (int)settingsDebugTmp; // see above
                myConfig.settingsDebug = settingsDebugTmp;
                printf("\tsettingsDebug \t = %d\n",myConfig.settingsDebug); 
                
            }
            
        }else{
            myConfig.settingsDebug = settingsDebugDefault;
            printf("\tsettingsDebug \t = %d (not set, using default value)\n",myConfig.settingsDebug);
        }       
        
        // now we can set the short, handier, global: 
        DEBUG = myConfig.settingsDebug;
        
//      the user should not touch this        
//         /* settingsNEQS */  
//         int settingsNEQSTmp = 0; 
//         if (config_lookup_int(&cfg, "settings.settingsNEQS",&settingsNEQSTmp)){
//             myConfig.settingsNEQS = settingsNEQSTmp;
//             printf("\tsettingsNEQS \t = %d\n",myConfig.settingsNEQS);    
//         }else{
//             myConfig.settingsNEQS = settingsNEQSDefault;
//             printf("\tsettingsNEQS \t = %d (not set, using default value)\n",myConfig.settingsNEQS);
//         }  
        
        /* settingsRMax */  
        if (config_lookup_float(&cfg, "settings.settingsRMax", &myConfig.settingsRMax)){
            printf("\tsettingsRMax \t = %f\n", myConfig.settingsRMax);
        }else{
            myConfig.settingsRMax = settingsRMaxDefault;
            printf("\tsettingsRMax \t = %f (not set, using default value)\n",myConfig.settingsRMax);
        }
        
        /* settingsRStart */
        if (config_lookup_float(&cfg, "settings.settingsRStart", &myConfig.settingsRStart)){
            printf("\tsettingsRStart \t = %f\n", myConfig.settingsRStart);
        }else{
            myConfig.settingsRStart = settingsRStartDefault;
            printf("\tsettingsRStart \t = %f (not set, using default value)\n",myConfig.settingsRStart);
        }
        
        /* settingsDeltaR */
        if (config_lookup_float(&cfg, "settings.settingsDeltaR", &myConfig.settingsDeltaR)){
            printf("\tsettingsDeltaR \t = %f\n", myConfig.settingsDeltaR);
        }else{
            myConfig.settingsDeltaR = settingsDeltaRDefault;
            printf("\tsettingsDeltaR \t = %f (not set, using default value)\n",myConfig.settingsDeltaR);
        }
        
        /* settingsDeltaT */     
        if (config_lookup_float(&cfg, "settings.settingsDeltaT", &myConfig.settingsDeltaT)){            
            printf("\tsettingsDeltaT \t = %f\n", myConfig.settingsDeltaT);            
        }else{
            myConfig.settingsDeltaT = settingsDeltaTDefault;
            printf("\tsettingsDeltaT \t = %f (not set, using default value)\n",myConfig.settingsDeltaT);            
        }  
        
        /* settingsWriteT */     
        if (config_lookup_float(&cfg, "settings.settingsWriteT", &myConfig.settingsWriteT)){            
            printf("\tsettingsWriteT \t = %f\n", myConfig.settingsWriteT);
        }else{
            myConfig.settingsWriteT = settingsWriteTDefault;
            printf("\tsettingsWriteT \t = %f (not set, using default value)\n",myConfig.settingsWriteT);
        }          
        
        printf("\n");          
      
    
  } /* end of else - Load file success */

  /* Free the configuration */
  config_destroy(&cfg);

}
