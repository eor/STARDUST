
/***************************************************************
 * libraries
 ***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/***************************************************************
 * SD headers
 ***************************************************************/
#include "prototype.h"
#include "allvars.h"
#include "constants.h"
#include "log.h"

/***************************************************************
 * routines for  opening and closing our main & debug logs
 ***************************************************************/
void log_start(){
     
    char *tmpPath = NULL;         
         
    if(DEBUG){

        if( !asprintf(&tmpPath, "%s/%s_M%.3f_z%.3f_%s", myConfig.pathOutDir, myConfig.pathID, myConfig.haloMass, myConfig.redshiftLow, debugLogFile ) ){

            log_error("asprintf failed");
            log_close();
            exit(1);             
        }
                
        fpDebugLog = fopen(tmpPath, "w");
        log_debug("Starting debug log.");
        log_time(1);
    }

    tmpPath = NULL;
    
    if( !asprintf(&tmpPath, "%s/%s_M%.3f_z%.3f_%s", myConfig.pathOutDir, myConfig.pathID, myConfig.haloMass, myConfig.redshiftLow, mainLogFile ) ){
        log_error("asprintf failed");
        log_close();
        exit(1);             
    }
    
    fpMainLog = fopen(tmpPath, "w");
    log_main("Starting main log");
    log_main("Welcome to STARDUST V%.3f", SD_VERSION);
    log_time(0);    
    
    // dump parameters
    log_main("Using the following parameters:");
    log_main("-----------------------------------------------------------------" );
    log_main("\t pathOutDir   \t= %s",   myConfig.pathOutDir     );
    log_main("\t pathSED      \t= %s",   myConfig.pathSED        );
    log_main("\t pathDensity  \t= %s",   myConfig.pathDensity    );
    log_main("\t pathID       \t= %s",   myConfig.pathID         );
        
    log_main("\t sourceELow   \t= %e",   myConfig.sourceELow     );
    log_main("\t sourceEHigh  \t= %e",   myConfig.sourceEHigh    );
    log_main("\t sourceLifetime = %.3f", myConfig.sourceLifetime ); 
    
    log_main("\t haloMass     \t= %e",   myConfig.haloMass       );
    log_main("\t redshiftLow  \t= %.3f", myConfig.redshiftLow    );
    log_main("\t redshiftHigh \t= %.3f", myConfig.redshiftHigh   ); 
    log_main("\t redshiftStride = %.3f", myConfig.redshiftStride ); 
    
    
    log_main("\t cosmoOmegaM  \t= %.3f", myConfig.cosmoOmegaM    );
    log_main("\t cosmoOmegaL  \t= %.3f", myConfig.cosmoOmegaL    );
    log_main("\t cosmoOmegaB  \t= %.3f", myConfig.cosmoOmegaB    ); 
    log_main("\t cosmoH0      \t= %.3f", myConfig.cosmoH0        ); 
    log_main("\t cosmoH100    \t= %.3f", myConfig.cosmoH100      );
    log_main("\t cosmoSigma8  \t= %.3f", myConfig.cosmoSigma8    );
    log_main("\t cosmoTauThom \t= %.3f", myConfig.cosmoTauThom   ); 
    log_main("\t cosmoTCMB    \t= %.3f", myConfig.cosmoTCMB      ); 
    
   
    log_main("\t settingsDebug\t= %d",   myConfig.settingsDebug  ); 
    log_main("\t settingsRMax \t= %.3f", myConfig.settingsRMax   ); 
    log_main("\t settingsRStart = %.3f", myConfig.settingsRStart );
    log_main("\t settingsDeltaR = %.3f", myConfig.settingsDeltaR );
    log_main("\t settingsDeltaT = %.3f", myConfig.settingsDeltaT ); 
    log_main("\t settingsWriteT = %.3f", myConfig.settingsWriteT );    
    
    log_main("-----------------------------------------------------------------" ); 
    
    free(tmpPath);    
}


void log_close(){
 
    if(DEBUG){
        log_time(1);
        log_debug("Exiting here.");
        fclose(fpDebugLog);
    }
    log_time(0);
    log_main("Exiting here.");
    fclose(fpMainLog);
    
}


void log_time(int type){
    
    time_t timeNow;
    struct tm * timeInfo;

    time(&timeNow);
    
    timeInfo = localtime(&timeNow);
    
    if (type==1){
        // type == 1 is the debug case
        log_debug("Time stamp: %s", asctime (timeInfo) );
    }   
    else{
        // default to main log
        log_main("Time stamp: %s", asctime (timeInfo) );
    }
}
