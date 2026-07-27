/***************************************************************
 * libraries
 ***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/***************************************************************
 * SD headers
 ***************************************************************/
#include "constants.h"
#include "allvars.h"
#include "prototype.h"
#include "log.h"
#include "config_defaults.h"


/***************************************************************
 * set up redshift list
 ***************************************************************/
void density_read_file(char *densityFileName){

    /* function to read the over densities from file and
     * store them in the over_densities array. We assume, for now,
     * that the provided density profile has the correct length and
     * resolution and that there is no file header.
     *
     * TODO: add interpolation routine
     *
     * If no file is provided, fill the array with constant values.
     */


    if(DEBUG)log_debug("Entering density_read_file(...)");

    /* local variables */
    FILE *fp;
    int i, densityLineCount = 0;
    double tmp;


    /* fallback case if no density file is provided */
    if (densityFileName[0]== '\0'){
        for (i=0; i<numGridPoints; i++){
            over_densities[i] = OverDensity;
        }
        return;
    }


    /* open file */
    printf(" Trying to read over densities from file \'%s\' \n",densityFileName);
    fp = fopen(densityFileName,"r");

    if(!fp){
        printf("ERROR: Could not read density file %s \n", densityFileName);
        log_error("Could not read density file '%s'. Exiting.", densityFileName);
        log_close();
        exit(1);
    }


    /* Count data lines.
     * NOTE: the original used while(!feof(fp)), which over-counts by one -- it
     * reads one extra time past the last line before feof trips. That inflated
     * densityLineCount, so a correctly-sized profile failed the size check
     * below and density-profile input was effectively unusable. Testing the
     * fscanf return value directly is correct. */
    while( fscanf(fp,"%le",&tmp) == 1 ){
        densityLineCount++;
    }

    if(DEBUG)log_debug("densityLineCount = %d", densityLineCount);


    if (densityLineCount!=numGridPoints){
        printf(" Error: Number of provided over density grid points does not match size of computing grid.\n");
        printf("        Over density grid points =  %d\n", densityLineCount);
        printf("        Computing grid points    =  %d\n", numGridPoints);
        log_error("Density grid points (%d) != computing grid points (%d). Exiting.", densityLineCount, numGridPoints);
        log_close();
        exit(1);
    }



    /* Go back to the beginning */
    rewind(fp);


    for(i=0; i<densityLineCount; i++){
        if( fscanf( fp,"%le", &over_densities[i] ) != 1 ){
            printf(" Error. Could not read from file '%s'. Exiting.\n", densityFileName);
            log_error("Could not read from file '%s'. Exiting.", densityFileName);
            log_close();
            exit(1);
        }
    }

    fclose(fp);


}
