/***************************************************************
 * libraries
 ***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/***************************************************************
 * GSL header files
 ***************************************************************/
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>

/***************************************************************
 * SD headers
 ***************************************************************/
#include "constants.h"
#include "allvars.h"
#include "prototype.h"
#include "table_settings.h"
#include "log.h"

/***************************************************************
 * Etc
 ***************************************************************/
#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))

/***************************************************************
 * global variables needed only in the scope of this file
 ***************************************************************/
int SEDLineCount    = 0;
double SEDNorm      = 0;

/***************************************************************
 * Module-level cached interpolation objects.
 *
 * sed_get_value() and sed_get_log_value() are called from inside
 * the GSL integration loops in table_ion.c, table_temp.c and
 * table_compton.c - potentially hundreds of thousands of times per
 * run. The original code allocated, initialised and freed a fresh
 * gsl_spline + gsl_interp_accel on every single call. Instead we
 * build them once (sed_init_spline, from sed_read_file) and reuse
 * them, freeing at shutdown (sed_free_spline, from memory_free_all).
 *
 * Two splines are kept because the two accessors deliberately use
 * different interpolation types (linear vs. cspline), matching the
 * original per-call behaviour exactly.
 ***************************************************************/
static gsl_interp_accel *sed_accel_lin = NULL;   /* for sed_get_value     */
static gsl_spline       *sed_spline_lin = NULL;  /* linear interpolation  */
static gsl_interp_accel *sed_accel_log = NULL;   /* for sed_get_log_value */
static gsl_spline       *sed_spline_log = NULL;  /* cspline interpolation */

static void sed_init_spline(void){

    sed_accel_lin  = gsl_interp_accel_alloc();
    sed_spline_lin = gsl_spline_alloc(gsl_interp_linear, SEDLineCount);
    gsl_spline_init(sed_spline_lin, photonEnergy, sedLuminosity, SEDLineCount);

    sed_accel_log  = gsl_interp_accel_alloc();
    sed_spline_log = gsl_spline_alloc(gsl_interp_cspline, SEDLineCount);
    gsl_spline_init(sed_spline_log, photonEnergy, sedLuminosity, SEDLineCount);
}

void sed_free_spline(void){

    if(sed_spline_lin){ gsl_spline_free(sed_spline_lin);       sed_spline_lin = NULL; }
    if(sed_accel_lin){  gsl_interp_accel_free(sed_accel_lin);  sed_accel_lin  = NULL; }
    if(sed_spline_log){ gsl_spline_free(sed_spline_log);       sed_spline_log = NULL; }
    if(sed_accel_log){  gsl_interp_accel_free(sed_accel_log);  sed_accel_log  = NULL; }
}


/***************************************************************
 * for a given energy return (interpolated) SED value
 ***************************************************************/
double sed_get_value(double E){

    if (myConfig.settingsStroemgrenTest){
        /* Strömgren test: top-hat ("delta function") source over
         * [stroemgrenPeakE, stroemgrenPeakE + stroemgrenWidth) eV. */
        if (E >= myConfig.stroemgrenPeakE && E < myConfig.stroemgrenPeakE + myConfig.stroemgrenWidth)
            return SEDNorm * 1.0;
        else
            return 1e-100;
    }

    /* interpolate the spectral luminosity (sedLuminosity) at photon energy E (photonEnergy) */
    return gsl_spline_eval(sed_spline_lin, E, sed_accel_lin);
}

/***************************************************************
 * SED log interpolation
 ***************************************************************/
double sed_get_log_value(double E, void *pars){

    /***************************************************************
     * For a given energy (passed in log) this function returns an
     * interpolated SED value.
     ***************************************************************/

    E = exp(E);
    return gsl_spline_eval(sed_spline_log, E, sed_accel_log);
}


/***************************************************************
 * read spectrum of the radiating source
 ***************************************************************/
void sed_read_file(const char *SEDFileName){

    /***************************************************************
     * This function reads the provided SED file.
     * It assumes that the first line is a header, i.e.
     * it will be ignored.
     *
     * File format (whitespace separated, one entry per line):
     *   photon energy [eV]   spectral luminosity [eV/s/eV]
     ***************************************************************/

    if(DEBUG)log_debug("Entering sed_read_file(...)");

    FILE *fp;
    char buf[256];
    double lam, egy;
    int i;

    /* open file */
    printf(" Trying to read SED from file \'%s\' \n",SEDFileName);
    fp = fopen(SEDFileName,"r");

    if(!fp){
        printf("ERROR: Could not read SED file %s \n", SEDFileName);
        log_error("Could not read SED file '%s'. Exiting.", SEDFileName);
        log_close();
        exit(1);
    }

    /* Reading Header */
    if( !fgets(buf, sizeof buf, fp) ){
        printf(" Error. Could not read header from file '%s'. Exiting.\n", SEDFileName);
        log_error("Could not read header from file '%s'. Exiting.", SEDFileName);
        log_close();
        exit(1);
    }

    if(DEBUG) log_debug("File header: %s", buf);

    /* Count data lines.
     * NOTE: the original used while(!feof(fp)), which is a well known
     * antipattern - it reads one extra (garbage) line past the last real
     * one and over-counts SEDLineCount by one. Testing the fscanf return
     * value directly is correct. */
    while( fscanf(fp,"%le %le", &lam, &egy) == 2 ){
        SEDLineCount++;
    }

    if(DEBUG)log_debug("SEDLineCount = %d", SEDLineCount);

    if(SEDLineCount < 2){
        printf(" Error. SED file '%s' contains fewer than 2 data points. Exiting.\n", SEDFileName);
        log_error("SED file '%s' contains fewer than 2 data points.", SEDFileName);
        log_close();
        exit(1);
    }

    /* Go back to the beginning */
    rewind(fp);

    /* skip header again */
    if( !fgets(buf, sizeof buf, fp) ){
        printf(" Error. Could not read header from file '%s'. Exiting.\n", SEDFileName);
        log_error("Could not read header from file '%s'. Exiting.", SEDFileName);
        log_close();
        exit(1);
    }

    /* allocate the global SED arrays:
     *   photonEnergy [eV]      - the photon-energy grid (SED x-axis)
     *   sedLuminosity [eV/s/eV] - the spectral luminosity  (SED y-axis) */
    photonEnergy = (double*) malloc(SEDLineCount * sizeof(double));
    sedLuminosity = (double*) malloc(SEDLineCount * sizeof(double));

    /* read data */
    for(i=0; i<SEDLineCount; i++){
        if( fscanf(fp,"%le %le", &photonEnergy[i], &sedLuminosity[i]) != 2 ){
            printf(" Error. Could not read from file '%s'. Exiting.\n", SEDFileName);
            log_error("Could not read from file '%s'. Exiting.", SEDFileName);
            log_close();
            exit(1);
        }
    }

    if(DEBUG){
        log_debug("Printing SED:");
        log_debug("# Photon E [eV]\tSpectral luminosity [eV/s/eV]");
        for(i=0;i<SEDLineCount;i++){ log_debug(" %e\t%e",photonEnergy[i], sedLuminosity[i]); }
        log_debug("End of SED file");
    }


    /* check energy range of the SED. part 1.
     * If the user-provided range [ELow, EHigh] is out of bounds, throw an error. */
    if(myConfig.sourceELow< photonEnergy[0] || myConfig.sourceEHigh> photonEnergy[SEDLineCount-1]){
        printf("ERROR: Provided energy range [sourceELow, sourceEHigh] is out of bounds. Exiting.\n");
        log_error("Provided energy range [sourceELow, sourceEHigh] is out of bounds.");
        log_close();
        exit(1);
    }

    /* Check energy range of the SED. part 2.
     * If the provided E range is smaller than the range in the SED file, the
     * unnecessary entries are removed.
     *
     * We use an el-cheapo method to resize the dynamic SED arrays.
     * Not elegant but it seems to do the job:
     *
     * 1. find indexes of energies that come closest to the provided limits
     */

    double  absDeltaELow, absDeltaEHigh;
    double  tmpDeltaELow    = 0.0;
    double  tmpDeltaEHigh   = 0.0;
    int     indexELow       = 0;
    int     indexEHigh      = 0;

    absDeltaELow  = fabs(myConfig.sourceELow  - photonEnergy[0]);
    absDeltaEHigh = fabs(myConfig.sourceEHigh - photonEnergy[0]);

    for(i=0; i<SEDLineCount; i++){

        tmpDeltaELow  = fabs(myConfig.sourceELow - photonEnergy[i]);
        tmpDeltaEHigh = fabs(myConfig.sourceEHigh - photonEnergy[i]);

        if (tmpDeltaELow<=absDeltaELow && photonEnergy[i]<HIONIZEeV){
            absDeltaELow = tmpDeltaELow;
            indexELow = i;
        }

        if(tmpDeltaEHigh<=absDeltaEHigh){
            absDeltaEHigh = tmpDeltaEHigh;
            indexEHigh = i;
        }

    }
    printf(" Energy range closest to user input is [%e, %e]\n", photonEnergy[indexELow], photonEnergy[indexEHigh]);

    /* 2. count elements  */
    int newSEDLineCount = indexEHigh - indexELow+1;
    if(newSEDLineCount != SEDLineCount){printf(" Resizing SED arrays accordingly\n");}

    /* 3. allocate memory for tmp E and lambda arrays here */
    double *tmpPhotonEnergy = (double*) malloc(newSEDLineCount * sizeof(double));
    double *tmpSedLuminosity = (double*) malloc(newSEDLineCount * sizeof(double));

    /* 4. copy content */
    for(i=0; i<newSEDLineCount; i++){
        tmpPhotonEnergy[i] = photonEnergy[i+indexELow];
        tmpSedLuminosity[i] = sedLuminosity[i+indexELow];
    }

    /* 5. resize original arrays and copy content back */
    photonEnergy = (double*)realloc(photonEnergy, newSEDLineCount * sizeof(double));
    sedLuminosity = (double*)realloc(sedLuminosity, newSEDLineCount * sizeof(double));

    for(i=0; i<newSEDLineCount; i++){
        photonEnergy[i] = tmpPhotonEnergy[i];
        sedLuminosity[i] = tmpSedLuminosity[i];
    }

    SEDLineCount = newSEDLineCount;

    /* build the cached interpolation splines over the (now final) photonEnergy/sedLuminosity
     * arrays. Must happen before sed_compute_norm(), which calls
     * sed_get_log_value() through the integrator. */
    sed_init_spline();

    /* compute normalization  */
    SEDNorm = sed_compute_norm();
    printf(" Integral over SED (normalization) : %le eV/sec\n", SEDNorm);

    /* clean up */
    free(tmpPhotonEnergy);
    free(tmpSedLuminosity);
    fclose(fp);
}


/***************************************************************
 * compute SED normalization
 ***************************************************************/
double sed_compute_norm(){

    /***************************************************************
     * This function computes the normalization of the SED
     * in the (approximated) user-provided limits.
     * It uses the GSL QAGS integrator. If the desired relative
     * tolerance cannot be achieved, it will be increased and then
     * integration will be re-run.
     ***************************************************************/

    /* logging */
    if(DEBUG){
        log_debug("Entering sed_compute_norm() to compute SEDNorm");
        log_debug("Integrating provided SED in the limits of (%le,%le) [eV] ", myConfig.sourceELow, myConfig.sourceEHigh);
    }

    if (myConfig.settingsStroemgrenTest){
        /* here we by-pass the integration for the Strömgren sphere test (delta-function SED). */
        double maxE = 0.0;
        int i;
        for(i=0;i<SEDLineCount;i++){
            if (sedLuminosity[i]>maxE)
                maxE = sedLuminosity[i];
        }

        return maxE;
    }

    /* save error handler to old_handler, then switch it off */
    gsl_error_handler_t * oldErrorHandler = gsl_set_error_handler_off();

    /* set up integration workspace, variables and function */
    gsl_integration_workspace *wnorm = gsl_integration_workspace_alloc(1000);

    double result, error;
    double alpha      =  0.;                            // alpha is only here for legacy reasons, TODO: get rid of it!
    double absError   = 1e-4;
    double relError   = 1e-6;                           // relative error
    double lowerLimit = log(photonEnergy[0]);
    double upperLimit = log(photonEnergy[SEDLineCount-1]);

    gsl_function  myFunction;
    myFunction.function     = &sed_get_log_value;
    myFunction.params       = &alpha;

    int status = 1;

    while(status){

        /* if successful, 0 is returned. If the relError is too small, status will be 13 */
        status     = gsl_integration_qags(&myFunction, lowerLimit , upperLimit, absError, relError, 1000, wnorm, &result, &error);
        relError  *= 5;

        if(status){
            if(DEBUG){
                    log_debug("gsl_integration_qags status = %d", status);
                    log_debug("SEDNorm integration: Increased relError to %e",relError);
            }
        }
    }

    gsl_set_error_handler(oldErrorHandler); // reset error handler (might be unnecessary).
    gsl_integration_workspace_free (wnorm);

    if(DEBUG){
        log_debug("Integrated SED = %e", result);
    }

    return result;

    /* references for:
     * the integrator: https://www.gnu.org/software/gsl/manual/html_node/QAGS-adaptive-integration-with-singularities.html
     * error handling: https://lists.gnu.org/archive/html/help-gsl/2006-04/msg00055.html
     */
}


/***************************************************************
 * compute an estimated radius of the Stroemgren's sphere
 ***************************************************************/
double sed_estimate_stroemgren_radius( double redshift, double time ){

        // find density
      double density = n_H0  * OverDensity * pow3(1 + redshift);

      // find number of emitted ionizing photons
      double ionPhotonN = (SEDNorm/13.6) * (time)*MYR;

      // find ionized volume
      double volume = ionPhotonN / density;

      // volume [cm^3] --> radius [kpc]
      double radius = pow(0.75 * volume / M_PI, 1./3.) / KPC;

      return radius;

}
