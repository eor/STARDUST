/***************************************************************
 * libraries
 ***************************************************************/
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <stdlib.h>

/***************************************************************
 * GSL header files
 ***************************************************************/
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>

/***************************************************************
 * SD headers
 ***************************************************************/
#include "table_settings.h"
#include "constants.h"
//#include "ODE_solver.h"
//#include "Patterson.h"
#include "prototype.h"
#include "allvars.h"
#include "log.h"

/***************************************************************
 * Etc
 ***************************************************************/
#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))

#define find_max(a,b) ((a) > (b) ? (a) : (b))   /* from GSL */
#define find_min(a,b) ((a) < (b) ? (a) : (b))



/* The dimension of the set of differential equation  */
#define NEQS  4  


/***************************************************************
 * Boost headers etc
 ***************************************************************/
#include <iostream>
#include <fstream>
#include <utility>

#include <boost/numeric/odeint.hpp>

#include <boost/phoenix/core.hpp>

#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>

using namespace std;
using namespace boost::numeric::odeint;
namespace phoenix = boost::phoenix;


typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;


#include "ode_derivatives.h" 
#include "ode_jacobian.h"


void rt_main_run(){
    
    int i;
    
    /***************************************************************
     * set up redshift list
     ***************************************************************/
    float zMin  = (float) myConfig.redshiftLow;     
    float zMax  = (float) myConfig.redshiftHigh;    
    float dz    = (float) myConfig.redshiftStride; // TODO: catch dz=0 and dz>(zMax-zMin) if zMin!=zMax
    int   zNum  = roundf((zMax - zMin)/dz) + 1;
    float *zList= (float*) malloc( sizeof(float) * zNum );
    
    if(DEBUG)log_debug("zList contains %d element(s):", zNum); 
    
    for(i=0; i<zNum; i++){
          zList[i] = zMin + i* dz;
          if(DEBUG)log_debug(" zList[%d]\t= %f", i, zList[i]);  
    }      
             
    /***************************************************************
     * GSL definitions  // TODO: refactor these
     ***************************************************************/
    int statusFindRedshift;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T_Findz;
    gsl_root_fsolver *s_Findz;
    double z_lo, z_hi;
    gsl_function F_Findz;
    struct Findz_params findzparams;

    F_Findz.function= &Findz;
    F_Findz.params  = &findzparams;     
    T_Findz         = gsl_root_fsolver_brent;             

    
    /***************************************************************
     * local variables
     ***************************************************************/
    
    /* physics */
//     double l1 = log(HIONIZEeV);
//     double l2 = log(He1IONIZEeV);
//     double l3 = log(He2IONIZEeV);
//     double l4 = log(myConfig.sourceEHigh);

    double hnu21_k = 0.068;         // this should go into constants.h
    double C_H, C_e, C_p;
    double kappa, gamma_e, gamma_p;   
    double ypara;    
        
    
    /* time & redshift */
    int     zIndex;
    double  srcAge;                                 
    double  zSrcTurnOn;   
    
    double  exeTimeCount = 0.0;                     // measure exec time
    double  tmpTimeCount = 0.0;
    clock_t startTime;  
    

    /* file writing */
    FILE    *fpOutFile      = NULL;  
    int     writeCounter    = 1;    
    double  writeInterval   = myConfig.settingsWriteT;         // in [Myr]
    int     isLastTimestep  = 0;
    
    
    /* parameters */ 
    double  timeStep    = myConfig.settingsDeltaT;          // ODE time step in [Myr]
    double  radialStep  = myConfig.settingsDeltaR;          // Spatial resolution of the simulation in [kpc] 
    double  maxRadius   = myConfig.settingsRMax;            // Maximum radius radially away from the radiating source in [kpc]
    double  startRadius = myConfig.settingsRStart;          // in [kpc]. Everything before this radius is assumed to be ionized.
    double  T_CMB0      = myConfig.cosmoTCMB;
    

     /* ODE-related */ 
    double absErrorODE  = 1e-5;
    double relErrorODE  = 1e-5;
    double tStartODE    = 0.0;
    double tEndODE      = timeStep;
    double dtODE        = 0.1 * tEndODE;
   

        
    /***************************************************************
     * Main redshift loop starts here
     ***************************************************************/
    
    for (zIndex = 0; zIndex < zNum; zIndex++){                       
 
        z = zPrevious = zSrcTurnOn = zList[zIndex];        // initialize current, previous (both used when moving a long the grid), source turn on redshifts
                                                           //  at first, z= zPrevious, zPrevious is changed below.  While moving away from the source, we correct for redshift
                                                            
        printf(" Starting run for z = %.3f - %d of %d redshift(s)\n", z, zIndex+1, zNum);
        log_main("Starting run for z = %.3f - %d of %d redshift(s)", z, zIndex+1, zNum); 
        

        /* initializing computing grids */
        rt_initialize_grids(zSrcTurnOn);

        
        srcAge = 0.;        
        
        /* A constant used in the temperature evolution part of the ODE. Used in the Compton cooling term. eq B29 in [2],  */
        lamcons = 4 * k_BeV * (M_PI * M_PI / 15) * pow3( (k_B * T_CMB0 * (1 + z) * 2 * M_PI / (PLANCKCONSTCGS * LIGHTVEL)) ) 
                    * (k_B * T_CMB0 * (1 + z) / (MASSELECTRON * LIGHTVEL * LIGHTVEL)) * LIGHTVEL * THOMSON_ELEC_CROSS;
        
                    //what about n_e and (T - T_CMB*(1+z)) ???? --> those are multiplied in the solver equation      // TODO: maybe clean this up?
                    

        /* 
         * Due to the fact that light can only propagate as far as d=(\Delta t)*c,
         * we introduce one additional grid index, iGridTmpLimit, and only solve the RT only for radii
         * corresponding to indices <= iGridTmpLimit. As the sources ages, i.e.
         * srcAge increases, iGridTmpLimit will increase as well.
         */ 

        int iGridTmpLimit = 0;  
    
        int deltaGrid = (int) ( LIGHTVEL *timeStep*MYR / (radialStep * MPC /1e3 ) ); //3.153e13 * 3e10 * 1e-3 / 3.086e20
        if(DEBUG)printf("DEBUG: deltaGrid = %d\n", deltaGrid);
    

        /***************************************************************
         * Time loop starts here
         ***************************************************************/
        while (isLastTimestep == 0 || srcAge >= myConfig.sourceLifetime){

            startTime = clock();

            iGridTmpLimit += deltaGrid;
            if (iGridTmpLimit > numGridPoints) iGridTmpLimit = numGridPoints;    /* sanity check */

            if(DEBUG)printf("DEBUG: z = %.3f \t srcAge = %.3f Myr\t iGridTmpLimit = %d \t numGridPoints = %d \n", z, srcAge, iGridTmpLimit, numGridPoints);


            /***************************************************************
             * Loop over radii starts here
             ***************************************************************/
            for(iGrid=0; iGrid < iGridTmpLimit; iGrid++){

                /* calculate current radius */
                radius = startRadius + (iGrid) * radialStep;                          
            
                /* compute the cumulative sums of the neutral hydrogen & helium fractions (column densities) */
                nHX1  = 1e-50;            
                nHeX2 = 1e-50;

                 for (i = 0; i < iGrid; i++){                
                    nHX1  += n_H1[i]  * radialStep;
                    nHeX2 += n_He1[i] * radialStep ;
                }               

                nHX13 = (nHX1 + 64. * nHeX2);       /* The factor of 64 is basically sigma_HeII/sigma_HI */

                nHX1  = log(nHX1);                  /* values are tabulated as LOG in the tables, therefore  */    
                nHeX2 = log(nHeX2); 
                nHX13 = log(nHX13); 
                

                if( nHX1 > (UPLIM1 - 3 * TABLERES) ) nHX1  = UPLIM1 - 4 * TABLERES;  
                
                if( nHeX2 > (UPLIM1 - 3 * TABLERES) ) nHeX2 = UPLIM1 - 4 * TABLERES;
                
                if( nHX13 > (UPLIM2 - 3 * TABLERES) ) nHX13 = UPLIM2 - 4 * TABLERES;                
                
                if ( nHX1 <= LOWLIM ) nHX1 = LOWLIM + 2 * TABLERES;  
                
                if ( nHeX2 <= LOWLIM ) nHeX2 = LOWLIM + 2 * TABLERES;
                
                if ( nHX13 <= LOWLIM ) nHX13 = LOWLIM + 2 * TABLERES;
                
                
                // Interpolate the table to obtain the params.* values below 
                interpolation(nHX1, nHeX2, nHX13, radius); 

#ifdef STROEMGRENTEST                
                // experimental: override for small radii
                // limit max value of integrals for small radii
                //if (fuku_e1h1[iGrid] > 8e-13) fuku_e1h1[iGrid] = 8e-13;
                
#else                
                if (fuku_e1h1[iGrid] > 1e-11) fuku_e1h1[iGrid] = 1e-11;
                if (fuku_ehe1[iGrid] > 1e-11) fuku_ehe1[iGrid] = 1e-11;
                if (fuku_ehe2[iGrid] > 1e-11) fuku_ehe2[iGrid] = 1e-11;
                
                
                
#endif          
                      
                // load values of fuku integrals into global struct (needed within the ODE solver)    
                params.fe1h1   = fuku_e1h1[iGrid];
                params.fehe1   = fuku_ehe1[iGrid];
                params.fehe2   = fuku_ehe2[iGrid];
                params.intH1   = integral_H1[iGrid];
                params.intHe1  = integral_He1[iGrid];
                params.intHe2  = integral_He2[iGrid];
                params.localOD = over_densities[iGrid];

                
                // setting up initial conditions for the ODE solver 
                vector_type y(4);
                y[0] = x_HII[iGrid];     
                y[1] = x_HeII[iGrid];
                y[2] = x_HeIII[iGrid];
                y[3] = T_e[iGrid];
                

#ifdef ROSENBROCK4
               
                size_t numStepsODE = integrate_const( make_dense_output< rosenbrock4<double> >( absErrorODE , relErrorODE ) ,
                                                       make_pair( ode_system_derivatives() , ode_system_jacobi() ) , 
                                                       y , tStartODE , tEndODE , dtODE 
                                                     );
#endif
                
               
                
#ifdef RKDOPRI5
                
                size_t numStepsODE = integrate_adaptive( make_dense_output < runge_kutta_dopri5 <vector_type>>( absErrorODE , relErrorODE ),
                                                        ode_system_derivatives() , y , tStartODE , tEndODE , dtODE 
                                                     );
#endif       
                
#ifdef RK78
//                 typedef runge_kutta_fehlberg78< double > stepper_RK78
                size_t numStepsODE = integrate_adaptive( make_controlled<runge_kutta_fehlberg78<vector_type> >( absErrorODE , relErrorODE ) ,
                                                        ode_system_derivatives() , y , tStartODE , tEndODE , dtODE 
                                                     );
#endif                
                
                



#ifdef BS
//                 typedef runge_kutta_fehlberg78< double > stepper_RK78
                size_t numStepsODE = integrate_adaptive( bulirsch_stoer<vector_type>( absErrorODE , relErrorODE ) ,
                                                        ode_system_derivatives() , y , tStartODE , tEndODE , dtODE 
                                                     );
#endif                
                
                
                
#ifdef EULER
                
                typedef implicit_euler< double > stepper_IE;

                size_t numStepsODE = integrate_const( stepper_IE() , std::make_pair( ode_system_derivatives() , ode_system_jacobi() ) ,
                                                       y , tStartODE , tEndODE , dtODE
                                                    );
                
#endif                
                

                
                
#ifdef FIREHOSEDEBUG
                
                
                //if(iGrid>102 && iGrid<120){
                cout << "DEBUG: number of steps:  " << numStepsODE 
                << " at iGrid: " << iGrid
                << "\ty0_in = "  << x_HII[iGrid]    << "\ty0_out = " <<y[0] 
                << "\ty1_in = "  << x_HeII[iGrid]   << "\ty1_out = " <<y[1]
                << "\ty2_in = "  << x_HeIII[iGrid]  << "\ty2_out = " <<y[2]
                << "\ty3_in = "  << T_e[iGrid]      << "\ty3_out = " <<y[3]
                << endl;
                //}              
                
#endif                
               
                // store ODE output 
                x_HII[iGrid]    = y[0]; 
                x_HeII[iGrid]   = y[1];  
                x_HeIII[iGrid]  = y[2];  
                T_e[iGrid]      = y[3];
                
                
#ifdef STROEMGRENTEST    
                if (isnan(x_HII[iGrid]) || isinf(x_HII[iGrid]) ){
                
                    // check if gas should be ionized. If yes, set x_HII = 1                    
                    double radiusEst = sed_estimate_stroemgren_radius( z , srcAge + timeStep );
                    
                    if (radius < 0.8 * radiusEst) x_HII[iGrid] = 1.0;
                   
                }     
                
#endif                

                
                // check for numerical instabilities, to regularize the solutions if necessary 
                x_HII[iGrid]    = GSL_MAX_DBL( x_HII[iGrid], 1.e-55 );
                x_HII[iGrid]    = GSL_MIN_DBL( x_HII[iGrid], 1.0 );
                
                x_HeII[iGrid]   = GSL_MAX_DBL( x_HeII[iGrid], 1.e-15 );
                x_HeII[iGrid]   = GSL_MIN_DBL( x_HeII[iGrid], 1.0 );
                
                x_HeIII[iGrid]  = GSL_MAX_DBL( x_HeIII[iGrid], 1.e-15 );                
                x_HeIII[iGrid]  = GSL_MIN_DBL( x_HeIII[iGrid], 1.0 );  
                
                x_HI[iGrid]     = 1.0 - x_HII[iGrid];
                
                x_HI[iGrid]     = GSL_MAX_DBL( x_HI[iGrid], 1.e-55 );
                x_HI[iGrid]     = GSL_MIN_DBL( 1.0, x_HI[iGrid] );

                x_HeI[iGrid]    = 1.0 - x_HeII[iGrid] - x_HeIII[iGrid];
                x_HeI[iGrid]    = GSL_MAX_DBL( x_HeI[iGrid], 1.e-15 );

                // update n_H1 and n_He1:                
                n_H1[iGrid]  = n_H0  * over_densities[iGrid] * pow3(1 + z) * x_HI[iGrid] ;
                n_He1[iGrid] = n_He0 * over_densities[iGrid] * pow3(1 + z) * x_HeI[iGrid];
                

#ifdef STROEMGRENTEST

                x_HeI[iGrid]   = 0.0;
                x_HeII[iGrid]  = 0.0;
                x_HeIII[iGrid] = 0.0;        
                T_e[iGrid]     = STROEMGRENTEMP; 
#endif  
       
            }  
            /***************************************************************
            * Loop over radii ends here
            ***************************************************************/  

            
            
            
            
            
            /***************************************************************
            * File writing happens here
            ***************************************************************/            
            // note that 
            // 1) after having passed the ODE solver,  the current source age is now srcAge + timeStep
            // 2) without the /2.0 in the if statement files are written one timestep prior to the desired time.            
            
            isLastTimestep = rt_check_is_last_time_step(srcAge, timeStep);
            
            if( ( writeCounter*writeInterval - (srcAge + timeStep) ) < timeStep/2.0 || isLastTimestep == 1 ){               
                
        
                writeCounter++;
                
                // write the fuku integrals to a log file
                if(DEBUG){           
                    
                    char *tmpFukuLog = (char*) malloc( sizeof(char) * 256);
                    sprintf(tmpFukuLog, "log_fuku_integrals_M%.3f_z%.3f_t%.3f.dat", myConfig.haloMass, zSrcTurnOn, (srcAge+ timeStep));    
                    tmpFukuLog = utils_concat_path(myConfig.pathOutDir, tmpFukuLog);
                
                    FILE *fpFukuLog = fopen(tmpFukuLog, "w");
                    fprintf(fpFukuLog, "# R[kpc]  fe1h1 \t  fehe1 \t   fehe2 \t  intH1 \t  intHe1 \t intHe2\n");
                
                    for(i=0;i<iGridTmpLimit;i++){
                        radius = startRadius + (i) * radialStep;    
                        fprintf(fpFukuLog, "%.2f\t%e\t%e\t%e\t%e\t%e\t%e\n", radius, fuku_e1h1[i], fuku_ehe1[i], fuku_ehe2[i], integral_H1[i], integral_He1[i], integral_He2[i]  );  
                    
                    }  
                    fclose(fpFukuLog);                    
                    log_debug("Wrote '%s' at t = %.3f Myr", tmpFukuLog, (srcAge + timeStep) ) ;
                    free(tmpFukuLog);
 
                }
                
                                
                // generate output file name and open file 
                sprintf(outFile, "profile_M%.3f_z%.3f_t%.3f.dat", myConfig.haloMass, zSrcTurnOn, (srcAge+ timeStep));           
                
                if (myConfig.pathID[0] == '\0')
                    outFile = utils_concat_path_noID(myConfig.pathOutDir, outFile);
                else
                    outFile = utils_concat_path(myConfig.pathOutDir, outFile);
                
                
                
                fpOutFile  = fopen(outFile, "w");
                                
                printf(" Writing '%s' at t = %.3f Myr\n", outFile, (srcAge + timeStep) ) ;
                
                
                double tmpR = 0.0;

#ifdef SHORTPROFILES
                for(i=0;i<iGridTmpLimit;i++){
#else
                for(i=0;i<numGridPoints;i++){
#endif

                    /* Conversion to fractions */
                    //x_HII[i]        = n_H2[i]  / ( n_H0  * over_densities[i] * pow3(1 + z) );
                    #ifdef STROEMGRENTEST
                    x_HeI[i]       = 0.0;           // necessary because n_He0 = 0.
                    x_HeII[i]      = 0.0;
                    x_HeIII[i]     = 0.0;
                    #else
                    /*x_HeI[i]       = n_He1[i] / ( n_He0 * over_densities[i] * pow3(1 + z) );
                    x_HeII[i]       = n_He2[i] / ( n_He0 * over_densities[i] * pow3(1 + z) );
                    x_HeIII[i]       = n_He3[i] / ( n_He0 * over_densities[i] * pow3(1 + z) );
                    */
                    
                    ne[i] = (1.0 - x_HI[i]) * n_H0 * over_densities[i] * pow3(1+z)
                            + (x_HeII[i] + 2 * x_HeIII[i]) * n_He0*over_densities[i] * pow3(1+z);
                    #endif
                    
                    /* compute spin and brightness temperatures */
                    kappa   =  3.1e-11 * pow (T_e[i], .357) * exp (-32. / T_e[i]);
                
                    C_H     = x_HI[i]*n_H0*over_densities[i]*pow3(1+z)*kappa;//(n_H0 * over_densities[i] * pow3((1 + z)) - n_H2[i]) * kappa;

                    gamma_e = pow (10., (-9.607 + .5 * log10 (T_e[i]) * exp(-log10(pow (T_e[i], 4.5) / 1800))));
                    
                    C_e     = ne[i] * gamma_e;

                    gamma_p = 3.2 * kappa;
                    
                    C_p     = x_HII[i]*(n_H0*over_densities[i]*pow3(1+z)) * gamma_p;

                    ypara   = hnu21_k * (C_H + C_e + C_p) / (EINSTEIN_A * T_e[i]);
                    
                    T_spin[i] = ( hnu21_k + T_CMB0 * (1 + z) + ypara * T_e[i] ) / (1 + ypara);

                    T_brig[i] = 16 * x_HI[i] * ( 1 - T_CMB0 * (1 + z) / T_spin[i] ) * pow((1 + z) / 10., .5);
                        
                    
                    /* write line by line */
                    tmpR = startRadius + (i) * radialStep;
                    fprintf (fpOutFile, "%le %f %le %le %le %le %le %le %le %le \n", 
                            z, tmpR, x_HI[i], 
                            x_HII[i], x_HeI[i],  x_HeII[i], x_HeIII[i], 
                            T_e[i], T_spin[i], T_brig[i]);          
                    
                }
                
                fclose(fpOutFile);
                
                // some output
                tmpTimeCount = exeTimeCount + ((double) clock() - startTime) / (60.0 * CLOCKS_PER_SEC);                
                log_main("Wrote '%s' at t = %.3f Myr after %.2f min", outFile, (srcAge + timeStep), tmpTimeCount ) ;
                printf(" Time elapsed after %.3f Myr: %.2f min\n", (srcAge + timeStep), tmpTimeCount);
            } 
            
            /* update source age and store previous redshift value */
            srcAge   += timeStep;
            zPrevious = z;           
            
            
            
            /***************************************************************
             * find the new redshift TODO: put in extra function!
             ***************************************************************/              
            
            s_Findz = gsl_root_fsolver_alloc (T_Findz);

            z_lo = 4.0, z_hi = zSrcTurnOn;          /* limits within which to look for z  */

            findzparams.cumulative_source_lifetime = srcAge;
            findzparams.switchon_source_redshift = zSrcTurnOn;

            gsl_root_fsolver_set (s_Findz, &F_Findz, z_lo, z_hi);
            
            do{
                iter++;
                statusFindRedshift = gsl_root_fsolver_iterate (s_Findz);
                z    = gsl_root_fsolver_root (s_Findz);
                z_lo = gsl_root_fsolver_x_lower (s_Findz);
                z_hi = gsl_root_fsolver_x_upper (s_Findz);
                statusFindRedshift = gsl_root_test_interval (z_lo, z_hi, 0, 0.001);
                
                    
            }
            while (statusFindRedshift == GSL_CONTINUE && iter < max_iter);
        
            gsl_root_fsolver_free (s_Findz);
            

            
            /* for the Strömgen sphere test we want a static Universe */
            #ifdef STROEMGRENTEST
                z = zPrevious = zSrcTurnOn = zList[zIndex];
            #endif    
         
            /* determine time spent for this step and add it to the total */
            exeTimeCount += ((double) clock() - startTime) / (60.0 * CLOCKS_PER_SEC);
            
        }       
        /***************************************************************
        * Loop over source age ends here
        ***************************************************************/  
        

        printf (" Total run time for z=%.3f : %.2f min \n", zSrcTurnOn, exeTimeCount);
        log_main("Total run time for z=%.3f : %.2f min ", zSrcTurnOn, exeTimeCount ) ;

    }  
    /***************************************************************
     * redshift loop ends here
     ***************************************************************/  

 
}


/***************************************************************
 * Initialize arrays
 ***************************************************************/  
void rt_initialize_grids(double zSrcTurnOn){
    
    int i;
    double T_CMB0 = myConfig.cosmoTCMB;
    
    for (i=0; i<numGridPoints; i++){
        
        x_HeII[i]       = 0.0;
        x_HeIII[i]      = 0.0; 
        ne[i]           = 0.0;
        n_H1[i]         = n_H0  * pow3(1. + zSrcTurnOn) * over_densities[i]; // TODO NFW or other density profiles here
#ifdef STROEMGRENTEST
        
        x_HI[i]         = 1.0 - 1.2e-3; 
        x_HII[i]        = 1.2e-3;           // Illiev 1 test
        x_HeI[i]        = 0.0;
        T_e[i]          = STROEMGRENTEMP;
        n_He1[i]        = 0.0;
#else
        x_HI[i]         = 1.0;            
        x_HII[i]        = 0.0;        
        n_He1[i]        = n_He0 * pow3(1. + zSrcTurnOn) * over_densities[i];
        x_HeI[i]        = 1.0; 
        T_e[i]          = T_CMB0 * pow2(1. + zSrcTurnOn) / (1. + zTkinEQTCMB);
#endif            

        T_spin[i]       = 0.0;
        T_brig[i]       = 0.0;
        
        fuku_e1h1[i]    = 0.0;
        fuku_ehe1[i]    = 0.0;
        fuku_ehe2[i]    = 0.0;
        
        integral_H1[i]  = 0.0;
        integral_He1[i] = 0.0;
        integral_He2[i] = 0.0;            

        }
    
}

/***************************************************************
 * Check if we are in the final time step
 ***************************************************************/
int rt_check_is_last_time_step(double srcAge, double timeStep){
    
    // after ODE solver the time is srcAge + timeStep
    // We are in the last time step, if srcAge+ 2* timeStep > myConfig.sourceLifetime
    
    if (srcAge+ 2*timeStep >= myConfig.sourceLifetime)
        return 1;    
    
    return 0;
}

