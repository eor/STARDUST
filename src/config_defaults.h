

#ifndef CONFIG_DEFAULTS_H

#define CONFIG_DEFAULTS_H
    

/* paths */
#define CONFIG_DEFAULT_PATH_OUTDIR      "."
#define CONFIG_DEFAULT_PATH_SED         "SAMPLE_SED.dat"
#define CONFIG_DEFAULT_PATH_DENSITY     "\0"
#define CONFIG_DEFAULT_PATH_ID          ""


/* simulation */
#define CONFIG_DEFAULT_SOURCE_ELOW      13.6        // float, in eV
#define CONFIG_DEFAULT_SOURCE_EHIGH     1.0e4       // float, in EV
#define CONFIG_DEFAULT_SOURCE_LIFETIME  10.0        // float, in Myr

#define CONFIG_DEFAULT_HALOMASS         9.0         // float, log10(M/M_sol)

#define CONFIG_DEFAULT_REDSHIFT_LOW     6.0         // float
#define CONFIG_DEFAULT_REDSHIFT_HIGH    12.0        // float
#define CONFIG_DEFAULT_REDSHIFT_STRIDE  1.0         // float


/* cosmology - defaults are from Planck results XIII, 2015, TT+lowP (Table 4) */
#define CONFIG_DEFAULT_COSMO_OMEGA_M    0.315        //  Total matter density 
#define CONFIG_DEFAULT_COSMO_OMEGA_L    0.685        //  Dark energy / cosmological constant
#define CONFIG_DEFAULT_COSMO_OMEGA_B    0.0491       //  Baryon density
#define CONFIG_DEFAULT_COSMO_H0         67.31        //  Hubble parameter [km/sec/Mpc] at z = 0
#define CONFIG_DEFAULT_COSMO_H100       0.6731       //  Hubble parameter / 100
#define CONFIG_DEFAULT_COSMO_SIGMA8     0.829        //  power spectrum normalization    
#define CONFIG_DEFAULT_COSMO_TAU_THOM   0.078        //  Optical depth of Thomson scattering
#define CONFIG_DEFAULT_COSMO_TCMB0      2.731        //  CMB temperature at z = 0


/* general settings */
#define CONFIG_DEFAULT_SETTINGS_DEBUG   0           // int 0,1
#define CONFIG_DEFAULT_SETTINGS_RMAX    1000.0      // float, kpc
#define CONFIG_DEFAULT_SETTINGS_RSTART  1.0         // float, kpc
#define CONFIG_DEFAULT_SETTINGS_DELTA_R 0.1         // float, kpc
#define CONFIG_DEFAULT_SETTINGS_DELTA_T 0.01        // float, Myr
#define CONFIG_DEFAULT_SETTINGS_WRITE_T 1.0         // float, Myr (write interval)

    
#endif    
