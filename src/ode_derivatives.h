struct ode_system_derivatives
{
    
    
    void operator()( const vector_type &y , vector_type &dydt , double /* t */ )
    {

        
        
        double integral_ionisation_rate_H_I     = (params.fe1h1);    //    fuku_e1h1[iGrid]
        double integral_ionisation_rate_He_I    = (params.fehe1);    //    fuku_ehe1[iGrid]
        double integral_ionisation_rate_He_II   = (params.fehe2);    //    fuku_ehe2[iGrid]
        double integral_heating_rate_H_I        = (params.intH1);    //    integral_H1[iGrid]
        double integral_heating_rate_He_I       = (params.intHe1);   //    integral_He1[iGrid]
        double integral_heating_rate_He_II      = (params.intHe2);   //    integral_He2[iGrid]
        double local_over_density               = (params.localOD);  //     over_densities[iGrid]
        
        double y0 = y[0]; 
        double y1 = y[1];
        double y2 = y[2];
        double y3 = y[3];


        y3 = find_max( y3, 0.0 );   // T_e sanity check. should not be negative
       
        double n_H  = n_H0   * local_over_density * pow3(1 + z);
        double n_He = n_He0  * local_over_density * pow3(1 + z);
        double ne  = y0 * n_H + (y1 + 2*y2) * n_He;
        
        double n_Hn, n_Hen;

        double a2h2, R1c;

        double bhe1, bhe2, b1h1;
        double ahe2, ahe3;
        double zhe2;

        double zeta_H1, zeta_He1, zeta_He2;
        double eta_H2, eta_He2, eta_He3;
        double psi_H1, psi_He1, psi_He2;
        double theta_ff, w_He2, mu = 1.24;       // mu = mean molecular weight
        
        double T_CMB0 = myConfig.cosmoTCMB;
        
        // TODO change names, do so also for the jacobian



        /* eq(B6) in Ref [2] */
        a2h2 = 2.59e-13 * pow(y3 / 1.e4, -0.8);                                                    // alpha2HII
        
        /* eq(B1) in Ref [2] */
        b1h1 = 5.85e-11 * pow(y3, .5) * pow(1 + pow(y3 / 1.e5, .5), -1.) * exp(-1.578e5 / y3);    // beta1HI
        
        /* eq(B3) in Ref [2] */
        bhe1 = 2.38e-11 * pow(y3, .5) * pow(1 + pow(y3 / 1e5, .5), -1.) * exp(-2.853e5 / y3);     // betaHeI    
        
        /* eq(B4) in Ref [2] */
        bhe2 = 5.68e-12 * pow(y3, .5) * pow(1 + pow(y3 / 1e5, .5), -1.) * exp(-6.315e5 / y3);     // betaHeII
        
        /* eq(B7) in Ref [2] */
        ahe2 = 1.50e-10 * pow(y3, -.6353);
        
        /* eq(B11) in Ref [2] */
        ahe3 = 3.36e-10 * pow(y3, -.5) * pow(y3 / 1e3, -.2) * pow(1 + pow(y3 / 4e6, .7), -1.);
        
        /* eq(B10) in Ref [2] */
        zhe2 = 1.9e-3 * pow(y3, -1.5) * exp(-4.7e5 / y3) * (1 + 0.3 * exp(-9.4e4 / y3));
        
        /* eq(B28) in Ref [2] */
        theta_ff = 1.42e-27 * 1.1 * pow(y3, .5); // for T between 1e4-1e8K, gaunt factor is between 1.1 and 1.5
        
        /* eq(B17) in Ref [2] */
        zeta_H1 = 1.27e-21 * pow(y3, .5) * pow(1 + pow(y3 / 1.e5, .5), -1.) * exp(-1.58e5 / y3);    

        /* eq(B18) in Ref [2] */
        zeta_He1 = 9.38e-22 * pow(y3, .5) * pow(1 + pow(y3 / 1.e5, .5), -1.) * exp(-2.85e5 / y3);
        
        /* eq(B20) in Ref [2] */
        zeta_He2 = 4.95e-22 * pow(y3, .5) * pow(1 + pow(y3 / 1.e5, .5), -1.) * exp(-6.31e5 / y3);
        
        /* eq(B21) in Ref [2] */
        eta_H2 = 6.5e-27 * pow(y3, .5) * pow(y3 / 1e3, -.2) * pow(1 + pow(y3 / 1.e6, .7), -1.);
        
        /* eq(B22) in Ref [2] */
        eta_He2 = 1.55e-26 * pow(y3, .3647);
        
        /* eq(B24) in Ref [2] */
        eta_He3 = 3.48e-26 * pow(y3, .5) * pow(y3 / 1e3, -.2) * pow(1 + pow(y3 / 4.e6, .7), -1.);
        
        /* eq(B23) in Ref [2] */
        w_He2 = 1.24e-13 * pow(y3, -1.5) * exp(-4.7e5 / y3) * (1 + .3 * exp(-9.4e4 / y3));
        
        /* eq(B25) in Ref [2] */
        psi_H1 = 7.5e-19 * pow(1 + pow(y3 / 1e5, .5), -1.) * exp(-1.18e5 / y3);
        
        /* eq(B26) in Ref [2] */
        psi_He1 = 9.1e-27 * pow(y3, -.1687) * pow(1 + pow(y3 / 1e5, .5),-1.) * exp(-1.31e4 / y3) * ne * y1*n_He;
        
        /* eq(B27) in Ref [2] */
        psi_He2 =  5.54e-17 * pow(y3, -.397) * pow(1 + pow(y3 / 1e5, .5), -1.) * exp(-4.73e5 / y3);
        

 
        /* ---Explicitly regularize ---------- */
        #ifdef STROEMGRENTEST
        n_Hen       = 0.0;
        psi_He1     = 0.0; 
        psi_He2     = 0.0; 
        integral_heating_rate_He_I = 0.0;
        integral_heating_rate_He_II = 0.0;
        integral_ionisation_rate_He_I       = 0.0;
        integral_ionisation_rate_He_II       = 0.0;
        #endif    

        /* --------------------------------- */


        R1c = b1h1 * ne + integral_ionisation_rate_H_I;

        /* eq(26) in Ref [2] */
        dydt[0] = R1c * (1-y0) - a2h2 * ne*ne/n_H;
  
        /* eq(29) in Ref [2] */
        dydt[1] =   (1.0-y1-y2) * integral_ionisation_rate_He_I
                    + bhe1 * ne * (1.0-y1-y2) 
                    - bhe2 * ne * y1 
                    - ahe2 * ne * y1 
                    + ahe3 * ne * y2 
                    - zhe2 * ne * y1;
        
        /* eq(30) in Ref [2] */
        dydt[2] =   y1 * integral_ionisation_rate_He_II
                    + bhe2 * ne * y1 
                    - ahe3 * ne * y2;
        
        
        
        double srH  = 1.0/(1.0 + 4*(0.15/1.9));  // species ratio
        double srHe = 1.0/((1.9/0.15) + 4);      // species ratio

        // the species ratios correspond to the fraction n_H/n_B and n_He/n_B,
        // where n_H and n_He are the number densities of all hydrogen and helium, respectively
        // and n_B is the number density of all baryons. The redshift dependence for these values cancels out.

        
        /* eq(36) in Ref [2] */
        dydt[3] = (  ( (1-y0) * srH        * integral_heating_rate_H_I
                      + (1.0-y1-y2) * srHe * integral_heating_rate_He_I
                      + y1 * srHe          * integral_heating_rate_He_II
                      )
                    - (ne * (zeta_H1 * (1.0-y0) * srH + zeta_He1 * (1.0-y1-y2) * srHe + zeta_He2 * y1 *srHe)) 
                    - (ne * (eta_H2 * y0 * srH + eta_He2 * y1 *srHe + eta_He3 * y2 * srHe)) 
                    - w_He2 * y2 * srHe * ne  
                    - ( ne * (psi_H1 * (1.0-y0)*srH + psi_He1 * (1.0-y1-y2) * srHe + psi_He2 * y1 * srHe) ) 
                    - (y3 - T_CMB0 * (1 + z)) * lamcons * ne 
                    - theta_ff * (y0*srH + y1*srHe + 4 * y2*srHe) * ne
                    - 7.5 * hubble (z) * k_BeV * y3  / mu
                 ) * mu / (k_BeV  * 1.5);
 
                 
        // use [dt] = Myr instead of seconds
        dydt[0] *= MYR;      
        dydt[1] *= MYR;    
        dydt[2] *= MYR;   
        dydt[3] *= MYR;  
        
        
        // testing:
        // if (isnan(dydt[0]) || isinf(dydt[0])) dydt[0]=0.0;

        
        
        
        #ifdef STROEMGRENTEST
        // no Helium (1,2) & constant Temperature (3)
        dydt[1] = dydt[2] = dydt[3] = 0.0;
        #endif             
        
    
        
    }
};
