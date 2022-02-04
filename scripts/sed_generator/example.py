#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sed import sed



if __name__ == "__main__":

    case = 1


    # -----------------------------------------------------------------
    # SED comparison for BEARS pipeline paper
    # -----------------------------------------------------------------
    # if case == 0:
    #
    #
    #     dir_BB_8  = 'BB_8'
    #
    #     #dir_IMF_10 = 'IMF_redux'
    #     dir_IMF_10 = 'IMF_redux_fEsc0.10_tmp'
    #     dir_IMF_05 = 'IMF_redux_fEsc0.05'
    #
    #     dir_PL_1_5 = 'PL_alpha1.5'
    #     dir_PL_1_0 = 'PL_alpha1.0'
    #
    #     dir_CE_1_0 = 'PL_alpha1.0_below_z8.5+IMF_redux'
    #     dir_CE_1_5 = 'PL_alpha1.5+IMF_redux_fEsc0.05'
    #
    #     M           = 10
    #     z           = 8.0
    #
    #     sedFile    ='sed_%s_M%.3f_z%.3f.dat'%('BB',M,z)
    #
    #     fileList =  [dir_BB_8   +'/'+ 'sed_%s_M%.3f_z%.3f.dat'%('BB',M,z) ,
    #     dir_IMF_10 +'/'+ 'sed_%s_M%.3f_z%.3f.dat'%('IMF',M,z),
    #     dir_IMF_05 +'/'+ 'sed_%s_M%.3f_z%.3f.dat'%('IMF',M,z),
    #     dir_PL_1_5 +'/'+ 'sed_%s_M%.3f_z%.3f.dat'%('PL',M,z),
    #     dir_PL_1_0 +'/'+ 'sed_%s_M%.3f_z%.3f.dat'%('PL',M,z),
    #     dir_CE_1_5 +'/'+ 'sed_%s_M%.3f_z%.3f.dat'%('IMF+PL',M,z),
    #     dir_CE_1_0 +'/'+ 'sed_%s_M%.3f_z%.3f.dat'%('IMF+PL',M,z),
    #     ]
    #
    #     legend=['BB 8', 'IMF 10',  'IMF 05',  'PL 1.5', 'PL 1.0', 'CE 1.5', 'CE 1.0']
    #
    #     colors = ['black','blue','blue','green', 'green', 'red', 'red', 'gray']
    #     ls = ['-', '-', '--', '-', '--', '-', '--', '-', '--']
    #
    #     outFile =  'SED_comparison_M%.3f_z%.3f.pdf'%(M, z)
    #
    #     sed.plot_SED(fileName=fileList,
    #                  eHigh=1e4, eLow=10.4,
    #                  logX=True, logY=True,
    #                  xLimit=[12,1e4], yLimit=[2e48,8e53],
    #                  altFileName=outFile,
    #                  legend=legend,
    #                  legendLoc=1,
    #                  colors=colors,
    #                  lineStyles=ls)

    # -----------------------------------------------------------------
    # BB comparison
    # -----------------------------------------------------------------
    if case == 1:

        # variable star mass, z = 9, log M_halo = 11
        z = 9.
        logM = 11.0

        starM = 10.0
        f1 = 'sed_BB_Mhalo%.3f_z%.3f_Mstar%.3f.dat' %(logM, z, starM )
        sed.generate_SED_stars_BB(haloMass=10**logM, redshift=z, starMass=starM, eHigh=1.e4, eLow=10.4, fileName=f1, N=1000, logGrid=True)

        starM = 200.0
        f2 = 'sed_BB_Mhalo%.3f_z%.3f_Mstar%.3f.dat' %(logM, z, starM )
        sed.generate_SED_stars_BB(haloMass=10**logM, redshift=z, starMass=starM, eHigh=1.e4, eLow=10.4, fileName=f2, N=1000, logGrid=True)

        starM = 700
        f3 = 'sed_BB_Mhalo%.3f_z%.3f_Mstar%.3f.dat' %(logM, z, starM )
        sed.generate_SED_stars_BB(haloMass=10**logM, redshift=z, starMass=starM, eHigh=1.e4, eLow=10.4, fileName=f3, N=1000, logGrid=True)

        sed.plot_SED(fileName=[f1, f2, f3],
                     logX=True, logY=True,
                     yLimit=[2e51, 1e57], xLimit=[1e1, 1e4],
                     legend=['BB 10', 'BB 200', 'BB 700'],
                     legendLoc=1,
                     altFileName="SED_comparison_BB_1.pdf")

        # variable host halo mass
        z = 7.
        logM = 8.0

        f1 = 'sed_BB_M%.3f_z%.3f.dat' %(logM, z )
        sed.generate_SED_stars_BB(haloMass=10**logM, redshift=z, starMass=200, eHigh=1.e4, eLow=10.4, fileName=f1, N=1000, logGrid=True)

        logM = 10.0
        f2 = 'sed_BB_M%.3f_z%.3f.dat' %(logM, z )
        sed.generate_SED_stars_BB(haloMass=10**logM, redshift=z, starMass=200, eHigh=1.e4, eLow=10.4, fileName=f2, N=1000, logGrid=True)

        logM = 12.0
        f3 = 'sed_BB_M%.3f_z%.3f.dat' %(logM, z )
        sed.generate_SED_stars_BB(haloMass=10**logM, redshift=z, starMass=200, eHigh=1.e4, eLow=10.4, fileName=f3, N=1000, logGrid=True)

        logM = 14.0
        f4 = 'sed_BB_M%.3f_z%.3f.dat' %(logM, z )
        sed.generate_SED_stars_BB(haloMass=10**logM, redshift=z, starMass=200, eHigh=1.e4, eLow=10.4, fileName=f4, N=1000, logGrid=True)


        sed.plot_SED(fileName=[f1 ,f2 ,f3 ,f4], logX=True, logY=True, yLimit=[2e47, 5e60], xLimit=[1e1 ,1e4], legend=['BB M08' ,'BB M10' ,'BB M12' ,'BB M14'] ,legendLoc=1, altFileName="SED_comparison_BB_2.pdf")

    # -----------------------------------------------------------------
    # Comparison of PL, BB, IMF, IMF+PL, const z, const halo mass
    # -----------------------------------------------------------------
    if case == 2:

        z = 7.
        logM = 10.0

        f1 = 'sed_PL_M%.3f_z%.3f.dat' %(logM, z )
        f2 = 'sed_IMF_M%.3f_z%.3f.dat' %(logM, z )
        f3 = 'sed_BB_M%.3f_z%.3f.dat' %(logM, z )
        f4 = 'sed_IMF+PL_M%.3f_z%.3f.dat' %(logM, z )


        sed.generate_SED_PL(haloMass=10**logM, eHigh=1.e4, eLow=10.4, fileName=f1, alpha=1.5, N=1000, logGrid=False, qsoEfficiency=0.1)


        sed.generate_SED_stars_IMF(haloMass=10**logM, redshift=z, eLow=10.4, eHigh=1.e4, N=1000,  logGrid=True,
                               starMassMin=5, starMassMax=100, imfBins=99, imfIndex=2.35, fileName=f2)

        sed.generate_SED_stars_BB(haloMass=10**logM, redshift=z, starMass=150, eHigh=1.e4, eLow=10.4, fileName=f3, N=1000, logGrid=True)


        sed.generate_SED_IMF_PL(haloMass=10**logM, redshift=z, eLow=10.4, eHigh=1.e4, N=1000,  logGrid=True,
                            starMassMin=30, starMassMax=100, imfBins=99, imfIndex=2.35,
                            alpha=1.0, qsoEfficiency=0.1,
                            fileName=f4)

        sed.plot_SED(fileName=[f1 ,f2 ,f3 ,f4], logX=True, logY=True, yLimit=[2e47, 5e54], xLimit=[1e1 ,1e4], legend=['PL' ,'IMF' ,'BB' ,'PL + IMF'] ,legendLoc=1, altFileName="SED_comparison_4_types.pdf")

    # -----------------------------------------------------------------
    # Test of the IMF SED
    # -----------------------------------------------------------------
    if case == 3:
        z = 7.
        logM = 8.0

        a = 'sed_IMF_M%.3f_z%.3f'%(logM, z)

        sed.generate_SED_stars_IMF(haloMass=10**logM, redshift=z, eLow=10.4, eHigh=1.e4, N=1000,  logGrid=False,
                                   starMassMin=5, starMassMax=100, imfBins=99, imfIndex=2.35, fileName=a)

        sed.plot_SED(fileName=a,
                     logX=True, logY=True,
                     yLimit=[1e25, 1e55],
                     legend='IMF',
                     legendLoc=1,
                     altFileName="SED_IMF_test.pdf")


