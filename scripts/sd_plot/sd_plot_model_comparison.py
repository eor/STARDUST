#!/usr/bin/python
# -*- coding: utf-8 -*-


from sd_plot import multi_plot

fileDir    = '../../tests/'
outFile    = './test_model_comparison_M11_z7_log.pdf'        # or png?                 
inputFiles = [ 'test_1_PL_source/test_1_profile_M11.000_z7.000_t5.000.dat', 
               'test_1_PL_source/test_1_profile_M11.000_z7.000_t10.000.dat',
               'test_2_IMF_source/test_2_profile_M11.000_z7.000_t5.000.dat',
               'test_2_IMF_source/test_2_profile_M11.000_z7.000_t10.000.dat',
               'test_3_IMF+PL_source/test_3_profile_M11.000_z7.000_t5.000.dat',
               'test_3_IMF+PL_source/test_3_profile_M11.000_z7.000_t10.000.dat'
              ] 


labels     = ['PL\ source, t = 5 Myr',
              'PL\ source, t = 10 Myr',
              'IMF\ source, t = 5 Myr',
              'IMF\ source, t = 10 Myr',
              'IMF+PL\ source, t = 5 Myr',
              'IMF+PL\ source, t = 10 Myr',
             ]


colors = ['black','black','red','red','blue','blue']
ls     = ['-', '-.', '-', '-.', '-', '-.']


rLimits    = [1.0, 750.0]       # lower and upper bound on radius

ylimits    = [[1e-6,1.1], [1e-6,1.1], [1e-6,1.1], [1e-6,1.1], [1e-6,1.1], [-0.2,6], [1.2,2.0], [None,None] ]

for i in range(0,len(inputFiles)):
  inputFiles[i] = fileDir + inputFiles[i]

print inputFiles

multi_plot(inputFiles, outFile, rLimits, ylimits, logT=True, logFractions=True, legends=labels,showLegend=True, ncolsLegend=3, colors=colors, ls=ls)

print "all done"



