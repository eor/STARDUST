#!/usr/bin/python
# -*- coding: utf-8 -*-


from sd_plot import multi_plot

fileDir    = '../tests/'
outFile    = './test_1_PL_evolution_M11_z7_log.pdf'        # or png?                 
inputFiles = [ 'test_1_PL_source/test_1_profile_M11.000_z7.000_t1.000.dat', 
               'test_1_PL_source/test_1_profile_M11.000_z7.000_t3.000.dat', 
               'test_1_PL_source/test_1_profile_M11.000_z7.000_t5.000.dat', 
               'test_1_PL_source/test_1_profile_M11.000_z7.000_t7.000.dat', 
               'test_1_PL_source/test_1_profile_M11.000_z7.000_t9.000.dat' 
              ] 


labels     = []

colors = ['black','black','black','black','black','black','red', 'gray']
ls     = ['-.'    ,  '-.', '-.', '-.', '-.', '-.', '-.', '-.', '-.',]


rLimits    = [1.0, 750.0]       # lower and upper bound on radius   

#ylimits    = [[0,1.05], [0,1.05], [0,1.05], [0,1.05], [0,1.05], [None,None], [None,None], [None,None] ]
ylimits    = [[1e-6,1.1], [1e-6,1.1], [1e-6,1.1], [1e-6,1.1], [1e-6,1.1], [-0.2,6], [1.2,2.0], [None,None] ]


for i in range(0,len(inputFiles)):
  inputFiles[i] = fileDir + inputFiles[i]

print inputFiles

multi_plot(inputFiles, outFile, rLimits, ylimits, logT=True, logFractions=True, legends=labels,showLegend=False, colors=colors, ls=ls)

print "all done"
