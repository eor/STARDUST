#!/usr/bin/python
# -*- coding: utf-8 -*-

from sd_plot import multi_plot

profileFile = 'profile_M11.000_z7.000_t10.000.dat'



inputFiles = [ 'BB_8/'+profileFile, 'IMF_10/'+ profileFile,'IMF_05/'+ profileFile , 'PL_1.5/'+ profileFile ,'PL_1.0/'+ profileFile ,'CE_1.5/'+ profileFile , 'CE_1.0/'+ profileFile ] 


#labels     = []
#labels    += ['PL: \mbox{    }   $z=7.000$,\ $\mathrm{M}_{\mathrm{halo}}=10^8 \mathrm{M}_{\odot}$,\ $t=10.0\mathrm{Myr}$']
#labels    += ['IMF:              $z=7.000$,\ $\mathrm{M}_{\mathrm{halo}}=10^8 \mathrm{M}_{\odot}$,\ $t=10.0\mathrm{Myr}$']
#labels    += ['BB:\mbox{ }\, \, \,  \, \, \,$z=7.000$,\ $\mathrm{M}_{\mathrm{halo}}=10^8 \mathrm{M}_{\odot}$,\ $t=10.0\mathrm{Myr}$']
#labels    += ['IMF+PL:           $z=7.000$,\ $\mathrm{M}_{\mathrm{halo}}=10^8 \mathrm{M}_{\odot}$,\ $t=10.0\mathrm{Myr}$']              
              
labels = ['BB 8', 'IMF 10',  'IMF 05',  'PL 1.5', 'PL 1.0', 'CE 1.5', 'CE 1.0']       

outFile    = 'SD_model_comparison_M10_z8_10.0Myr_log.pdf'        # or png?

#colors = ['black','blue','blue','green','green','red','red', 'gray']
#ls     = ['-'    ,  '-', '--', '-', '--', '-', '--', '-', '--',]
# BW friendly
ls     = ['-.'    ,  '-', '--', '-', '--', '-', '--', '-', '--',]
colors = [ 'black','#d73027','#d73027','#fdb863','#fdb863','#42bcf4','#42bcf4']


#inputFiles = ['PL/profile_M8.000_z7.000_t5.000.dat' ,'BB/profile_M8.000_z7.000_t5.000.dat', 'IMF/profile_M8.000_z7.000_t5.000.dat', 'IMF+PL/profile_M8.000_z7.000_t5.000.dat'  ] 
#labels     = ['PL z=7.000, M=8.000, t=5.000Myr', 'BB z=7.000, M=8.000, t=5.000Myr', 'IMF z=7.000, M=8.000, t=5.000Myr', 'IMF+PL z=7.000, M=8.000, t=5.000Myr'  ] 
#outFile    = 'compare_M8_z7_5.0Myr.pdf'        # or png?




rLimits    = [1.0, 450.0]       # lower and upper bound on radius   

#ylimits   = [[0,1], [0,1], [0,1], [0,None], [0,1], [None,None], [None,None], [None,None] ]
#ylimits       = [[0,1.05], [0,1.05], [0,1.05], [0,1.05], [0,1.05], [None,None], [None,None], [None,None] ]
ylimits    = [[1e-6,1.1], [1e-6,1.1], [1e-6,1.1], [1e-6,1.1], [1e-6,1.1], [-0.2,6], [1.2,2.0], [None,None] ]


for i in range(0,len(inputFiles)):
  inputFiles[i] = fileDir + inputFiles[i]

print inputFiles
#multi_plot_SD(inputFiles, outFile, rLimits, ylimits, logT=True, logFractions=True, legends=labels)
multi_plot_SD(inputFiles, outFile, rLimits, ylimits, logT=True, logFractions=True, legends=labels,showLegend=True, colors=colors, ls=ls)
print "all done"
