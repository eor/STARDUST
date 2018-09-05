#!/usr/bin/python
# -*- coding: utf-8 -*-

def multi_plot(inputFiles, outFile, xLimits, ylimits, logT=False, logFractions=False, legends=[], showLegend=True, ncolsLegend=4, doColor=True, colors=None, ls=None):
  
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    
    
    # enable these if you wannt nice looking plots
    #from matplotlib import rc
    #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    
    ## for Palatino and other serif fonts use:
    ##rc('font',**{'family':'serif','serif':['Palatino']})
    ##rc('text', usetex=True)
    
    rows = 4
    columns = 2
    stylechange = 5
    plotdistX = 0.1
    plotdistY = 0.1
    labelFontSize = 15
  
    if colors:        
        color = colors
    else:    
        color = ['black','red','blue','green','purple','orange','yellow']        
        
        
    if not ls:
        ls     = ['-'    ,  '-', '--', '-', '--', '-', '--', '-', '--',]
        
        
    f, axarr = plt.subplots(rows, columns, sharex=True, sharey=False, figsize=(10,15))
    xlabel = '$r\, [\mathrm{kpc}]$'


    # labels for fractions
    if logFractions == True:
        ylabelFrac = [ r'$\log\left(x_\mathrm{H_{I}}\right)$',r'$\log\left(x_\mathrm{H_{II}}\right)$',r'$\log\left(x_\mathrm{He_{I}}\right)$',r'$\log\left(x_\mathrm{He_{II}}\right)$',r'$\log\left(x_\mathrm{He_{III}}\right)$']
    else:
        ylabelFrac = [ r'$x_\mathrm{H_{I}}$',r'$x_\mathrm{H_{II}}$',r'$x_\mathrm{He_{I}}$',r'$x_\mathrm{He_{II}}$',r'$x_\mathrm{He_{III}}$',]

    # labels for temperatures
    if logT == True:
        ylabelTemp = [ r'$\log\left(T_{\mathrm{kin}}/\mathrm{K}\right)$', r'$\log\left(T_{\mathrm{Spin}}/\mathrm{K}\right)$', r'$\log\left(\delta T_{\mathrm{B}}/\mathrm{K}\right)$']
    else:
        ylabelTemp = [ r'$T_{\mathrm{kin}}$',r'$T_{\mathrm{Spin}}$',r'$\delta T_{\mathrm{B}}\, [K]$']
  
    #ylabel = ['x_H','x_H_II','x_He','x_He_II','x_He_III','T_e','T_Spin','T_Brig']
    ylabel = ylabelFrac + ylabelTemp




    if logFractions == True:
        for i in range(0,5):
            ylimits[i][0] =  math.log10(ylimits[i][0])
            ylimits[i][1] =  math.log10(ylimits[i][1])  
  
    # loop over all input
    for k in range(0,len(inputFiles)):
    
        fileName = inputFiles[k].split('/')[-1]
        data = np.transpose(np.genfromtxt(inputFiles[k],dtype='float'))
        
        if( len(legends) == len(inputFiles) ):
            tmpLegend = legends[k]
        else:
            tmpLegend = fileName
        
        if k>= stylechange:
            if doColor:
                lineColor = color[k-stylechange]
            else:
                lineColor = 'black'
            linestyle = '--'
        else:
            if doColor:
                lineColor = color[k]
            else:
                lineColor = 'black'
            linestyle = '-'
          
         # if color or ls were set by user, override  the above:
        if ls:
             linestyle = ls[k]
         
        if colors:
             lineColor = color[k]
    

    
        for i in range(len(data[0])-1,-1,-1):
            if np.isnan(data[9,i]):
                data[2:10,i] = data[2:10,i+1]
        
            #print len(data[0])
            #for i in range(len(data[0])-1,-1,-1):
            #if np.any(data[2:10,i]<0):
            ##print i
            #data = np.delete(data, [i], axis=1)
            ##np.delete(data[:,i])
            #print len(data[0])

        if logT == True:
            for i in range(0,len(data[0])):
                data[7,i] = math.log10(data[7,i])
                data[8,i] = math.log10(data[8,i])
                #data[9,i] = math.log10(data[9,i])
            
        if logFractions == True:
            #axarr[0,0].set_yscale("log")
            #axarr[0,1].set_yscale("log")
            #axarr[1,0].set_yscale("log")
            #axarr[1,1].set_yscale("log")           
            #axarr[2,0].set_yscale("log")
            for ii in range(2,7):
              for i in range(0,len(data[0])):
                  value = math.fabs(data[ii,i])
                  if value == 0.0:
                      value = 1.0e-15                       
                  data[ii,i] = math.log10(value)    # some of Rajat's old profiles have negative fractions
            
        # plot stuff
        for i in range(0,rows):
            for j in range(0,columns):
                
                if linestyle == '-':
                    axarr[i,j].plot( data[1], data[(i+1)*columns+j], linestyle, lw=1.8, color=lineColor, label=r'$\mathrm{%s}$'%(tmpLegend) )

                if linestyle=='--':    
                    axarr[i,j].plot( data[1], data[(i+1)*columns+j], linestyle, lw=1.8, dashes=(9, 2), color=lineColor, label=r'$\mathrm{%s}$'%(tmpLegend) )
                    
                if linestyle=='-.':
                    axarr[i,j].plot( data[1], data[(i+1)*columns+j], linestyle, lw=1.8, dashes=(1, 2, 8, 2), color=lineColor, label=r'$\mathrm{%s}$'%(tmpLegend) )



    
        for i in range(0,rows):
            for j in range(0,columns):
                #print i*columns+j
                axarr[i, j].set_xlim(xLimits[0],xLimits[1])
                if ylimits[i*columns+j][0] == None:
                    ymin = axarr[i, j].get_ylim()[0]
                else:
                    ymin = ylimits[i*columns+j][0]
                if ylimits[i*columns+j][1] == None:
                    ymax = axarr[i, j].get_ylim()[1]
                else:
                    ymax = ylimits[i*columns+j][1]            
                
                axarr[i, j].set_ylim(ymin,ymax)
                if i == rows-1:
                    axarr[i, j].set_xlabel(xlabel)#,fontsize=20)
                    axarr[i, j].xaxis.label.set_size(labelFontSize)
                if j == columns-1:
                    axarr[i, j].yaxis.tick_right()
                    axarr[i, j].yaxis.set_ticks_position('both')
                    axarr[i, j].yaxis.set_label_position("right")
                    axarr[i, j].set_ylabel(ylabel[i*columns+j])
                    axarr[i, j].yaxis.label.set_size(labelFontSize)
                else:
                    axarr[i, j].yaxis.tick_left()
                    axarr[i, j].yaxis.set_ticks_position('both')
                    axarr[i, j].yaxis.set_label_position("left")
                    axarr[i, j].set_ylabel(ylabel[i*columns+j])
                    axarr[i, j].yaxis.label.set_size(labelFontSize)
                
      
    f.subplots_adjust(hspace=plotdistY)
    f.subplots_adjust(wspace=plotdistX)
    if showLegend:
        lgd = plt.legend(loc='upper center', bbox_to_anchor=(-0.1, -0.20), ncol=ncolsLegend, fontsize=10)#, fancybox=True, shadow=True)
        #plt.minorticks_on()
        #plt.savefig('example_plot.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.savefig(outFile, bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
        plt.savefig(outFile, bbox_inches='tight')
