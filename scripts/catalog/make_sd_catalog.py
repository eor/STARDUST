#!/usr/bin/env python
# -*- coding: utf-8 -*-
#-------------------------------------------------------------
# This file contains functions to set up and run
# STARDUST to generate spectral catalogs for BEARS
#
# catalogs for the following source types are currently available
#   (1) A halo of given mass populated with a pop III stars of the same mass 
#   (2) QSO-like (power law) for a given host halo mass,
#   (3) A halo of given mass populated with pop III stars following an IMF
#   (4) A combined SED from (2) and (3)
#
#-------------------------------------------------------------

import numpy as np
import subprocess as sub
import multiprocessing
import os, string, re, sys
import random
from functools import partial

sys.path.append('../sed_generator/')
import sed

#-------------------------------------------------------------
# Default values & parameters
#-------------------------------------------------------------

# number of CPUs to be used to run SD
nCPU                    = 2            

# script related
sdExe                   = 'STARDUST'
sdExePath               = '/path/to/your/exe_file/STARDUST'
sdConfigTemplate        = './SD_CONFIG_TEMPLATE' 
logFile                 = 'log_catalog' 


# SED related
sourceTypes             = ['BB','PL','IMF','PL+IMF'] # (1), (2), (3), (4)

# Note: the following are not yet processed by the script
# If you want to change those, add the functionality to the script or
# change the value in the SD_CONFIG_TEMPLATE file.

eHighDefault            = 1.0e4          # eV
eLowDefault             = 13.6

# SD related            
sourceLifetimeDefault   = 10.0           # Myr
settingsDebugDefault    = 0
settingsRStartDefault   = 0.01           # IGM is considered ionized from 0 to this radius in kpc 
settingsRMaxDefault     = 1500.0         # Maximum radius radially away from source in kpc
settingsDeltaRDefault   = 0.1            # Spatial resolution of the simulation in kpc
settingsDeltaTDefault   = 0.01           # Time resolution in the simulation 
settingsWriteTDefault   = 1.0


# Cosmology
cosmoOmegaM             = 0.3089          # Total matter density 
cosmoOmegaL             = 0.6911          # Dark energy / cosmological constant 
cosmoOmegaB             = 0.0486          # Baryon density
cosmoH0                 = 67.74           # Hubble parameter at z = 0
cosmoH100               = 0.6774          # Hubble / 100
cosmoSigma8             = 0.8159          # power spectrum normalization    
cosmoTauThom            = 0.078           # Optical depth of Thomson scattering
cosmoTCMB               = 2.731           # CMB temperature [K] at z = 0



#-------------------------------------------------------------
# Functions
#-------------------------------------------------------------

def set_up_sd_dir(directory):
  
  # check if directory already exists
  if os.path.exists(directory):
      
    ## if it exists rename it with the current time stamp
    ## first, get the time
    from time import strftime, gmtime
    timestamp = strftime("%Y%m%d_%H:%M", gmtime())
    
    ## rename the directory accordingly
    if not os.path.exists(directory+'_backup_'+timestamp):
      os.rename(directory,directory+'_backup_'+timestamp)
    else:
      print 'Error: Backup directory already exists'
      exit(1)
      
    ## create the directory new
    os.makedirs(directory)
    print 'Renamed '+directory+' to '+directory+'_backup_'+timestamp
    
  # if the directory does not exist yet
  else:
    ## create it
    os.makedirs(directory)
    print 'Created empty directory '+directory



def build_catalog(directory, Ncpu, zList, haloMassList, name='default', **kwargs):

    #-1. parameters is an object/ dictionary that needs to evaluated before the main work can start    
    sourceType = kwargs['source']
    del kwargs['source']
  
    # 0. set up output dir, move old stuff to ackup dir, if necessary.
    set_up_sd_dir(directory+name)
    
    # 1. 
    configFileList = []
    
    # 2. for each grid point
    for redshift in zList:
        for logMhalo in haloMassList:
            
            ## generate the SED
            if sourceType == 'BB'    : sed.generate_SED_stars_BB(  pow(10,logMhalo), redshift, fileName = directory+name+'/sed_%s_M%.3f_z%.3f.dat'%( sourceType, logMhalo, redshift), **kwargs)
            if sourceType == 'PL'    : sed.generate_SED_PL(        pow(10,logMhalo),           fileName = directory+name+'/sed_%s_M%.3f_z%.3f.dat'%( sourceType, logMhalo, redshift), **kwargs)
            if sourceType == 'IMF'   : sed.generate_SED_stars_IMF( pow(10,logMhalo), redshift, fileName = directory+name+'/sed_%s_M%.3f_z%.3f.dat'%( sourceType, logMhalo, redshift), **kwargs)
            if sourceType == 'IMF+PL': sed.generate_SED_IMF_PL(    pow(10,logMhalo), redshift, fileName = directory+name+'/sed_%s_M%.3f_z%.3f.dat'%( sourceType, logMhalo, redshift), **kwargs)
            
            ## generate the config
            ### open and read template config file
            #f = open('./SD_CONFIG_TEMPLATE')
            #sedTemplate = f.read()
            #f.close()
            ### replace the sed file name
            #sedTemplateNew = sedTemplate.replace()
            #sedTemplateNew = sedTemplateNew.replace('P.ID','id_M%.3f_z%.3f'%(logMhalo, redshift))
            #sedTemplateNew = sedTemplateNew.replace('P.haloMass','%f'%(logMhalo))
            #sedTemplateNew = sedTemplateNew.replace('P.zLow','%f'%(redshift))
            #sedTemplateNew = sedTemplateNew.replace('P.zHigh','%f'%(redshift))
            
            configNew = []
            with open(sdConfigTemplate,'r') as f:                
                for line in f.readlines():    
                    line = line.replace('P.pathSED'  , 'sed_%s_M%.3f_z%.3f.dat'%(sourceType, logMhalo, redshift) )       
                    #line = line.replace('P.ID'       ,'id_M%.3f_z%.3f'%(logMhalo, redshift))
                    line = line.replace('P.haloMass' , '%f'%(logMhalo)) 
                    line = line.replace('P.zLow'     , '%f'%(redshift))
                    line = line.replace('P.zHigh'    , '%f'%(redshift))
                    configNew.append(line)
            
            
            
            ### write out config file
            configTarget =  directory+name+'/config_%s_M%.3f_z%.3f.dat'%(sourceType, logMhalo, redshift)
            
            with open(configTarget, 'w') as f:
                for line in configNew:
                    f.write(line)
            
            
            
            ### append config file name to array
            configFileList.append('config_%s_M%.3f_z%.3f.dat'%(sourceType, logMhalo, redshift))




    # 3. copy sdExe to catalog dir, set up run commands, see # 4. in old_SD_run_skript
    sub.call(' cp %s %s '%(sdExePath, directory+name), shell=True)
    
    # 4. this is a workaround to be able to pass more (non-iterable) arguments to the run function
    run_SD_additional_arguments = partial(run_SD,directory+name,configFileList)
        
    # run stardust in multiprocessing with Ncpu cores
    pool = multiprocessing.Pool(processes=Ncpu)
    pool.map(run_SD_additional_arguments, range(0, len(configFileList)))

    
    # 4. after finishing, write metadata to some log file in the catalog containing all parameters, date, total runtime, plus whatever else cold be useful
    # tba

    # 5. delete tables
    cleanup_SD_dir(directory+name)
  


def run_SD(OutDir, fileList, i):
     
    wait = random.randint(0, nCPU*4) 
    # SD is I/O heavy in the beginning, which causes high load when launching lots of instaces at once, 
    # therefor a random wait to spread out the load ... until we fix the table issue in SD ¯\_(ツ)_/¯
  
    sub.call('sleep %d && cd %s && ./%s  %s '%(wait, OutDir, sdExe, fileList[i]), shell=True)
      
      
      
def drange(start, stop, step):
    while start <= stop:
        yield start
        start += step
        
        
        
        
def cleanup_SD_dir(outDir):
        print 'deleting table files'
        sub.call(' cd %s && rm *table_i* '%( outDir ), shell=True)
        sub.call(' cd %s && rm *table_c* '%( outDir ), shell=True)
        sub.call(' cd %s && rm *table_t* '%( outDir ), shell=True)
    


#-----------------------------------------------------------------  
# Testing 
#-----------------------------------------------------------------
if __name__ == "__main__":
 

    myRedshiftList  = []
    myMassList      = []
 
    # desired redshift range (min, max, step)
    for z in drange(6.0, 8.0, 0.5):
        myRedshiftList.append(z)

    # desired log_10 (M_halo/M_sol) range (min, max, step)
    for m in drange(9, 13, 0.5):
        myMassList.append(m)
        

    # example for a power-law type source catalog    
    #build_catalog('./', nCPU, zList=myRedshiftList, haloMassList=myMassList, name='myCatalog_PL', eHigh=eHighDefault, eLow=eLowDefault, source='PL', alpha=1.5)    
    

    # example black body type catalog
    #build_catalog('./', nCPU, zList=myRedshiftList, haloMassList=myMassList, name='myCatalog_BB_200',     eHigh=eHighDefault, eLow=eLowDefault, source='BB', starMass=200.0, fEsc=.10)
    

    # example IMF source catalog
    #build_catalog('./', nCPU, zList=myRedshiftList, haloMassList=myMassList, name='IMF_fEsc0.10', eHigh=eHighDefault, eLow=eLowDefault, source='IMF', starMassMin=5.0, starMassMax=100.0, fEsc=.10,logGrid=True)
   
    # example for a co-evolution catalog
    #build_catalog('./', nCPU, zList=myRedshiftList, haloMassList=myMassList, name='PL_alpha1.5+IMF_redux_fEsc0.05', eHigh=eHighDefault, eLow=eLowDefault, source='IMF+PL', alpha=1.5, starMassMin=5.0, starMassMax=100.0, fEsc=0.05, logGrid=True)
























































