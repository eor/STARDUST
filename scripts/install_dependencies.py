#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import subprocess as sub
import sys

#------------------------------------------------------------------------------
# User defined settings
#------------------------------------------------------------------------------
#
# 1. Set install path, select *one*
#

# a) location given relative to  home dir
installDir      = os.path.expanduser('~/local')        

# b) location given relative to exec path of this script
#installDir      = os.path.abspath("../../local1") 

# c) absolute path
#installDir      = os.path.abspath("/home/tester/Software/local")

#
# 2. Select packages to install
#

installLibConfig = True
installGSL       = True
installBoost     = True


#------------------------------------------------------------------------------
# Global variables
#------------------------------------------------------------------------------

softwareDir     = installDir+'/dl' 
pythonVersion   = '.'.join(sys.version.split('.')[0:2])

baseURL         = 'http://www.astro.rug.nl/~krause/static/STARDUST/'
packageList     = [
                    'libconfig-1.5.tar.gz',         
                    'gsl-1.16.tar.gz', 
                    'boost_1_59_0_redux.tar.gz'
                  ]

installPackage = [installLibConfig, installGSL, installBoost]
    
    
#------------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------------
def check_all_dependencies():
    print("Checking for dependencies")
    for dep in ["gcc", "make", "wget"]:
        if not check_dependency_tool(dep):
            print("\nERROR: %s not found. Exiting."%dep )
            exit(1)

def check_dependency_tool(toolName):
    print("Checking for %s"%toolName)
    try:
        f = open(os.devnull, "w")
        sub.call(toolName, stdout=f, stderr=f)
        return True
    except OSError:
        return False
    
def check_dependency_lib(libName):
    print( "Checking for %s"%libName )
    # TODO ... what do we need?

def create_directories(path):
    # function creates path directory and all necessary sub-directories
    for d in [ 'bin', 'data', 'dl', 'include', 'lib', 'man', 'python', 'share', 'src', 'python/lib/python%s/site-packages'%(pythonVersion), 'python/lib64/python%s/site-packages'%(pythonVersion)]:
        tmpDir = path+'/'+d
        if not os.path.exists(tmpDir):
            print( ' Creating %s '%tmpDir )
            os.makedirs(tmpDir)
            
def run_on_shell(cmd):
    print("> " + cmd)
    os.system(cmd)    
    
def download_packages():
    print( 'Checking for missing packages' )
    for i in range(len(packageList)):
        if installPackage[i]:
            tmpDir = softwareDir + '/' + packageList[i]
            if not os.path.exists(tmpDir):
                print( ' Downloading %s:'%packageList[i] )
                tmpURL = baseURL + '/' + packageList[i]
                run_on_shell('cd %s && wget -c %s '%(softwareDir, tmpURL))
                # check if download was successful
                if not os.path.exists(tmpDir):
                    print("\nERROR: Could not download file '%s'. Exiting."%packageList[i])
                    exit(1)
    

def install_typical_autoconf(srcTarFile, srcDir, configureOptions="", postInstallCommand=None):
    print('Installing %s'%srcDir)
    srcTarPath  = softwareDir + '/' + srcTarFile    
    run_on_shell('tar -xzvf %s -C ${SD}/src'%( srcTarPath))
    run_on_shell('cd ${SD}/src/%s/ && ./configure %s --prefix=${SD} && make && make install'%(srcDir, configureOptions))
    if postInstallCommand is not None:
        run_on_shell('cd %s/src/%s/ && %s'%(installDir, srcDir, postInstallCommand))

    
def install_tools(srcTarFile, srcDir):
    toolList = ['make_dummy_glass', 'mean_X_HI_mass_weighted', 'mean_X_HI', 'powerspec']         
    print('Installing various small tools')
    srcTarPath  = softwareDir + '/' + srcTarFile 
    run_on_shell('tar -xzvf %s -C %s/src/'%( srcTarPath, installDir))
    run_on_shell('cd %s/src/%s/ && make'%(installDir, srcDir))
    for tool in toolList:
        run_on_shell('cd %s/bin/ && ln -s %s/src/%s/%s'%(installDir, installDir, srcDir, tool ) )   
    

def install_unpack_only(srcTarFile, srcDir, pkgName):
    print('Installing %s'%pkgName)
    srcTarPath  = softwareDir + '/' + srcTarFile  
    run_on_shell('tar -xzvf %s -C %s/src/'%( srcTarPath, installDir))
    
    
def install_boost(srcTarFile, srcDir, pkgName):
    print('Installing %s'%pkgName)
    srcTarPath  = softwareDir + '/' + srcTarFile  
    run_on_shell('tar -xzvf %s -C %s/'%( srcTarPath, installDir))

#def install_python_module(srcTarFile,srcDir,pkgName):
    #print('Installing python module %s'%pkgName)
    #srcTarPath  = softwareDir + '/' + srcTarFile
    #run_on_shell('tar -xzvf %s -C %s/src/'%( srcTarPath, installDir))
    #run_on_shell('cd %s/src/%s/ && python setup.py build'%( installDir, srcDir))
    #run_on_shell('cd %s/src/%s/ && python setup.py install --prefix=%s/python/'%( installDir, srcDir, installDir ))
    
#def set_up_configs(srcTarFile='configs.tar.gz'):
    #print('Copying configs and Makefiles')
    #srcTarPath  = softwareDir + '/' + srcTarFile
    #run_on_shell('tar -xzvf %s -C %s/configs/'%( srcTarPath, installDir))   
 
 
def set_paths():    
    os.environ["PATH"] = os.environ.get("PATH") + ":%s/bin"%installDir
    os.environ["PATH"] = os.environ.get("PATH") + ":%s/sbin"%installDir
    if os.environ.get("LD_LIBRARY_PATH") != None:
        os.environ["LD_LIBRARY_PATH"] = os.environ.get("LD_LIBRARY_PATH", "") + ":%s/lib/"%installDir
    else:
        os.environ["LD_LIBRARY_PATH"] = "%s/lib/"%installDir
    
    os.environ["SD"] = "%s"%installDir
    if os.environ.get("PYTHONPATH") != None:
        os.environ["PYTHONPATH"] = '%s/python/lib/python%s/site-packages:%s/python/lib64/python%s/site-packages:'%( installDir, pythonVersion, installDir, pythonVersion  ) + os.environ.get("PYTHONPATH")
    else:
        os.environ["PYTHONPATH"] = '%s/python/lib/python%s/site-packages:%s/python/lib64/python%s/site-packages'%( installDir, pythonVersion, installDir, pythonVersion  ) 
    #os.environ["HDF5_DIR"] = "%s"%installDir
    
    
def user_message():    
    print('\n All done here.\n\n')
    
    print(' To compile STARDUST edit its Makefile and set the path variables for the libraries accordingly: \n')
    
    if installPackage[1]:
        print('\tGSL_INCL\t= -I%s/include'%installDir)
        print('\tGSL_LIB \t= -L%s/lib\n'%installDir)
    if installPackage[0]:
        print('\tCONF_INCL\t= -I%s/include'%installDir)
        print('\tCONF_LIB \t= -L%s/lib\n'%installDir)
    if installPackage[2]:
        print('\tBOOST_INCL\t= -I%s/boost_1_59_0_redux\n'%installDir)
    
    
    #print(' If you wish to use these libraries with other programs as well, you might want to consider adding')
    #print(' the necessary environment variables permanently.\n')
    #print('\tTo do so, edit your .bashrc/.cshrc/.tcshrc and add the following:\n')
    
    #print('\tTo the PATH variable:')
    #print('\t\t%s/bin'%installDir)
    #print('\t\t%s/sbin\n'%installDir)
    
    #print('\tTo the LD_LIBRARY_PATH variable:')
    #if installPackage[0] or installPackage[1]:
    #    print('\t\t%s/lib'%installDir)
    #if installPackage[1]:  
    #    print('\t\t%s/boost_1_59_0_redux\n'%installDir)
    
    #print('\tTo the PYTHONPATH variable:')
    #if os.listdir('%s/python/lib/python%s/site-packages'%(installDir, pythonVersion)):
    #  print('\t\t%s/python/lib/python%s/site-packages\n'%(installDir, pythonVersion))
    #if os.listdir('%s/python/lib64/python%s/site-packages'%(installDir, pythonVersion)):
    #  print('\t\t%s/python/lib64/python%s/site-packages\n'%(installDir, pythonVersion))
      
      
    #print('\t\tTODO!!\n')    
    #export PYTHONPATH="${HOME}/local/lib/python2.6/site-packages:${PYTHONPATH}"
    #export PYTHONPATH="${HOME}/local/lib64/python2.6/site-packages:${PYTHONPATH}"

    #export PYTHONPATH="/net/hathor/data4/users/krause/pyBEARS/python:${PYTHONPATH}"
  
    #print('\tAdd the SD variable. E.g. for Bash:')
    #print("\t\texport SD='%s'\n"%(installDir))
    

#------------------------------------------------------------------------------
# Main script starts here
#------------------------------------------------------------------------------
if __name__ == "__main__":

    # 0. SETUP
    check_all_dependencies()
    create_directories(installDir)
    set_paths()                     # sets environment variables  
    download_packages()

    # 1. AUTOCONF COMPONENTS 
    print("Handling autoconf components")
    
    # gsl
    if not os.path.exists( '%s/lib/libgsl.a'%installDir ):
        install_typical_autoconf(srcTarFile='gsl-1.16.tar.gz', srcDir='gsl-1.16')

    # libconfig
    if not os.path.exists( '%s/lib/libconfig.a'%installDir ):
        install_typical_autoconf(srcTarFile='libconfig-1.5.tar.gz', srcDir='libconfig-1.5')

    # 2. OTHER COMPONENTS 
    print("Handling other components")
 
    # boost_1_59_0_redux
    if not os.path.exists( '%s/boost_1_59_0_redux'%installDir ):
        install_boost(srcTarFile='boost_1_59_0_redux.tar.gz', srcDir='../boost_1_59_0_redux', pkgName='Boost (reduced)')


    # 3. FINISH
    user_message()




