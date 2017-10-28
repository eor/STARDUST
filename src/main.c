
/***************************************************************
 * This program simulates the evolution of HI, HII, HeI, HeII, 
 * HeIII fractions and the kinetic (gas), spin, and brightness 
 * temperatures around a central source (defined by some 
 * input SED).  
 * 
 * Originally written by Rajat Mani Thomas (2005). 
 * Currently maintained by someone else.
 * 
 ***************************************************************/


/***************************************************************
 * libraries
 ***************************************************************/
#include <stdio.h>

#include "prototype.h"
#include "allvars.h"
#include "log.h"

/***************************************************************
 * start of main - fasten seatbelts!
 ***************************************************************/
int main(int argc, char **argv){
    
    if(argc!=2){ utils_print_help();}     

    utils_print_hello();
    
   /***************************************************************
    * read config file & start logging
    ***************************************************************/   
   
    int  fileNameLength = strlen(argv[1]);
    char fileName[fileNameLength];
    
    strncpy(fileName,argv[1],fileNameLength); 
    fileName[fileNameLength] = '\0';  
    
    config_load_from_file(fileName);     
    
    log_start();
    
    /***************************************************************
     * read SED         
     ***************************************************************/

    sed_read_file(myConfig.pathSED);
    
    // TODO: read UVB spectrum file     
    // TODO: read background density file

    /***************************************************************
     * pre-compute integrals         
     ***************************************************************/    

     table_ion();
     table_temp();
     table_compton(); 
     
     read_tables();

    /***************************************************************
     * main run happens here
     ***************************************************************/    
          
    rt_main_run();            

    /***************************************************************
     * clean up
     ***************************************************************/   
    
    memory_free_all();
    
    log_close();
    
    return 0;    
        
}


