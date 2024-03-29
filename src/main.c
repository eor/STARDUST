
/***************************************************************
 * This program simulates the evolution of HI, HII, HeI, HeII, 
 * HeIII fractions and the kinetic (gas), spin, and brightness 
 * temperatures around a central source (defined by some 
 * input SED).  
 * 
 * Originally written by Rajat Mani Thomas (2005). 
 * Currently maintained by Fabian Krause.
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
   
    size_t  len = strlen(argv[1]);
    char    file_name[len];
    
    strncpy(file_name, argv[1], len);
    file_name[len] = '\0';
    
    config_load_from_file(file_name);
    config_set_grid_points();
    log_start();


    /***************************************************************
     * allocate required memory for the computing grids
     ***************************************************************/

    memory_allocate_all();
    
    /***************************************************************
     * read SED & over densities
     ***************************************************************/

    sed_read_file(myConfig.pathSED);
    density_read_file(myConfig.pathDensity);

     // TODO: read UVB spectrum file

    /***************************************************************
     * pre-compute integrals         
     ***************************************************************/    

     table_ion();
     table_temp();
     table_compton(); 
     
     table_read_all();

    /***************************************************************
     * main run happens here
     ***************************************************************/    

    rt_main_run();

    /***************************************************************
     * clean up
     ***************************************************************/   
    
    memory_free_all();
    
    log_close();
    
    printf("\n Run completed successfully. Exiting.\n\n");
    
    return 0;    
        
}


