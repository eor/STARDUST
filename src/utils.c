 
/***************************************************************
 * libraries
 ***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/***************************************************************
 * SD headers
 ***************************************************************/
#include "constants.h"
#include "allvars.h"
#include "prototype.h"
#include "table_settings.h"
#include "log.h"

/***************************************************************
 * Etc
 ***************************************************************/
#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))

/***************************************************************
 * given directory and filename this function returns a full path
 ***************************************************************/
char* utils_concat_path(char *dir, char *file){
    
    // before the file we add a prefix
    // also add "_M%.3f_z%.3f" --> "M8.123_z9.123_"
    
    int lenPrefix = 14;
    if (myConfig.haloMass>=10.0) lenPrefix +=1;
    if (myConfig.redshiftLow>=10)lenPrefix +=1;
    
    
    char *tmpPrefix = (char*) malloc( sizeof(char) * lenPrefix);
    sprintf(tmpPrefix, "M%.3f_z%.3f_", myConfig.haloMass, myConfig.redshiftLow );
    
    
    char *result = (char*) malloc( sizeof(char) * (strlen(dir)+1+strlen(file)+1 +lenPrefix) ); // +1 for "/" and +1 for the zero-terminator
    
    strcpy(result, dir);
    strcat(result, "/");
    strcat(result,  tmpPrefix);
    strcat(result, file);
    
    free(tmpPrefix);   
    
    return result;
    
}


/***************************************************************
 * given directory and filename this function returns a full path
 ***************************************************************/
char* utils_concat_path_noID(char *dir, char *file){
    
    
    char *result = (char*) malloc( sizeof(char) * (strlen(dir)+1+strlen(file)+1 ) ); // +1 for "/" and +1 for the zero-terminator
    
    strcpy(result, dir);
    strcat(result, "/");
    
    strcat(result, file);
    
    return result;    

}


/***************************************************************
 * print help, then exit
 ***************************************************************/
void utils_print_help(){
    
    printf(" \n");
    printf("STARDUST help page\n-------------------------\n");
    printf("\n Usage:\n\t ./STARDUST <path/to/config file>\n\n");
    exit(0);    
    
}

/***************************************************************
 * print hello, 
 ***************************************************************/
void utils_print_hello(){

// print a friendly hello message at startup   
//   ___ _____ _   ___ ___  _   _ ___ _____ 
//  / __|_   _/_\ | _ \   \| | | / __|_   _|
//  \__ \ | |/ _ \|   / |) | |_| \__ \ | |  
//  |___/ |_/_/ \_\_|_\___/ \___/|___/ |_|  

    
    printf(" \n\n");
    printf(" Welcome to \n");     
    printf("  ___ _____ _   ___ ___  _   _ ___ _____      \n");
    printf(" / __|_   _/_\\ | _ \\   \\| | | / __|_   _|  \n");
    printf(" \\__ \\ | |/ _ \\|   / |) | |_| \\__ \\ | |  \n");
    printf(" |___/ |_/_/ \\_\\_|_\\___/ \\___/|___/ |_| (version %.3f)\n", SD_VERSION);
    printf(" -------------------------------------------------------\n\n");
 
    
}
