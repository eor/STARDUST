 
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
 * given directory and filename this function returns a full path
 ***************************************************************/
char* utils_concat_path(char *dir, char *file){
    
    // we insert the config parameter given in pathID in front of
    // the filename.
    
    int lenPrefix = strlen(myConfig.pathID);
    
    char *result = (char*) malloc( sizeof(char) * ( strlen(dir) + 1 + lenPrefix +1 + strlen(file) + 1 ) );
    // +1 for "/" 
    // +1 for "_"
    // +1 for the zero-terminator
    
    strcpy(result, dir);
    strcat(result, "/");
    strcat(result,  myConfig.pathID);
    strcat(result, "_");
    strcat(result, file);

    return result;
    
}


/***************************************************************
 * given directory and filename this function returns a full path
 ***************************************************************/
char* utils_concat_path_noID(char *dir, char *file){


    char *result = (char*) malloc( sizeof(char) * ( strlen(dir) + 1 + strlen(file) + 1 ) );
    // +1 for "/" and +1 for the zero-terminator
    
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
