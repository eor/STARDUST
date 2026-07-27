 
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
 * Build "dir/file", or "dir/pathID_file" path if
 * myConfig.pathID is set.
 *
 * Caller owns the returned buffer and must free() it.
 ***************************************************************/
char* utils_concat_path(const char *dir, const char *file){

    int havePathID = (myConfig.pathID[0] != '\0');

    size_t len = strlen(dir) + 1 /* "/" */ + strlen(file) + 1 /* '\0' */;
    if (havePathID){
        len += strlen(myConfig.pathID) + 1; /* pathID + "_" */
    }

    char *result = (char*) malloc(sizeof(char) * len);
    if (!result){
        printf("ERROR: utils_concat_path: malloc failed for '%s/%s'. Exiting.\n", dir, file);
        exit(1);
    }

    strcpy(result, dir);
    strcat(result, "/");
    if (havePathID){
        strcat(result, myConfig.pathID);
        strcat(result, "_");
    }
    strcat(result, file);

    return result;
}


/***************************************************************
 * Build "dir/file" path, ignoring myConfig.pathID entirely.
 *
 * Kept only for the existing calls in rt.c that explicitly want
 * the no-pathID behaviour regardless of whether pathID is set. New
 * code should just call utils_concat_path(), which now does the
 * right thing in both cases.
 * TODO: remove!
 ***************************************************************/
char* utils_concat_path_noID(const char *dir, const char *file){

    size_t len = strlen(dir) + 1 /* "/" */ + strlen(file) + 1 /* '\0' */;

    char *result = (char*) malloc(sizeof(char) * len);
    if (!result){
        printf("ERROR: utils_concat_path_noID: malloc failed for '%s/%s'. Exiting.\n", dir, file);
        exit(1);
    }

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
