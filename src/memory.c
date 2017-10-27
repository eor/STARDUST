/*! File allocate.c allocates and determines all the memory required in the 
program. It also frees the memory once the program is about to terminate.

*/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <stdlib.h>

#include "table_settings.h"
#include "prototype.h"
#include "allvars.h"


#define pow2(x) ((x)*(x))

void memory_allocate_all(void){

    long int totalMemoryAlloc = 0;
    int i;


    ne           = (double*) malloc(sizeof(double) * numGridPoints);    
    n_H1         = (double*) malloc(sizeof(double) * numGridPoints);    
    n_He1        = (double*) malloc(sizeof(double) * numGridPoints);    
    n_H2         = (double*) malloc(sizeof(double) * numGridPoints);    
    n_He2        = (double*) malloc(sizeof(double) * numGridPoints);    
    n_He3        = (double*) malloc(sizeof(double) * numGridPoints);    
    T_e          = (double*) malloc(sizeof(double) * numGridPoints);    
    T_spin       = (double*) malloc(sizeof(double) * numGridPoints);    
    T_brig       = (double*) malloc(sizeof(double) * numGridPoints);
    
    x_HI         = (double*) malloc(sizeof(double) * numGridPoints);        
    x_HII        = (double*) malloc(sizeof(double) * numGridPoints);
    x_HeI        = (double*) malloc(sizeof(double) * numGridPoints);
    x_HeII       = (double*) malloc(sizeof(double) * numGridPoints);
    x_HeIII      = (double*) malloc(sizeof(double) * numGridPoints);
    

    fuku_e1h1    = (double*) malloc(sizeof(double) * numGridPoints); 
    fuku_ehe1    = (double*) malloc(sizeof(double) * numGridPoints); 
    fuku_ehe2    = (double*) malloc(sizeof(double) * numGridPoints);    

    integral_H1  = (double*) malloc(sizeof(double) * numGridPoints); 
    integral_He1 = (double*) malloc(sizeof(double) * numGridPoints); 
    integral_He2 = (double*) malloc(sizeof(double) * numGridPoints); 

    comp_integ1  = (double*) malloc(sizeof(double) * numGridPoints); 
    comp_integ2  = (double*) malloc(sizeof(double) * numGridPoints);
    

    xa           = (double*) malloc(sizeof(double) * INTERPOINTS); 
    nHx1a        = (double*) malloc(sizeof(double) * INTERPOINTS); 
    nHex2a       = (double*) malloc(sizeof(double) * INTERPOINTS); 
    nHx13a       = (double*) malloc(sizeof(double) * INTERPOINTS); 
    ya           = (double*) malloc(sizeof(double) * INTERPOINTS);
    

    yP2a = (double**) malloc(sizeof(double*) * INTERPOINTS);
    for(i=0;i<INTERPOINTS;i++)
        yP2a[i] = (double*) malloc(sizeof(double) * INTERPOINTS);


    yP3a = (double**) malloc(sizeof(double*) * INTERPOINTS);
    for(i=0;i<INTERPOINTS;i++)
        yP3a[i] = (double*) malloc(sizeof(double) * INTERPOINTS);
    
    outFile = (char*) malloc( sizeof(char) * 256);
    
    
    totalMemoryAlloc += 22 * sizeof(double) * numGridPoints;  
    totalMemoryAlloc += 5  * sizeof(double) * INTERPOINTS;   
    totalMemoryAlloc += 2  * sizeof(double) * INTERPOINTS * INTERPOINTS;
    totalMemoryAlloc +=      sizeof(char)   * 256;

#ifdef ODEJENSSTYLE
  mem_jens = new void*[numGridPoints];
  for(i=0; i<numGridPoints;i++) 
      mem_jens[i] = NULL;
#endif

  printf(" Total allocated memory: %.2f KB\n",  (double) ( totalMemoryAlloc / (1024.0) )  );

}

/***************************************************************
 * Freeing all the memory allocated above
 ***************************************************************/
void memory_free_all(void){

    int i;


    free(xa);
    free(nHx1a);
    free(nHex2a);
    free(nHx13a);
    free(ya);

    for(i = 0; i < INTERPOINTS; i++)
        free(yP2a[i]);
    free(yP2a);

    for(i = 0; i < INTERPOINTS; i++)
        free(yP3a[i]);
    free(yP3a);

    free(ne);
    free(n_H1);
    free(n_He1);
    free(n_H2);
    free(n_He2);
    free(n_He3);
    free(T_e);
    free(T_brig);
    free(T_spin);
    
    free(x_HI); 
    free(x_HII);  
    free(x_HeI); 
    free(x_HeII); 
    free(x_HeIII); 

    free(fuku_e1h1);
    free(fuku_ehe1);
    free(fuku_ehe2);

    free(integral_H1);
    free(integral_He1);
    free(integral_He2);

    free(comp_integ1);
    free(comp_integ2);

    free(Lambda);
    free(Energy);

    free(outFile);

}
