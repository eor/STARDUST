/***************************************************************
 * Libraries and header files
 ***************************************************************/
#include <stdio.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "allvars.h"
#include "prototype.h"
#include "table_settings.h"
#include "log.h"

void table_read_all(void){
 
    /***************************************************************
    * NOTE: 
    * The tables generated by table_ion.c, table_temp.c and 
    * table_compton.c are read in the following lines. The tables
    * are stored in the following format:
    *     column 1: n_HeI fractions i.e., n_He times x_HeI 
    *     column 2: n_HeII times x_HeII +n_HI times x_HI
    * 
    * Column 2 could be combined in such a way only becaue the
    * cross sections of n_HI and n_HeII have the same function form
    * except that they are scaled by a factor of 64.        
    ***************************************************************/    

    int i,j; 
    long int array_len_x1a, array_len_x2a, array_len_x13a; 
    
    array_len_x1a  = (long int)( (UPLIM1-LOWLIM)/TABLERES+1 ); 
    array_len_x2a  = (long int)( (UPLIM1-LOWLIM)/TABLERES+1 ); 
    array_len_x13a = (long int)( (UPLIM2-LOWLIM)/TABLERES+1 );
  
    FILE *fp_e1h1_p1,*fp_e1h1_p2,*fp_e1h1_p3;
    FILE *fp_ehe1_p2,*fp_ehe1_p3,*fp_ehe2_p3;

    FILE *fp_temp_e1h1_p1,*fp_temp_e1h1_p2,*fp_temp_e1h1_p3;
    FILE *fp_temp_ehe1_p2,*fp_temp_ehe1_p3;
    FILE *fp_temp_ehe2_p3;

    FILE *fp_comp1_p1,*fp_comp1_p2,*fp_comp1_p3;
    FILE *fp_comp2_p1,*fp_comp2_p2,*fp_comp2_p3;

    if(DEBUG)log_debug("Reading integral files");
    
    
    fp_e1h1_p1 = fopen( utils_concat_path( myConfig.pathOutDir, (char *) "table_ion_e1h1_p1.dat" ), "r" );
    fp_e1h1_p2 = fopen( utils_concat_path( myConfig.pathOutDir, (char *) "table_ion_e1h1_p2.dat" ), "r" );
    fp_e1h1_p3 = fopen( utils_concat_path( myConfig.pathOutDir, (char *) "table_ion_e1h1_p3.dat" ), "r" );    
    fp_ehe1_p2 = fopen( utils_concat_path( myConfig.pathOutDir, (char *) "table_ion_ehe1_p2.dat" ), "r" );
    fp_ehe1_p3 = fopen( utils_concat_path( myConfig.pathOutDir, (char *) "table_ion_ehe1_p3.dat" ), "r" );    
    fp_ehe2_p3 = fopen( utils_concat_path( myConfig.pathOutDir, (char *) "table_ion_ehe2_p3.dat" ), "r" );

    fp_temp_e1h1_p1 = fopen( utils_concat_path( myConfig.pathOutDir, (char *) "table_temp_e1h1_p1.dat" ), "r" );
    fp_temp_e1h1_p2 = fopen( utils_concat_path( myConfig.pathOutDir, (char *) "table_temp_e1h1_p2.dat" ), "r" );
    fp_temp_e1h1_p3 = fopen( utils_concat_path( myConfig.pathOutDir, (char *) "table_temp_e1h1_p3.dat" ), "r" );
    fp_temp_ehe1_p2 = fopen( utils_concat_path( myConfig.pathOutDir, (char *) "table_temp_ehe1_p2.dat" ), "r" );
    fp_temp_ehe1_p3 = fopen( utils_concat_path( myConfig.pathOutDir, (char *) "table_temp_ehe1_p3.dat" ), "r" );
    fp_temp_ehe2_p3 = fopen( utils_concat_path( myConfig.pathOutDir, (char *) "table_temp_ehe2_p3.dat" ), "r" );

    fp_comp1_p1 = fopen( utils_concat_path( myConfig.pathOutDir, (char *) "table_comp1_p1.dat" ), "r" );
    fp_comp1_p2 = fopen( utils_concat_path( myConfig.pathOutDir, (char *) "table_comp1_p2.dat" ), "r" );
    fp_comp1_p3 = fopen( utils_concat_path( myConfig.pathOutDir, (char *) "table_comp1_p3.dat" ), "r" );
    fp_comp2_p1 = fopen( utils_concat_path( myConfig.pathOutDir, (char *) "table_comp2_p1.dat" ), "r" );
    fp_comp2_p2 = fopen( utils_concat_path( myConfig.pathOutDir, (char *) "table_comp2_p2.dat" ), "r" );
    fp_comp2_p3 = fopen( utils_concat_path( myConfig.pathOutDir, (char *) "table_comp2_p3.dat" ), "r" );
    

    e1h1_p1 = gsl_vector_alloc( array_len_x1a);
    e1h1_p2 = gsl_vector_alloc( (array_len_x1a*array_len_x2a)  );
    e1h1_p3 = gsl_vector_alloc( (array_len_x2a*array_len_x13a) );

    ehe1_p2 = gsl_vector_alloc( (array_len_x1a*array_len_x2a)  );
    ehe1_p3 = gsl_vector_alloc( (array_len_x2a*array_len_x13a) );

    ehe2_p3 = gsl_vector_alloc( (array_len_x2a*array_len_x13a) );
 
    temp_e1h1_p1=gsl_vector_alloc(array_len_x1a);
    temp_e1h1_p2=gsl_vector_alloc((array_len_x1a*array_len_x2a));
    temp_e1h1_p3=gsl_vector_alloc((array_len_x2a*array_len_x13a));
    
    temp_ehe1_p2=gsl_vector_alloc((array_len_x1a*array_len_x2a));
    temp_ehe1_p3=gsl_vector_alloc((array_len_x2a*array_len_x13a));
    
    temp_ehe2_p3=gsl_vector_alloc((array_len_x2a*array_len_x13a));

    comp1_p1=gsl_vector_alloc(array_len_x1a);
    comp1_p2=gsl_vector_alloc((array_len_x1a*array_len_x2a));
    comp1_p3=gsl_vector_alloc((array_len_x2a*array_len_x13a));
    
    comp2_p1=gsl_vector_alloc(array_len_x1a);
    comp2_p2=gsl_vector_alloc((array_len_x1a*array_len_x2a));
    comp2_p3=gsl_vector_alloc((array_len_x2a*array_len_x13a));

 /**** Only temporary ***/ 

    gsl_matrix *tmpe1h1_p1=gsl_matrix_alloc(array_len_x1a,2);
    gsl_matrix *tmpe1h1_p2=gsl_matrix_alloc((array_len_x1a*array_len_x2a),3);
    gsl_matrix *tmpe1h1_p3=gsl_matrix_alloc((array_len_x2a*array_len_x13a),3);
    
    gsl_matrix *tmpehe1_p2=gsl_matrix_alloc((array_len_x1a*array_len_x2a),3);
    gsl_matrix *tmpehe1_p3=gsl_matrix_alloc((array_len_x2a*array_len_x13a),3);
    
    gsl_matrix *tmpehe2_p3=gsl_matrix_alloc((array_len_x2a*array_len_x13a),3);
    
    gsl_matrix *tmptemp_e1h1_p1=gsl_matrix_alloc(array_len_x1a,2);
    gsl_matrix *tmptemp_e1h1_p2=gsl_matrix_alloc((array_len_x1a*array_len_x2a),3);
    gsl_matrix *tmptemp_e1h1_p3=gsl_matrix_alloc((array_len_x2a*array_len_x13a),3);
    
    gsl_matrix *tmptemp_ehe1_p2=gsl_matrix_alloc((array_len_x1a*array_len_x2a),3);
    gsl_matrix *tmptemp_ehe1_p3=gsl_matrix_alloc((array_len_x2a*array_len_x13a),3);
    
    gsl_matrix *tmptemp_ehe2_p3=gsl_matrix_alloc((array_len_x2a*array_len_x13a),3);
    
    gsl_matrix *tmpcomp1_p1=gsl_matrix_alloc(array_len_x1a,2);
    gsl_matrix *tmpcomp1_p2=gsl_matrix_alloc((array_len_x1a*array_len_x2a),3);
    gsl_matrix *tmpcomp1_p3=gsl_matrix_alloc((array_len_x2a*array_len_x13a),3);
    
    gsl_matrix *tmpcomp2_p1=gsl_matrix_alloc(array_len_x1a,2);
    gsl_matrix *tmpcomp2_p2=gsl_matrix_alloc((array_len_x1a*array_len_x2a),3);
    gsl_matrix *tmpcomp2_p3=gsl_matrix_alloc((array_len_x2a*array_len_x13a),3);

 /* end of temporary*/ 

    gsl_matrix_fscanf(fp_e1h1_p1, tmpe1h1_p1);
    gsl_matrix_get_col(e1h1_p1, tmpe1h1_p1, 1);

    gsl_matrix_fscanf(fp_e1h1_p2,tmpe1h1_p2);  
    gsl_matrix_get_col(e1h1_p2,tmpe1h1_p2,2); 

    gsl_matrix_fscanf(fp_e1h1_p3,tmpe1h1_p3);
    gsl_matrix_get_col(e1h1_p3,tmpe1h1_p3,2); 
    
    gsl_matrix_fscanf(fp_ehe1_p2,tmpehe1_p2); 
    gsl_matrix_get_col(ehe1_p2,tmpehe1_p2,2); 

    gsl_matrix_fscanf(fp_ehe1_p3,tmpehe1_p3);
    gsl_matrix_get_col(ehe1_p3,tmpehe1_p3,2);

    gsl_matrix_fscanf(fp_ehe2_p3,tmpehe2_p3);
    gsl_matrix_get_col(ehe2_p3,tmpehe2_p3,2);

    if(DEBUG)log_debug(" Ionization integrals files read successfully");
 
/* ----------- For the second set of integrals -------------------*/


    gsl_matrix_fscanf(fp_temp_e1h1_p1,tmptemp_e1h1_p1);
    gsl_matrix_get_col(temp_e1h1_p1,tmptemp_e1h1_p1,1);

    gsl_matrix_fscanf(fp_temp_e1h1_p2,tmptemp_e1h1_p2);
    gsl_matrix_get_col(temp_e1h1_p2,tmptemp_e1h1_p2,2);

    gsl_matrix_fscanf(fp_temp_e1h1_p3,tmptemp_e1h1_p3);
    gsl_matrix_get_col(temp_e1h1_p3,tmptemp_e1h1_p3,2);

    gsl_matrix_fscanf(fp_temp_ehe1_p2,tmptemp_ehe1_p2);
    gsl_matrix_get_col(temp_ehe1_p2,tmptemp_ehe1_p2,2);

    gsl_matrix_fscanf(fp_temp_ehe1_p3,tmptemp_ehe1_p3);
    gsl_matrix_get_col(temp_ehe1_p3,tmptemp_ehe1_p3,2);

    gsl_matrix_fscanf(fp_temp_ehe2_p3,tmptemp_ehe2_p3);
    gsl_matrix_get_col(temp_ehe2_p3,tmptemp_ehe2_p3,2);

    if(DEBUG)log_debug(" Temperature integrals files read successfully");

    /*-------------For the third set of integrals ---------------------*/

    gsl_matrix_fscanf(fp_comp1_p1,tmpcomp1_p1);
    gsl_matrix_get_col(comp1_p1,tmpcomp1_p1,1);

    gsl_matrix_fscanf(fp_comp1_p2,tmpcomp1_p2); 
    gsl_matrix_get_col(comp1_p2,tmpcomp1_p2,2);

    gsl_matrix_fscanf(fp_comp1_p3,tmpcomp1_p3);
    gsl_matrix_get_col(comp1_p3,tmpcomp1_p3,2);

    gsl_matrix_fscanf(fp_comp2_p1,tmpcomp2_p1);
    gsl_matrix_get_col(comp2_p1,tmpcomp2_p1,1);

    gsl_matrix_fscanf(fp_comp2_p2,tmpcomp2_p2);
    gsl_matrix_get_col(comp2_p2,tmpcomp2_p2,2);

    gsl_matrix_fscanf(fp_comp2_p3,tmpcomp2_p3);
    gsl_matrix_get_col(comp2_p3,tmpcomp2_p3,2);

    if(DEBUG)log_debug(" Compton integrals files read successfully");

    fclose(fp_e1h1_p1);
    fclose(fp_e1h1_p2);
    fclose(fp_e1h1_p3);
    
    fclose(fp_ehe1_p2);
    fclose(fp_ehe1_p3);

    fclose(fp_ehe2_p3);


    fclose(fp_temp_e1h1_p1);
    fclose(fp_temp_e1h1_p2);
    fclose(fp_temp_e1h1_p3);

    fclose(fp_temp_ehe1_p2);
    fclose(fp_temp_ehe1_p3);

    fclose(fp_temp_ehe2_p3);


    fclose(fp_comp1_p1);
    fclose(fp_comp1_p2);
    fclose(fp_comp1_p3);

    fclose(fp_comp2_p1);
    fclose(fp_comp2_p2);
    fclose(fp_comp2_p3);

    if(DEBUG)log_debug(" Setting up matrices");

    matrix_e1h1_p2 = gsl_matrix_alloc(array_len_x1a,array_len_x2a);
    matrix_e1h1_p3 = gsl_matrix_alloc(array_len_x13a,array_len_x2a); 
    
    matrix_ehe1_p2 = gsl_matrix_alloc(array_len_x1a,array_len_x2a);
    matrix_ehe1_p3 = gsl_matrix_alloc(array_len_x13a,array_len_x2a);

    matrix_ehe2_p3 = gsl_matrix_alloc(array_len_x13a,array_len_x2a);
    

    temp_matrix_e1h1_p2 = gsl_matrix_alloc(array_len_x1a,array_len_x2a);
    temp_matrix_e1h1_p3 = gsl_matrix_alloc(array_len_x13a,array_len_x2a); 
    
    temp_matrix_ehe1_p2 = gsl_matrix_alloc(array_len_x1a,array_len_x2a);
    temp_matrix_ehe1_p3 = gsl_matrix_alloc(array_len_x13a,array_len_x2a);

    temp_matrix_ehe2_p3 = gsl_matrix_alloc(array_len_x13a,array_len_x2a);


    comp1_matrix_p2 = gsl_matrix_alloc(array_len_x1a,array_len_x2a);
    comp1_matrix_p3 = gsl_matrix_alloc(array_len_x13a,array_len_x2a);

    comp2_matrix_p2 = gsl_matrix_alloc(array_len_x1a,array_len_x2a);
    comp2_matrix_p3 = gsl_matrix_alloc(array_len_x13a,array_len_x2a); 

    if(DEBUG)log_debug(" Filling ionization matrices");
    
    for(i=0;i<array_len_x1a;i++)
        for(j=0;j<array_len_x2a;j++)
            gsl_matrix_set(matrix_e1h1_p2,i,j,gsl_vector_get(e1h1_p2,i*array_len_x2a + j));

    for(i=0;i<array_len_x13a;i++)
        for(j=0;j<array_len_x2a;j++)
            gsl_matrix_set(matrix_e1h1_p3,i,j,gsl_vector_get(e1h1_p3,i*array_len_x2a + j));

    
    for(i=0;i<array_len_x1a;i++)
        for(j=0;j<array_len_x2a;j++)
            gsl_matrix_set(matrix_ehe1_p2,i,j,gsl_vector_get(ehe1_p2,i*array_len_x2a + j));

    for(i=0;i<array_len_x13a;i++)
        for(j=0;j<array_len_x2a;j++)
            gsl_matrix_set(matrix_ehe1_p3,i,j,gsl_vector_get(ehe1_p3,i*array_len_x2a + j));

        
    for(i=0;i<array_len_x13a;i++)
        for(j=0;j<array_len_x2a;j++)
            gsl_matrix_set(matrix_ehe2_p3,i,j,gsl_vector_get(ehe2_p3,i*array_len_x2a + j));


    if(DEBUG)log_debug(" Filling temperature matrices");
 /*---------for the second set of integrals ---------------------*/


    for(i=0;i<array_len_x1a;i++)
        for(j=0;j<array_len_x2a;j++)
            gsl_matrix_set(temp_matrix_e1h1_p2,i,j,gsl_vector_get(temp_e1h1_p2,i*array_len_x2a + j));

    for(i=0;i<array_len_x13a;i++)
        for(j=0;j<array_len_x2a;j++)
            gsl_matrix_set(temp_matrix_e1h1_p3,i,j,gsl_vector_get(temp_e1h1_p3,i*array_len_x2a + j));


    for(i=0;i<array_len_x1a;i++)
        for(j=0;j<array_len_x2a;j++)
            gsl_matrix_set(temp_matrix_ehe1_p2,i,j,gsl_vector_get(temp_ehe1_p2,i*array_len_x2a + j));

    for(i=0;i<array_len_x13a;i++)
        for(j=0;j<array_len_x2a;j++)
            gsl_matrix_set(temp_matrix_ehe1_p3,i,j,gsl_vector_get(temp_ehe1_p3,i*array_len_x2a + j));



    for(i=0;i<array_len_x13a;i++)
        for(j=0;j<array_len_x2a;j++)
            gsl_matrix_set(temp_matrix_ehe2_p3,i,j,gsl_vector_get(temp_ehe2_p3,i*array_len_x2a + j));



 /*----------for the third set of integrals -----------------------*/
    if(DEBUG)log_debug(" Filling Compton matrices");
    
    for(i=0;i<array_len_x1a;i++)
        for(j=0;j<array_len_x2a;j++)
            gsl_matrix_set(comp1_matrix_p2,i,j,gsl_vector_get(comp1_p2,i*array_len_x2a + j));

    for(i=0;i<array_len_x13a;i++)
        for(j=0;j<array_len_x2a;j++)
            gsl_matrix_set(comp1_matrix_p3,i,j,gsl_vector_get(comp1_p3,i*array_len_x2a + j));

    
    for(i=0;i<array_len_x1a;i++)
        for(j=0;j<array_len_x2a;j++)
            gsl_matrix_set(comp2_matrix_p2,i,j,gsl_vector_get(comp2_p2,i*array_len_x2a + j));

    for(i=0;i<array_len_x13a;i++)
        for(j=0;j<array_len_x2a;j++)
            gsl_matrix_set(comp2_matrix_p3,i,j,gsl_vector_get(comp2_p3,i*array_len_x2a + j));


		    /************freeing the matrix allocations*******/
    if(DEBUG)log_debug(" Freeing memory");
    
    gsl_matrix_free(tmpe1h1_p1);
    gsl_matrix_free(tmpe1h1_p2);
    gsl_matrix_free(tmpe1h1_p3);

    gsl_matrix_free(tmpehe1_p2);
    gsl_matrix_free(tmpehe1_p3);

    gsl_matrix_free(tmpehe2_p3);

    gsl_matrix_free(tmptemp_e1h1_p1);
    gsl_matrix_free(tmptemp_e1h1_p2);
    gsl_matrix_free(tmptemp_e1h1_p3);

    gsl_matrix_free(tmptemp_ehe1_p2);
    gsl_matrix_free(tmptemp_ehe1_p3);

    gsl_matrix_free(tmptemp_ehe2_p3);

    gsl_matrix_free(tmpcomp1_p1);
    gsl_matrix_free(tmpcomp1_p2);
    gsl_matrix_free(tmpcomp1_p3);

    gsl_matrix_free(tmpcomp2_p1);
    gsl_matrix_free(tmpcomp2_p2);
    gsl_matrix_free(tmpcomp2_p3);

    printf(" Pre-computed integrals successfully loaded into memory\n");
}
