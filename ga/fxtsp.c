/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* fxn.c
   07.27.07.RN	Created.

		Traveling Salesman Problem

   Routines:    read_fxn_file()         07.23.98.AW
   (required)   init_function()         07.09.98.AW
                end_function()          09.16.98.AW
                eval_indv()             07.19.17 RN
                found_solution()        07.19.17 RN
                fprint_fxn()            09.16.98.AW
                fprint_genes()          09.23.98.AW
		fxn_fprint_gen_indv()	10.07.98.AW

   Routines:	(dummy added because RR functions have a call and the
		GA needs to be able to compile a call to this function
		line 164 in ga.c)
		trace_bb_data()		12.30.98.AW

   Each new function is expected to provide it's own set of routines.
   Need to provide one of the following:
		- routine for reading in the function file
		- routine for initializing the function
                - routine for finalizing the function
		- routine for evaluating one individual
		- routine to decide on stopping condition
                - routine that prints out all info about function
                - routine to print genes or bb to file (mainly for Trace)
		- routine to print one indv plus generation statistics (for
		  Print_fxn_best - bc floating bbs require two lines to display)
   (07.23.98.AW I will keep adding things to this list as needed)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "types.h"
#include "extern.h"
#include "google_call.h"
#include "fxtsp.h"

#include<Python.h>

/* file scope object */
static double **dist_matrix;    // matrix of distance values from previous API calls
static double **time_matrix;    // matrix of time values from previous API calls
static char *d_matrix_file;     // file for reading/writing distance matrix information
static char *t_matrix_file;     // file for reading/writing time matrix information
static int matrix_changed = 0;  // flag checked before writing dist_matrix to file
static double* dist_time;       // array to hold one distance and time from google call

/* internal routine prototypes */
int allocate_coords_space();
int compare(const void *left, const void *right);
int init_matrix();
int free_matrix();
int write_matrix();
int make_matrix_filename(char *fxn_file);


/********** read_fxn_file **********/
/* parameters:
   called by:   ga_init(), ga.c
   actions:     Initialize fitness function if needed.
*/
int read_fxn_file(char *fxn_file)
   {
   FILE *fp;
   char *aline;
   char name[INPUT_LINE_LEN];

#ifdef DEBUG
   printf(" ---in read_fxn_file---\n");
#endif

   fp = fopen(fxn_file, "r");
   if (fp == NULL)
      {
      printf(" Error(read_fxn_file): cannot open file: %s\n", fxn_file);
      return ERROR;
      }
   printf(" Reading from file: %s\n", fxn_file);

  /* allocate space for aline */
   aline = (char *)malloc(INPUT_LINE_LEN * sizeof(char));
 
  /* read function name */
   if (get_next_line(fp, aline) != ENDOFFILE)
      sscanf(aline, "%s %s", name, Function_name);

  /* read rest of function info */
   if (get_next_line(fp, aline) != ENDOFFILE)
      sscanf(aline, "%s %d", name, &tsp.num_cities);
   if (allocate_coords_space() == ERROR)  return ERROR;
   if (read_city_coordiantes(fp, aline) == ERROR)  return ERROR;
  
  /* call routine to read rest of function info */
   printf("       Nothing else to read for Traveling Salesman Problem.\n");
 
   free(aline);
   fclose(fp);
   
   Max_gen_len = tsp.num_cities;

#ifdef DEBUG
   printf(" ---end read_fxn_file---\n");
#endif
   return OK;
   }  /* read_fxn_file */


/****** make_matrix_filename *******/
/* parameters: char *fxn_file - input filename
   called by: init_function, fxtsp.c
   actions:   create filename for dist_matrix based on input filename
*/
int make_matrix_filename(char *fxn_file)
{
    int len = strlen(fxn_file);
    
    d_matrix_file = malloc((len - 5) * sizeof(char));
    t_matrix_file = malloc((len - 5) * sizeof(char));
    
    int i = 0;
    while(fxn_file[i+3] != '.')
    {
        d_matrix_file[i] = fxn_file[i+3];
        t_matrix_file[i] = fxn_file[i+3];
        i++;
    }
    d_matrix_file[i] = '\0';
    strcat(d_matrix_file, "_d.txt");
    t_matrix_file[i] = '\0';
    strcat(t_matrix_file, "_t.txt");
    
}


/********** init_matrix ************/
/* parameters:
   called by:   init_function(), fxtsp.c
   actions:     allocate space for and initialize dist_matrix and time_matrix
*/
int init_matrix()
{
   FILE *fp_d;
   FILE *fp_t;
   char *mline;
    
   // allocate space for the matix
   // it is num_cities + 1) x (num_cities + 1) due to inclusion of origin
   dist_matrix = malloc((tsp.num_cities + 1) * sizeof(double*));
   time_matrix = malloc((tsp.num_cities + 1) * sizeof(double*));
   for(int i = 0; i <= tsp.num_cities; i++)
   {
       dist_matrix[i] = malloc((tsp.num_cities + 1) * sizeof(double));
       time_matrix[i] = malloc((tsp.num_cities + 1) * sizeof(double));
   }
   printf("      dist_matrix and time_matrix allocated\n");

   // if the matrix_file exists, read distance data from it
   // if it doesn't exist, initialize the matrix to -1.0
   fp_d = fopen(d_matrix_file, "r");
   fp_t = fopen(t_matrix_file, "r");
   mline = malloc(INPUT_LINE_LEN * sizeof(char));
    
   if (fp_d != NULL)
   {
        for(int i = 0; i <= tsp.num_cities; i++)
        {
            for(int j = 0; j <= tsp.num_cities; j++)
            {
              // process an entry for the distance matrix
              if (get_next_line(fp_d, mline) == ENDOFFILE)
              {
                 printf(" Error(read_matrix_file): unexpected end of file - distance\n");
                 return ERROR;
              }  /* if */
              
              sscanf(mline, "%lf", &dist_matrix[i][j]);

            }

         }  /* for */
        
        fclose(fp_d);
        printf("      dist_matrix read from file\n");
   }
   else
   {
       // initialize all elements to -1.0
       for(int i = 0; i <= tsp.num_cities; i++)
       {
           for(int j = 0; j <= tsp.num_cities; j++)
           {
               dist_matrix[i][j] = -1.0;
               time_matrix[i][j] = -1.0;
           }
       }
       printf("      dist_matrix initialized to -1.0\n");
   }

   if (fp_t != NULL)
   {
        for(int i = 0; i <= tsp.num_cities; i++)
        {
            for(int j = 0; j <= tsp.num_cities; j++)
            {
              // process an entry for the distance matrix
              if (get_next_line(fp_t, mline) == ENDOFFILE)
              {
                 printf(" Error(read_matrix_file): unexpected end of file - time\n");
                 return ERROR;
              }  /* if */
              
              sscanf(mline, "%lf", &time_matrix[i][j]);

            }

         }  /* for */
        
        fclose(fp_t);
        printf("      time_matrix read from file\n");
   }
   else
   {
       // initialize all elements to -1.0
       for(int i = 0; i <= tsp.num_cities; i++)
       {
           for(int j = 0; j <= tsp.num_cities; j++)
           {
               time_matrix[i][j] = -1.0;
           }
       }
       printf("      time_matrix initialized to -1.0\n");
   }
                
   free(mline);
   return OK;
}


/********** free_matrix ************/
/* parameters:
   called by: end_function(), fxtsp.c
   actions: deallocate the dist_matrix
*/
int free_matrix()
{
    for(int i = 0; i < tsp.num_cities; i++)
    {
        free(dist_matrix[i]);
        free(time_matrix[i]);
    }
    free(dist_matrix);
    free(time_matrix);
    
    return OK;
}


/********** write_matrix ***********/
/* parameters:
   called by: end_function(), fxtsp.c
   actions: Write the dist_matrix and time_matrix to files
            for use in later runs
*/
int write_matrix()
{
    FILE *fp_d;
    FILE *fp_t;
    
   fp_d = fopen(d_matrix_file, "w");
   fp_t = fopen(t_matrix_file, "w");
   if(fp_d == NULL)
   {
      printf(" Error(write_matrix_file): cannot open file: %s\n", d_matrix_file);
      return ERROR;              
   }
   else if(fp_t == NULL)
   {
       printf(" Error(write_matrix_file): cannot open file: %s\n", t_matrix_file);
      return ERROR;       
   }
   else
   {
       for(int i = 0; i <= tsp.num_cities; i++)
       {
           for(int j = 0; j <= tsp.num_cities; j++)
           {
               fprintf(fp_d, "%f\n", dist_matrix[i][j]);
               fprintf(fp_t, "%f\n", time_matrix[i][j]);
           }
       }
   fclose(fp_d);
   fclose(fp_t);
   }
    
   return OK;
}

/********** init_function **********/
/* parameters:
   called by:	ga_init(), ga.c
   actions:	Initialize fitness function if needed.
*/
int init_function(char *fxn_file)
   {
#ifdef DEBUG
   printf(" ---in init_function---\n");
#endif
    
    printf("Initializing function: %s\n", Function_name);

    // required by the C Python library
    Py_Initialize();
    printf("   completed Py_Initialize\n\n");
    
    // make matrix filename from fxn_file
    make_matrix_filename(fxn_file);
    
    // create and initialize the distance matrix
    init_matrix();
    printf("   created dist_matrix and time_matrix\n");

    // allocate space for dist_time
    dist_time = malloc(2 * sizeof(double));
    
    printf("\n");
    
#ifdef DEBUG
   printf(" ---end init_function---\n");
#endif
   return OK;
   }  /* init_function */

/********** end_function **********/
/* parameters:
   called by:   ga_end(), ga.c
   actions:     Finalize fitness function if needed.
*/
void end_function()
   {
#ifdef DEBUG
   printf(" ---in end_function---\n");
#endif
 
   printf("Finalizing function: %s\n", Function_name);

   // cleanup memory used by the Python module
   Py_Finalize();
   printf("   completed Py_Finalize\n");
    
   // write the dist_matrix to a file
   if(matrix_changed)
   {
       write_matrix();
       printf("   wrote dist_matrix and time_matrix to file\n");
   }
   else
   {
       printf("   dist_matrix and time_matrix unchanged\n");
   }
    
   // deallocate the dist_matrix
   free_matrix();
   printf("   deallocated dist_matrix and time_matrix\n\n");
    
   // deallocate dist_time
   free(dist_time);
    
#ifdef DEBUG
   printf(" ---end end_function---\n");
#endif
   }  /* end_function */

/********** eval_indv **********/
/* parameters:	indv		to evaluate
   called by:	pop_eval(), pop.c
   actions:	given an individual, evaluates the fitness of the
		genome and assigns the appropriate value to the
		fitness field.
*/
void eval_indv(INDIVIDUAL *indv)
   {
   int i, src, dest;
   double distance, time;
   FILE *fp;
   
   // static int init = 0;             // flag used to avoid reinitializing the matrix
   static int google_count = 0;     // count of calls to the API
   // static double **dist_matrix;     // matrix of values from previous API calls
    
   /* decode if using random keys */
   if (Init_pop == 2)
      {
      decode(indv);
      }
   
   /* total distance traveled and time required */
   distance = 0;
   time = 0;

   // calculate distance using Google API
   // the loop calculates distances for lpcations in the genome
   // distance to/from the origin are calcultaed after the loop
   for (i = 0; i < indv->length - 1; i++)
   { 
      src = indv->genome[i];
      dest = indv->genome[i+1];
      COORDS *a = tsp.city_coordinates[src];
      COORDS *b = tsp.city_coordinates[dest];

      // Before making a call to the API, check to see if there is
      // a valid value in the matrix (matrix is intitialized to -1.0)
      if(dist_matrix[src][dest] < 0.0)
      {
          double s[] = {a->lat, a->lon};
          double d[] = {b->lat, b->lon};
          printf("google_count: %d\n", google_count);
          printf("%d -> %d google call: (%f, %f)  (%f, %f)\n", src, dest, a->lat, a->lon, b->lat, b->lon);
          google_dist(s, d, dist_time);
          printf("%d -> %d distance: %f\n", src, dest, dist_time[0]);
          printf("%d -> %d time: %f\n\n", src, dest, dist_time[1]);
          dist_matrix[src][dest] = dist_time[0];
          time_matrix[src][dest] = dist_time[1];
          google_count++;
          matrix_changed = 1;
      }
//      else
//      {
//           segment_d = dist_matrix[src][dest];
//           segment_t = time_matrix[src][dest];
//      }
       
   distance += dist_matrix[src][dest];
   time += time_matrix[src][dest];
   }
    
   // Get distance from origin to first location
   src = tsp.num_cities;
   dest = indv->genome[0];
   if(dist_matrix[src][dest] < 0.0)
   {
       COORDS *a = tsp.origin;
       COORDS *b = tsp.city_coordinates[dest];
       double s[] = {a->lat, a->lon};
       double d[] = {b->lat, b->lon};
       printf("google_count: %d\n", google_count);
       printf("origin -> %d google call: (%f, %f)  (%f, %f)\n", dest, a->lat, a->lon, b->lat, b->lon);
       google_dist(s, d, dist_time);
       printf("%d -> %d distance: %f\n", src, dest, dist_time[0]);
       printf("%d -> %d time: %f\n\n", src, dest, dist_time[1]);
       dist_matrix[src][dest] = dist_time[0];
       time_matrix[src][dest] = dist_time[1];
       google_count++;
       matrix_changed = 1;
   }
//   else
//   {
//       segment_d = dist_matrix[src][dest];
//       segment_t = time_matrix[src][dest];
//   }
   distance += dist_matrix[src][dest];
   time += time_matrix[src][dest];
    
   // get distance from last location back to origin
   src = indv->genome[indv->length-1];
   dest = tsp.num_cities;
   if(dist_matrix[src][dest] < 0.0)
   {
       COORDS *a = tsp.city_coordinates[src];
       COORDS *b = tsp.origin;
       double s[] = {a->lat, a->lon};
       double d[] = {b->lat, b->lon};
       printf("google_count: %d\n", google_count);
       printf("%d -> origin google call: (%f, %f)  (%f, %f)\n", src, a->lat, a->lon, b->lat, b->lon);
       google_dist(s, d, dist_time);
       printf("%d -> %d distance: %f\n", src, dest, dist_time[0]);
       printf("%d -> %d time: %f\n\n", src, dest, dist_time[1]);
       dist_matrix[src][dest] = dist_time[0];
       time_matrix[src][dest] = dist_time[1];
       google_count++;
       matrix_changed = 1;
   }
//   else
//   {
//       segment_d = dist_matrix[src][dest];
//       segment_t = time_matrix[src][dest];
//   }
   distance += dist_matrix[src][dest];
   time += time_matrix[src][dest];
// printf("\n");
   
   // fitness metric: 0 = distance, 1 = time
   if(Metric == 0)
       indv->fitness = distance;
   else
       indv->fitness = time;
}  /* eval_indv */


/********** found_solution **********/
/* parameters:
   called by:	ga_continue(), ga.c
   actions:	returns 1 if a good enough solution has been found,
		otherwise returns 0.
   07.19.17 RN: Converted to "integermax", solution is all 9s
*/
int found_solution()
   {
   if (Gen.best_fitness <= 2579.0)
      return 1;
   else  return 0;
   }  /* found_solution */

/********** fprint_fxn **********/
/* parameters:  fp      where to print
   called by:   which ever routine is being tested
   actions:     prints out all information about the function.
*/
void fprint_fxn(FILE *fp)
   {
   int i, j, k;

#ifdef DEBUG
   printf(" ---in fprint_fxn---\n");
#endif
 
   fprintf(fp, " Function: %s\n", Function_name);

#ifdef DEBUG
   printf(" ---end fprint_fxn---\n");
#endif
   }  /* fprint_fxn */

/********** fprint_genes **********/
/* parameters:
   called by:   eval_indv()
   actions:     print gene data to file
		mainly used for printing Trace data
*/
void fprint_genes(FILE *fp, INDIVIDUAL *indv)
   {
#ifdef DEBUG
   printf(" ---in fprint_genes---\n");
#endif
  /* tsp doesn't have genes */
   fprintf(fp, "0\n");
#ifdef DEBUG
   printf(" ---end fprint_genes---\n");
#endif
   }  /* fprint_genes */

/********** fxn_fprint_gen_indv **********/
/* parameters:
   called by:   gen_output(), output.C
   actions:	print generation statistics and one individual.
		print individual in function specific format -- for now, this
		mainly means that if the fxn has floating building blocks,
		then a second line marking to start of each building block
		is printed below each genome.
		For functions like onemax, this is the same as the standard
		printing format.
*/
void fxn_fprint_gen_indv(FILE *fp, int indv)
   {
#ifdef DEBUG
   printf(" ---in fxn_fprint_gen_indv---\n");
#endif

    fprintf(fp, " Gen %3d %8.3lf %8.3lf %8.3lf ", Gen.index, Gen.avg_fitness,
                 Gen.std_dev, Gen.best_fitness);
    fprint_genome(fp, Pop[indv], 1);
    
    // printf("%d\n", Pop[indv]->length);

#ifdef DEBUG
   printf(" ---end fxn_fprint_gen_indv---\n");
#endif
   }  /* fxn_fprint_gen_indv */

/********** function specific rqrd functions -- will be called by GA *********/
 
/********** trace_bb_data **********/
/* parameters:
   called by:   ga_loop(), ga.c
   actions:     traces construction and disruption of bb
   note:        expect Pop to contain the evaluated population (Gen.index)
                and Kids to contain the evaluated parents (Gen.index-1).

   12.30.98.AW	added because the GA would not compile fxn.c without
		this function call.  Calls it in ga.c ~ line 164.
		Was supposed to only be needed by RR functions.
*/
void trace_bb_data()
   {
   }  /* trace_bb_data */

void decode(INDIVIDUAL *indv) {
#ifdef DEBUG
   printf(" ---in decode---\n");
#endif
    
   int i;

   KEY_MAP maps[indv->length];

   for (i = 0; i < indv->length; i++) 
      {
      maps[i].key = indv->floats_genome[i];
      maps[i].value = i;
      }
    
   qsort(maps, indv->length, sizeof(KEY_MAP), compare);
 
   for (i = 0; i < indv->length; i++) {
       indv->genome[i] = maps[i].value;
   }
    
#ifdef DEBUG
   printf(" ---end decode---\n");
#endif
} /* decode */

/********** local functions -- not called from outside this file *********/

int compare(const void *left, const void *right)
   {
   const KEY_MAP a = *(KEY_MAP *)left; 
   const KEY_MAP b = *(KEY_MAP *)right;
   // printf("a: %.2f  b: %.2f  ", a, b);
   // printf("\n");
   if (a.key < b.key) return -1;
   if (a.key > b.key) return 1;
   return 0;

   } /* compare */

int allocate_coords_space()
   {
#ifdef DEBUG
   printf(" ---in allocate_coords_space---\n");
#endif
   
   tsp.origin = (COORDS *)malloc(sizeof(COORDS));
    
   if (tsp.origin == NULL)
      {
      printf(" Error(allocate_coords_space): cannot allocate: tsp.origin\n");
      return ERROR;
      }

   tsp.city_coordinates = (COORDS **)malloc(tsp.num_cities * sizeof(COORDS *));
   int i;
   for (i = 0; i < tsp.num_cities; i++)
      {
      tsp.city_coordinates[i] = (COORDS *)malloc(sizeof(COORDS)); 
      }
   
   if (tsp.city_coordinates == NULL)
      {
      printf(" Error(allocate_coords_space): cannot allocate: tsp.city_coordinate\n");
      return ERROR;
      }

#ifdef DEBUG
   printf(" ---end allocate_coords_space---\n");
#endif

   return OK;
   }  /* allocate_coords_space */
   
/********** read_basic_bb **********/
/* parameters:	fp, aline
   called by:	read_fxn_file()
   actions:	reads basic building blocks, skipping comment lines (#).
		Given the number of basic building blocks, no error
		checking is done to make sure this is correct.
*/
int read_city_coordiantes(FILE *fp, char *aline)
   {
   int i;

   // get the designated start/end location
   if (get_next_line(fp, aline) == ENDOFFILE)
   {
       printf(" Error(read_city_coordiantes): unexpected end of file\n");
       return ERROR;
   }  /* if */
   printf("%s", aline);
   sscanf(aline, "%lf %lf", &tsp.origin->lat, &tsp.origin->lon);

   // get the remaining locations
   i = 0;
   while (i < tsp.num_cities)
      {
      if (get_next_line(fp, aline) == ENDOFFILE)
         {
         printf(" Error(read_city_coordiantes): unexpected end of file\n");
         return ERROR;
         }  /* if */
         printf("%s", aline);
      sscanf(aline, "%lf %lf", &tsp.city_coordinates[i]->lat, &tsp.city_coordinates[i]->lon);
      i++;
      }  /* while */
   return OK;
   }  /* read_city_coordiantes */
   
