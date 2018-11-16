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

/* internal routine prototypes */
int allocate_coords_space();
int compare(const void *left, const void *right);

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

/********** init_function **********/
/* parameters:
   called by:	ga_init(), ga.c
   actions:	Initialize fitness function if needed.
*/
int init_function()
   {
#ifdef DEBUG
   printf(" ---in init_function---\n");
#endif

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
 
   printf(" Finalizing function: %s\n", Function_name);
 
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
   double distance, segment;
   FILE *fp;
   
   static int init = 0;             // flag used to avoid reinitializing the matrix
   static int google_count = 0;     // count of calls to the API
   static double **dist_matrix;     // matrix of values from previous API calls
    
   if(!init)
   {
       // allocate space for the matix
       dist_matrix = malloc(tsp.num_cities * sizeof(double*));
       for(int i = 0; i < tsp.num_cities; i++)
       {
           dist_matrix[i] = malloc(tsp.num_cities * sizeof(double));
       }

       // initialize all elements to -1.0
       for(int i = 0; i < tsp.num_cities; i++)
       {
           for(int j = 0; j < tsp.num_cities; j++)
               dist_matrix[i][j] = -1.0;
       }
       init = 1;
   }
    
   /* decode if using random keys */
   if (Init_pop == 2)
      {
      decode(indv);
      }
   
  /* total distance traveled */
   distance = 0;
      
   // distance using Google API
   for (i = 0; i < indv->length; i++)
   { 
      src = indv->genome[i];   
      COORDS *a = tsp.city_coordinates[indv->genome[i]];
      COORDS *b;
      if (i != indv->length-1) { 
         b = tsp.city_coordinates[indv->genome[i+1]];
         dest = indv->genome[i+1];
      } else {
         // last item in chromosome is connected to first item
         b = tsp.city_coordinates[indv->genome[0]];
         dest = indv->genome[0];
      }

      // usleep(200000);  // pause for 200ms
       
      // Before making a call to the API, check to see if there is
      // a valid value in the matrix (matrix is intitialized to -1.0)
      if(dist_matrix[src][dest] < 0.0)
      {
          double s[] = {a->lat, a->lon};
          double d[] = {b->lat, b->lon};
          segment = google_dist(s, d);
          dist_matrix[src][dest] = segment;
          google_count++;
      }
      else
      {
           segment = dist_matrix[src][dest];
      }
       
      distance += segment;
   }
    
   indv->fitness = distance;
}  /* eval_indv */


/********** init_dist_matrix **********/
/* parameters:	
   called by:	ga_init(), ga.c
   actions:	allocates space for the global distances matrix
            which is defined in globals.h.
   11.15.18 DM: Created function
*/
//double** init_dist_matrix()
//{
//    double **distances;
//    printf("In init distance matrix\n");
//    printf("   tsp.num_cities: %d\n", tsp.num_cities);
//    distances = malloc(tsp.num_cities * sizeof(double*));
//    printf("After first line\n");
//    for(int i = 0; i < tsp.num_cities; i++)
//    {
//        distances[i] = malloc(tsp.num_cities * sizeof(double));
//        printf("After an iteration of malloc loop\n");
//    }
//    
//    for(int i = 0; i < tsp.num_cities; i++)
//        for(int j = 0; j < tsp.num_cities; j++)
//            distances[i][j] = -1.0;
// 
//    
//printf("Leaving init distance matrix\n");        
//    return distances;
//}   /* init_dist_matrix */


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
 
   fprintf(fp, " Function: %s\n", Function_name);
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
   