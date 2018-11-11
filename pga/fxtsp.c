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
#include "types.h"
#include "extern.h"
#include "fxn.h"
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

   // Alphabet size depends on size of tensor
   Alphabet_size = tsp.num_cities;
 
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
   int i;
   double distance;
   FILE *fp;
   double floats[tsp.num_cities];
   KEY_MAP maps[tsp.num_cities];
   int ordering[tsp.num_cities]; 
   
   // get ordering from proportional ratios
   
   for (i=0; i < tsp.num_cities; i++) {
	   floats[i] = get_value(i, indv);
       
       // floats[i] = get_value(j, indv);
       
   }

   for (i = 0; i < tsp.num_cities; i++) {
      maps[i].key = floats[i];
      maps[i].value = i;
   }
   
   // printf("order: ");
   qsort(maps, tsp.num_cities, sizeof(KEY_MAP), compare);
   for (i = 0; i < tsp.num_cities; i++) {
       ordering[i] = maps[i].value;
       // printf("%d ", ordering[i]);
   }
   // printf("\n");
   
  /* total distance traveled */
   distance = 0;
   
   // printf("distances: ");
   
   // euclidean distance
   for (i = 0; i < tsp.num_cities; i++)
      { 
      COORDS *a = tsp.city_coordinates[ordering[i]];
      
      COORDS *b;
      if (i != tsp.num_cities-1) { 
         b = tsp.city_coordinates[ordering[i+1]];
      } else {
         // last item in chromosome is connected to first item
         b = tsp.city_coordinates[ordering[0]];
      }
           
      distance += hypot(a->x - b->x, a->y - b->y);
      // printf("%.1lf ", hypot(a->x - b->x, a->y - b->y));
      }
    
   indv->fitness = -1 * distance;
  
  // printf("\n");
  // printf(" --- Press any key to continue ---\n");  fgetc(stdin);
  
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
   if (Gen.best_fitness >= -423.741)
      return 1;
   return 0;
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

/********** fxn_gen_output **********/
/* parameters:
   called by:	gen_output(), output.c
   actions:	anything that needs to be printed at the end of
		each generation that is function specific.
*/
void fxn_gen_output()
   {
   // int i;
   // double value;

// #ifdef DEBUG
   // printf(" ---in fxn_gen_output---\n");
// #endif
 
   // if (match.opfile_values.on)
      // {
      // match.opfile_values.fp = fopen(match.opfile_values.filename, "a");
      // fprintf(match.opfile_values.fp, " %4d best:xter_prop", Gen.index);
      // for (i=0; i<match.target[match.current].num_val; i++)
         // {
// /*****         value = get_value(i, Pop[Gen.best_indv_index]); *****/
         // value = get_value(i, Pop[Gen.best_indv_index]);
        // /* print evolved value from indv followed by desired value */
         // fprintf(match.opfile_values.fp, " %s%s %lf %lf",
		// match.target[match.current].val[i].positive_xter,
		// match.target[match.current].val[i].negative_xter,
		// value, match.target[match.current].val[i].val);
         // }  /* for */
      // fprintf(match.opfile_values.fp, " avg_pop:xter_prop");
     // /* NOTE(010513AW): for the following to print correctly, each pair of
        // related letters must be specified next to each other */
      // for (i=0; i<Alphabet_size; i++)
         // {
         // fprintf(match.opfile_values.fp, " ");
         // putc(Xters[i], match.opfile_values.fp);
         // fprintf(match.opfile_values.fp, " %lf", Gen.avg_xter_proportion[i]);
// /*
                 // Gen.avg_xter_proportion[i]+Gen.avg_xter_proportion[i+1]);
// */
         // }  /* for i */
      // fprintf(match.opfile_values.fp, "\n");
      // fclose(match.opfile_values.fp);
      // }  /* if print values */

// #ifdef DEBUG
   // printf(" ---end fxn_gen_output---\n");
// #endif

   }  /* fxn_gen_output */

/********** fxn_gen_start **********/
/* parameters:  
   called by:   gen_start(), stats.c 
   actions:     anything that needs to be done at the start of
                each generation that is function specific.
*/
int fxn_gen_start()
   {
   return OK;
   }  /* fxn_gen_start */

/********** fxn_gen_end **********/
/* parameters:  
   called by:   gen_end(), stats.c 
   actions:     anything that needs to be done at the end of
                each generation that is function specific.
*/
int fxn_gen_end()
   {
// #ifdef DEBUG
   // printf(" ---in fxn_gen_end---\n");
// #endif

// #ifdef DEBUG_MULT
      // printf("          ...... gen %d gencount %d current %d\n",
		// Gen.index, match.gencount, match.current);
// #endif
   // match.gencount++;
   // if (match.gencount >= match.target[match.current].duration)
      // {
// #ifdef DEBUG_MULT
      // printf(" >>>>>>>>>>>>>>> gen %d gencount %d current %d\n",
		// Gen.index, match.gencount, match.current);
// #endif
      // match.current = (match.current + 1) % match.num_targets;
      // match.gencount = 0;

      // if (Print_targets)
         // printf(" Target #%d:  %s, duration %d, match.gencount %d\n",
		// match.current, match.target[match.current].name,
		// match.target[match.current].duration,
		// match.gencount);
      // }  /* if */

// #ifdef DEBUG
   // printf(" ---end fxn_gen_end---\n");
// #endif

   return OK;
   }  /* fxn_gen_end */
   
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
   
/********** get_value *********/
/* parameters:	indv
   called by: 	gen_output(), output.c
		fxn_gen_output()
   action:	calculate a value represented by the selected indv.
*/
double get_value(int which_val, INDIVIDUAL *indv)
   {
   int i;
   double pos_count = 0;
   double neg_count = 0;
   double sum_parts = 0;

  /* calculate value represented by individual */
   for (i=0; i<indv->length; i++)
      {
	  if (indv->genome[i] == which_val) 
	     {
		 pos_count++; 
		 }
	  else if (indv->genome[i] == -1 * which_val) 
	     {
		 neg_count++; 
		 }
      }
 
   sum_parts = pos_count + neg_count;

   if (sum_parts == 0)
	  return 0.0;
   else if (sum_parts > 0)
      return (pos_count/sum_parts);
   }  /* get_value */
   
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
      sscanf(aline, "%d %d", &tsp.city_coordinates[i]->x, &tsp.city_coordinates[i]->y);
      i++;
      }  /* while */
   return OK;
   }  /* read_city_coordiantes */
   