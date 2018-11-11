/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* stats.c
   07.22.98.AW	Created.  Routines that maintain and calculate
		statistics during a run.

   Routines:	run_start()	07.22.98.AW
		run_end()	07.22.98.AW
		gen_start()	07.22.98.AW
		gen_end()	07.22.98.AW
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include "types.h"
#include "extern.h"
#include "stats.h"
#include "indv.h"
#include "output.h"
#include "trace.h"
#include "params.h"
#include "fxn.h"
#include "xtrack.h"

#define DEBUG_COUNT 1
#undef DEBUG_COUNT

/********** run_start **********/
/* parameters:
   called by:	ga_init(), ga.c
   actions:	Things to be done once at start of run.
		Initialize gen.index to zero and initialize
		the rest of the gen.* values so that calculations
		(via gen_end) may be done on the initial generation
		(generation 0).
		See also, gen_start().
*/
int run_start()
   {
   int ptr;
   int i;

#ifdef DEBUG
   printf(" ---in run_start---\n");
#endif

   Gen.index = 0;
   Gen.fitness_sum = 0.0;

   Gen.avg_fitness = 0.0;
   Gen.std_dev = 0.0;
   Gen.best_indv_index = -1;
   Gen.best_fitness = 0.0;
   Gen.worst_indv_index = -1;
   Gen.worst_fitness = 0.0;

  /* trace stuff for start of each generation */
   if (Trace == 1)  if (trace_gen_start() == ERROR)  return ERROR;

  /* initialize best individual of run */
   if (init_run_best_indv() == ERROR)  return ERROR;

  /* print params if file is on */
   if (file_on("params"))
      {
      ptr = get_file_pointer("params");
      Output_file[ptr].fp = fopen(Output_file[ptr].filename, "a");
      print_params(Output_file[ptr].fp);
      fclose(Output_file[ptr].fp);
      }  /* if */

   if (!strcmp(Base, "multichar"))
      {
     /* 000401AW  for keeping track of proportion of xters */
      Gen.avg_xter_proportion = (double *)malloc(Alphabet_size*sizeof(double));
     /* 031227AW  for keeping track of count of xters */
      Gen.avg_xter_count = (double *)malloc(Alphabet_size*sizeof(double));
     /* 000401AW: init xter proportion to zero.  Will be summed in eval_indv()
        for each individual, then divide by pop_size in gen_end().
        031227AW: xter count initialized here and summed in eval_indv() */
      for (i=Alphabet_size-1; i>=0; i--)
         {
         Gen.avg_xter_proportion[i] = 0.0;
         Gen.avg_xter_count[i] = 0.0;
         }
      }  /* if */

  /* 031208AW */
  /* if xtrack file is on, allocate space to store xtrack data for best indv */
   if (file_on("xtrack"))
      {
      Gen.avg_xtrack = (XTER_TRACK *)malloc(sizeof(XTER_TRACK));
      Gen.avg_xtrack->dist_sum = 0;
      Gen.avg_xtrack->dist_avg = -1;
      Gen.avg_xtrack->dist_sd = -1;
      Gen.avg_xtrack->dist_max = INT_MIN;
      Gen.avg_xtrack->dist_min = INT_MAX;
      Gen.best_xtrack = (XTER_TRACK *)malloc(sizeof(XTER_TRACK));
      Gen.best_xtrack->dist_sum = 0;
      Gen.best_xtrack->dist_avg = -1;
      Gen.best_xtrack->dist_sd = -1;
      Gen.best_xtrack->dist_max = INT_MIN;
      Gen.best_xtrack->dist_min = INT_MAX;
      }  /* if */

#ifdef DEBUG
   printf(" ---end run_start---\n");
#endif
   return OK;
   }  /* run_start */

/********** run_end **********/
/* parameters:
   called by:
   actions:     Things to be done once at end of run.     
*/
int run_end()
   {
#ifdef DEBUG
   printf(" ---in run_end---\n");
#endif

   run_output();

#ifdef DEBUG
   printf(" ---end run_end---\n");
#endif
   return OK;
   }  /* run_end */

/********** gen_start **********/
/* parameters:
   called by:
   actions:     Things to be done at start of each generation
                Such as incrementing gen.index and initializing the
		rest of the gen.* values.
		Subset of run_start() -- generation 0 doesn't call
		gen_start(), but calls run_start() instead.
*/
int gen_start()
   {
   int i;

#ifdef DEBUG
   printf(" ---in gen_start---\n");
#endif

   Prev_gen = Gen;

   Gen.index++;
   Gen.fitness_sum = 0;

   Gen.avg_fitness = 0.0;
   Gen.std_dev = 0.0;
   Gen.best_indv_index = -1;
   Gen.best_fitness = 0.0;
   Gen.worst_indv_index = -1;
   Gen.worst_fitness = 0.0;

   if (!strcmp(Base, "multichar"))
      {
     /* 000401AW: init xter count to zero.  Will be summed in eval_indv()
        for each individual, then divide by pop_size in gen_end */
      for (i=Alphabet_size-1; i>=0; i--)  Gen.avg_xter_proportion[i] = 0.0;
      }

  /* trace stuff for start of each generation */
   if (Trace == 1)  if (trace_gen_start() == ERROR)  return ERROR;

  /* function specific stuff  030528.AW */
   if (fxn_gen_start() == ERROR)  return ERROR;

  /* initialize xtrack values for generation if flag on */
   if (file_on("xtrack"))  xtrack_gen_start();

#ifdef DEBUG
   printf(" ---end gen_start---\n");
#endif
   return OK;
   }  /* gen_start */
 
/********** gen_end **********/
/* parameters:
   called by:
   actions:     Things to be done at end of each generation.     
		Such as calculating statistics.
*/
int gen_end()
   {
   int i;
   double sum, sum2;
   int len_sum, len_sum2;

#ifdef DEBUG
   printf(" ---in gen_end---\n");
#endif

  /* init for fitness stats */
   Gen.best_fitness = -DBL_MAX;
   Gen.worst_fitness = DBL_MAX;
   Gen.best_indv_index = Gen.worst_indv_index = -1;
   sum = sum2 = 0.0;
  /* init for length stats */
   Gen.longest = INT_MIN;
   Gen.shortest = INT_MAX;
   Gen.longest_index = Gen.shortest_index = -1;
   len_sum = len_sum2 = 0;

   for (i=Pop_size-1; i>=0; i--)
      {
     /* check for best/worst fitnesses of population */
      if (Pop[i]->fitness > Gen.best_fitness)
         {
         Gen.best_indv_index = i;
         Gen.best_fitness = Pop[i]->fitness;
         }  /* if */
      if (Pop[i]->fitness < Gen.worst_fitness)
         {
         Gen.worst_indv_index = i;
         Gen.worst_fitness = Pop[i]->fitness;
         }  /* if */

     /* check for longest/shortest of population */
      if (Pop[i]->length > Gen.longest)
         {
         Gen.longest_index = i;
         Gen.longest = Pop[i]->length;
         }  /* if */
      if (Pop[i]->length < Gen.shortest)
         {
         Gen.shortest_index = i;
         Gen.shortest = Pop[i]->length;
         }  /* if */

     /* calculations */
      sum += Pop[i]->fitness;
      sum2 += Pop[i]->fitness * Pop[i]->fitness;
      len_sum += Pop[i]->length;
      len_sum2 += Pop[i]->length * Pop[i]->length;
      }  /* for i */

  /* fitness: calculate average and standard deviation */
   Gen.avg_fitness = sum/(double)Pop_size;
   Gen.std_dev = sqrt(fabs(sum2 - sum*sum/(double)Pop_size)/
			(double)(Pop_size - 1) );
   Gen.fitness_sum = sum;

  /* length: calculate average and standard deviation */
   Gen.len_avg = (double)len_sum/(double)Pop_size;
   Gen.len_std = sqrt(fabs((double)len_sum2 -
			(double)len_sum*(double)len_sum/(double)Pop_size)/
			(double)(Pop_size - 1) );

  /* check Run_best_indv */
   if (Gen.best_fitness > Run_best_indv->fitness)
      copy_indv(Pop[Gen.best_indv_index], Run_best_indv);

   if (!strcmp(Base, "multichar"))
      {
     /* 000401AW: init xter count to zero.  Will be summed in eval_indv()
        for each individual, then divide by pop_size in gen_end */
      for (i=Alphabet_size-1; i>=0; i--)
         {
         Gen.avg_xter_proportion[i] =
            Gen.avg_xter_proportion[i]/(double)Pop_size;
         Gen.avg_xter_count[i] =
            Gen.avg_xter_count[i]/(double)Pop_size;
         }

#ifdef DEBUG_COUNT
      printf(" Generation averaged -- %d\n", Gen.index);
      for (i=0; i<Alphabet_size; i++)
         printf(" %lf", Gen.avg_xter_proportion[i]);
      printf("\n");
#endif
      }  /* if */

  /* calculate xtrack data for entire generation if file flag on */
   if (file_on("xtrack"))  xtrack_gen_end();

  /* output for each generation */
   gen_output();

  /* trace output for each generation */
   if (Trace == 1)  if (trace_gen_end() == ERROR)  return ERROR;

  /* CoG output for each generation */
   if (CoG == 1)  if (cog_gen_end() == ERROR)  return ERROR;

  /* function specific stuff  030528.AW */
   if (fxn_gen_end() == ERROR)  return ERROR;

#ifdef DEBUG
   printf(" ---end gen_end---\n");
#endif
   return OK;
   }  /* gen_end */

