/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* stats.c
   07.22.98.AW	Created.  Routines that maintain and calculate
		statistics during a run.

   Routines:	run_start()	07.22.98.AW
		run_end()	07.22.98.AW
		gen_start()	07.22.98.AW
		gen_end()	07.22.98.AW
		gen_stats()	11.27.04.AW	Need to calculate generation
				statistics before calling reproduce().
				gen_stats() is bulk of what used to be in
				gen_end().  gen_end() mainly outputs stuff now.
*/

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "types.h"
#include "extern.h"
#include "stats.h"
#include "indv.h"
#include "params.h"
#include "output.h"

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

  /* initialize best individual of run */
   if (init_run_best_indv() == ERROR)  return ERROR;

  /* print params if file is on */
   if (file_on("params"))
      {
      ptr = get_file_pointer("params");
printf("run_start: %s\n", Output_file[ptr].filename);
      Output_file[ptr].fp = fopen(Output_file[ptr].filename, "a");
      print_params(Output_file[ptr].fp);
      fclose(Output_file[ptr].fp);
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

#ifdef DEBUG
   printf(" ---end gen_start---\n");
#endif
   return OK;
   }  /* gen_start */
 
/********** gen_stats **********/
/* parameters:
   called by:	ga_loop(), ga.c
   actions:     Calculating statistics for this generation.
*/
int gen_stats()
   {
   int i;
   double sum, sum2;
   int len_sum, len_sum2;

#ifdef DEBUG
   printf(" ---in gen_stats---\n");
#endif

  /* init for fitness stats */
   Gen.best_fitness = DBL_MAX;
   Gen.worst_fitness = -DBL_MIN;
   Gen.best_indv_index = Gen.worst_indv_index = -1;
   sum = sum2 = 0.0;
  /* init for length stats */
   Gen.longest = REALLY_BIG_NUMBER;
   Gen.shortest = REALLY_SMALL_NUMBER;
   Gen.longest_index = Gen.shortest_index = -1;
   len_sum = len_sum2 = 0;

   for (i=Pop_size-1; i>=0; i--)
      {
     /* check for best/worst fitnesses of population */
      if (Pop[i]->fitness < Gen.best_fitness)
         {
         Gen.best_indv_index = i;
         Gen.best_fitness = Pop[i]->fitness;
         }  /* if */
      if (Pop[i]->fitness > Gen.worst_fitness)
         {
         Gen.worst_indv_index = i;
         Gen.worst_fitness = Pop[i]->fitness;
         }  /* if */

     /* check for longest/shortest of population */
      if (Pop[i]->length < Gen.longest)
         {
         Gen.longest_index = i;
         Gen.longest = Pop[i]->length;
         }  /* if */
      if (Pop[i]->length > Gen.shortest)
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
   if (Gen.best_fitness < Run_best_indv->fitness)
      copy_indv(Pop[Gen.best_indv_index], Run_best_indv);

#ifdef DEBUG
   printf(" ---end gen_stats---\n");
#endif
   return OK;
   }  /* gen_stats */
 
/********** gen_end **********/
/* parameters:
   called by:	ga_loop(), ga.c
   actions:     Things to be done at end of each generation.     
		Mainly outputting stuff to files.
*/
int gen_end()
   {
#ifdef DEBUG
   printf(" ---in gen_end---\n");
#endif

  /* output for each generation */
   gen_output();

#ifdef DEBUG
   printf(" ---end gen_end---\n");
#endif
   return OK;
   }  /* gen_end */
