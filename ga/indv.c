/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* indv.c
   07.28.98.AW	Created.  Routines having to do with individuals.

   Routines:	copy_indv()		07.28.98.AW
		init_run_best_indv()	07.29.98.AW
		init_indv()		07.29.98.AW
		print_indv()		07.29.98.AW
		fprint_indv()		07.29.98.AW
		free_indv()		07.29.98.AW
*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "types.h"
#include "extern.h"
#include "indv.h"
#include "output.h"

/********** copy_indv **********/
/* parameters:  indv1
		indv2
   called by:	run_start(), stats.c
   actions:     copies indv1 to indv2.  Assumes everything
		is allocated alright.
*/
void copy_indv(INDIVIDUAL *indv1, INDIVIDUAL *indv2)
   {
   int i;
 
#ifdef DEBUG
   printf(" ---in copy_indv---\n");
#endif
   indv2->index = indv1->index;
   indv2->gen = indv1->gen;
   indv2->length = indv1->length;
   indv2->fitness = indv1->fitness;
   indv2->raw_fitness = indv1->raw_fitness;
   indv2->calc_num_offspring = indv1->calc_num_offspring;
   indv2->parent1_index = indv1->parent1_index;
   indv2->parent2_index = indv1->parent2_index;
   indv2->chosen = indv1->chosen;
   indv2->can_mate = indv1->can_mate;
   indv2->num_xover = indv1->num_xover;
   indv2->num_mut = indv1->num_mut;
   indv2->num_kids = indv1->num_kids;

   for (i=0; i<indv1->length; i++)
      {
      indv2->genome[i] = indv1->genome[i];
      if (Init_pop == 2) 
         {
         indv2->floats_genome[i] = indv1->floats_genome[i];
         }
      }
 
#ifdef DEBUG
   printf(" ---end copy_indv---\n");
#endif
   }  /* copy_indv */

/********** init_run_best_indv **********/
/* parameters:  indv
   called by:
   actions:     allocate space for Run_best_indv.
*/
int init_run_best_indv()
   {
   Run_best_indv = (INDIVIDUAL *)malloc(sizeof(INDIVIDUAL));
   if (Run_best_indv == NULL)
      {
      printf(" Error(malloc_indv): cannot allocate space: Run_best_indv \n");
      return ERROR;
      }

   if (Variable_gen_len == 0)
      {
      Run_best_indv->length = Max_gen_len;
      }  /* if */
   else if (Variable_gen_len > 1 && Variable_gen_len <= Max_gen_len)
      {
      Run_best_indv->length = Variable_gen_len;
      }  /* else if */
   else if (Variable_gen_len == 1)
      {
      Run_best_indv->length = Min_gen_len + uniform(Max_gen_len-Min_gen_len+1);
      }  /* else if */
   else
      {
      printf(" Error(init_run_best_indv): invalid value: Variable_gen_len: %d",
             Variable_gen_len);
      printf("\n");
      return ERROR;
      }  /* else */
   Run_best_indv->floats_genome = (double *)malloc(Max_gen_len * sizeof(double));
   Run_best_indv->genome = (double *)malloc(Max_gen_len * sizeof(int));
   Run_best_indv->marker = (char *)malloc(Max_gen_len * sizeof(char));

   init_indv(Run_best_indv, -1);
   Run_best_indv->fitness = DBL_MAX;
   Run_best_indv->raw_fitness = DBL_MAX;

   return OK;
   }  /* init_run_best_indv */

/********** init_indv **********/
/* parameters:  indv
   called by:	init_run_best_indv() 
		pop_eval(), pop.c
   actions:     initializes and allocates space for individual.
		Does not initialize parent indexes because this
		is always done by crossover.
*/
void init_indv(INDIVIDUAL *indv, int gen)
   {
   indv->gen = gen;
   indv->chosen = 0;
   indv->can_mate = 1;
   indv->fitness = 0.0;
   indv->raw_fitness = 0.0;
   indv->calc_num_offspring = -1.0;

   indv->num_kids = 0;
  /* num_xover and num_mut are set every generation in genops.c
   so don't need to init */
   }  /* init_indv */

/********** print_indv **********/
/* parameters:  indv
   called by:
   actions:     prints data from one individual
*/
void print_indv(INDIVIDUAL *indv)
   {
   if (indv->gen == -1)  printf(" Individual %d from gen %d\n", indv->index, Gen.index);
   else  printf(" Individual %d from gen %d\n", indv->index, indv->gen);

   printf("          length = %d\n", indv->length);
   printf("          chosen = %d\n", indv->chosen);
   printf("          can_mate = %d\n", indv->can_mate);
   printf("          fitness = %lf\n", indv->fitness);
   printf("          calc_num_offspring = %lf\n", indv->calc_num_offspring);
   printf("          parent1_index = %d\n", indv->parent1_index);
   printf("          parent2_index = %d\n", indv->parent2_index);
   }  /* print_indv */

/********** fprint_indv **********/
/* parameters:  fp
		indv
   called by:
   actions:     about the same as print_indv only prints to a file
		and prints a little bit more stuff.
		** must check if trace file is turned on before trying
		to print xover and mutation information.
*/
void fprint_indv(FILE *fp, INDIVIDUAL *indv)
   {
   int i;

#ifdef DEBUG
   printf(" ---in fprint_indv---\n");
#endif

   fprintf(fp," Individual %d from generation %d\n",indv->index,indv->gen);

   fprintf(fp, " length = %d\n", indv->length);
   fprintf(fp, " chosen = %d\n", indv->chosen);
   fprintf(fp, " can_mate = %d\n", indv->can_mate);
   fprintf(fp, " fitness = %lf\n", indv->fitness);
   fprintf(fp, " raw_fitness = %lf\n", indv->raw_fitness);
   fprintf(fp, " calc_num_offspring = %lf\n",indv->calc_num_offspring);
   fprintf(fp, " parent1_index = %d\n", indv->parent1_index);
   fprintf(fp, " parent2_index = %d\n", indv->parent2_index);
   fprintf(fp, " ");
   fprint_genome(fp, indv, 1);

#ifdef DEBUG
   printf(" ---end fprint_indv---\n");
#endif
   }  /* fprint_indv */

/********** free_indv **********/
/* parameters:  indv
   called by:	not called by close_pop(), pop.c, right now because
		trying to minimize the number of checks for "trace".
   actions:     frees up individual.
*/
void free_indv(INDIVIDUAL *indv)
   {
   free(indv->floats_genome);
   if (indv->marker != NULL)  free(indv->marker);
   free(indv);
   }  /* free_indv */
