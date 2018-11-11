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
#include <string.h>
#include <float.h>
#include "types.h"
#include "extern.h"
#include "indv.h"
#include "output.h"
#include "random.h"

#define DEBUG_COUNT 1
#undef DEBUG_COUNT

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
      indv2->genome[i] = indv1->genome[i];

  /* FXN_SPECIFIC_CALL */
  /* copy marker field only if FRR function */
   if (!strcmp(Function_name, "FRR"))
      for (i=0; i<indv1->length; i++)
         indv2->marker[i] = indv1->marker[i];
 
   if (Trace)
      {
      for (i=0; i<indv1->num_xover; i++)
         indv2->xover_pts[i] = indv1->xover_pts[i];
      for (i=0; i<indv1->num_mut; i++)
         indv2->mut_pts[i] = indv1->mut_pts[i];
      for (i=0; i<indv1->num_kids; i++)
         indv2->kids[i] = indv1->kids[i];
      }  /* if Trace */
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
   Run_best_indv->genome = (int *)malloc(Max_gen_len * sizeof(int));
   Run_best_indv->marker = (char *)malloc(Max_gen_len * sizeof(char));

  /* trace file */
   if (Trace)
      {
      Run_best_indv->num_xover = 0;
      Run_best_indv->num_mut = 0;
      Run_best_indv->num_kids = 0;

      Run_best_indv->xover_pts = (int *)malloc(Run_best_indv->length *
						sizeof(int));
      if (Run_best_indv->xover_pts == NULL)
         {
         printf(" Error(malloc_indv): cannot allocate space:");
         printf(" indv->xover_pts\n");
         return ERROR;
         }  /* if */
      Run_best_indv->mut_pts = (int *)malloc(Run_best_indv->length*sizeof(int));
      if (Run_best_indv->mut_pts == NULL)
         {
         printf(" Error(malloc_indv): cannot allocate space:");
         printf(" indv->mut_pts\n");
         return ERROR;
         }  /* if */
      Run_best_indv->kids = (int *)malloc(Pop_size * sizeof(int));
      if (Run_best_indv->kids == NULL)
         {
         printf(" Error(malloc_indv): cannot allocate space:");
         printf(" indv->kids\n");
         return ERROR;
         }  /* if */
      }  /* if Trace */

   init_indv(Run_best_indv, -1);
   Run_best_indv->fitness = -DBL_MAX;
   Run_best_indv->raw_fitness = -DBL_MAX;

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
   int i;

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
   if (Trace)
      {
     /* xover locations */
      fprintf(fp, " num_xover = %d", indv->num_xover);
      if (indv->num_xover > 0)
         {
         fprintf(fp, " at");
         for (i=0; i<indv->num_xover; i++)
            fprintf(fp, " %d", indv->xover_pts[i]);
         }
      fprintf(fp, "\n");
     /* mutation locations */
      fprintf(fp, " num_mut = %d", indv->num_mut);
      if (indv->num_mut > 0)
         {
         fprintf(fp, " at");
         for (i=0; i<indv->num_mut; i++)
            fprintf(fp, " %d", indv->mut_pts[i]);
         }
      fprintf(fp, "\n");
     /* offspring */
      fprintf(fp, " num_kids = %d", indv->num_kids);
      if (indv->num_kids > 0)
         {
         fprintf(fp, " indices");
         for (i=0; i<indv->num_kids; i++)
            fprintf(fp, " %d", indv->kids[i]);
         }
      fprintf(fp, "\n");
      }  /* if Trace */
   fprintf(fp, " ");
   fprint_genome(fp, indv->genome, indv->length, 1);

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
  /* trace file */
   if (Trace)
      {
      free(indv->xover_pts);
      free(indv->mut_pts);
      free(indv->kids);
      }  /* if Trace */

   free(indv->genome);
   if (indv->marker != NULL)  free(indv->marker);
   free(indv);
   }  /* free_indv */

/********** indv_xter_count **********/
/* parameters:	indv
   created:	000401AW
   called by:	eval_indv(), fxsim.c
   actions:	Count the number of each xter in the individual.
		Only called for multichar alphabets.
*/
void indv_xter_count(INDIVIDUAL *indv)
   {
   int i, j;

#ifdef DEBUG
printf(" ---in indv_xter_count()---\n");
#endif

   for (i=Alphabet_size-1; i>=0; i--)  indv->xter_count[i] = 0;

   for (j=indv->length-1; j>=0; j--)
      {
      for (i=Alphabet_size-1; i>=0; i--)
         {
         if (indv->genome[j] == Xters[i])
            {
            indv->xter_count[i]++;
            break;
            }  /* if */
         }  /* fir */
      }  /* for j */

#ifdef DEBUG_COUNT
   print_genome(indv->genome, indv->length, 1);
   for (i=0; i<Alphabet_size; i++)
      {
      printf(" %d", indv->xter_count[i]);
      }
   printf("\n");
#endif

#ifdef DEBUG
printf(" ---end indv_xter_count()---\n");
#endif

   }  /* indv_xter_count */
