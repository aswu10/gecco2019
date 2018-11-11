/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* pop.c
   07.03.98.AW	Created.
		Routines that have to do with populations.
		Such as setting up the data structure and initializing it.

   Routines:	init_pop()			07.03.98.AW
		pop_allocate()			07.07.98.AW
		pop_fill()			07.07.98.AW
		indv_fill()			07.08.98.AW
		pop_individuals_0()		07.07.98.AW
		pop_individuals_1()		07.07.98.AW
		pop_individuals_random()	07.07.98.AW
		pop_individuals_readin()	07.08.98.AW
		pop_individuals_random_alphabet() 10.09.98.AW
		pop_individuals_random_multichar() 03.08.01.AW
		pop_eval()			07.08.98.AW
   		close_pop()			07.03.98.AW
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "types.h"
#include "extern.h"
#include "pop.h"
#include "params.h"
#include "random.h"
#include "indv.h"
#include "fxn.h"
/********** init_pop **********/
/* parameters:
   called by:	ga_init(), ga.c
   actions:	allocates space for population and gets
		the initial population.
*/
int init_pop()
   {
#ifdef DEBUG
   printf(" ---in init_pop()---\n");
#endif

   if (pop_allocate() == ERROR)  return ERROR;
   if (pop_fill() == ERROR)  return ERROR;

#ifdef DEBUG
   printf(" ---end init_pop()---\n");
#endif
   return OK;
   }  /* init_pop */

/********** pop_allocate **********/
/* parameters:
   called by:
   actions:	allocate space for populations: pop and kids.
		if Variable_gen_len = -1, don't allocate genome space,
			length = 0;
		if Variable_gen_len = 0, use fixed genome length given by
			Max_gen_len.
		if Variable_gen_len > 1, the run uses variable length
			genomes that range between Min_gen_len and Max_gen_len,
			initialize the individuals of the initial population
			all to the same length -- Variable_gen_len.
		if Variable_gen_len = 1, the run uses variable length
			genomes that range between Min_gen_len and Max_gen_len,
			initialize the individuals of the initial population
			randomly to lengths within acceptable range.
*/
int pop_allocate()
   {
   int i;
#ifdef DEBUG
   printf(" ---in pop_allocate()---\n");
#endif

  /* current population */
   Pop = (POPULATION)malloc(Pop_size * sizeof(INDIVIDUAL *));
   if (Pop == NULL)
      {
      printf(" Error(pop_allocate): cannot allocate space: Pop\n");
      return ERROR;
      }  /* if */

  /* next population */
   Kids = (POPULATION)malloc(Pop_size * sizeof(INDIVIDUAL *));
   if (Kids == NULL)
      {
      printf(" Error(pop_allocate): cannot allocate space: Kids\n");
      return ERROR;
      }  /* if */

  /* parent population (pointers) */
   Parents = (POPULATION)malloc(Pop_size * sizeof(INDIVIDUAL *));
   if (Parents == NULL)
      {
      printf(" Error(pop_allocate): cannot allocate space: Parents\n");
      return ERROR;
      }  /* if */

   for (i=0; i<Pop_size; i++)
      {
      Pop[i] = (INDIVIDUAL *)malloc(sizeof(INDIVIDUAL));
      if (Pop[i] == NULL)
         {
         printf(" Error(pop_allocate): cannot allocate space: pop[%d]\n",i);
         return ERROR;
         }  /* if */
      Kids[i] = (INDIVIDUAL *)malloc(sizeof(INDIVIDUAL));
      if (Kids[i] == NULL)
         {
         printf(" Error(pop_allocate): cannot allocate space: kids[%d]\n",i);
         return ERROR;
         }  /* if */

      Pop[i]->index = i;		Kids[i]->index = i;
      Pop[i]->gen = -1;		 	Kids[i]->gen = -1;
      Pop[i]->chosen = 0;		Kids[i]->chosen = 0;
      Pop[i]->can_mate = 1;		Kids[i]->can_mate = 1;

      Pop[i]->fitness = -1.0;		Kids[i]->fitness = -1.0;
      Pop[i]->calc_num_offspring= -1.0;	Kids[i]->calc_num_offspring = -1.0;
      Pop[i]->parent1_index = -1;	Kids[i]->parent1_index = -1;
      Pop[i]->parent2_index = -1;	Kids[i]->parent2_index = -1;

      if (Variable_gen_len == 0)
         {
			 printf("%d\n", Max_gen_len);
         Pop[i]->length = Max_gen_len;
         Kids[i]->length = Max_gen_len;
         }  /* if */
      else if (Variable_gen_len > 1 && Variable_gen_len <= Max_gen_len)
         {
/* **********************************************************************

  06.07.01.IG
          
  ***********************************************************************

  Pop[i]->length = Min_gen_len + uniform(Variable_gen_len - Min_gen_len +1);
  Kids[i]->length = 0;


 *********************************************************************	  
*/
	 Pop[i]->length = Variable_gen_len;
         Kids[i]->length = Variable_gen_len;
		 
         }  /* else if */
      else if (Variable_gen_len == 1)
         {
         Pop[i]->length = Init_min_gen_len +
		uniform(Init_max_gen_len - Init_min_gen_len +1);
         Kids[i]->length = 0;
         }  /* else if */
      else
         {
         printf(" Error(pop_allocate): invalid value: how/Variable_gen_len: %d",
		Variable_gen_len);
         printf("\n");
         return ERROR;
         }  /* else */

      Pop[i]->genome = (int *)malloc(Max_gen_len * sizeof(int));
      Kids[i]->genome = (int *)malloc(Max_gen_len * sizeof(int));
     /* 000401AW */
      if (!strcmp(Base, "multichar"))
         {
         Pop[i]->xter_count = (int *)malloc(Alphabet_size * sizeof(int));
         Kids[i]->xter_count = (int *)malloc(Alphabet_size * sizeof(int));
         }

     /* trace file */
      if (Trace)
         {
         Pop[i]->num_xover = 0;		Kids[i]->num_xover = 0;
         Pop[i]->num_mut = 0;		Kids[i]->num_mut = 0;
         Pop[i]->num_kids = 0;		Kids[i]->num_kids = 0;
         Pop[i]->num_segments = 0;	Kids[i]->num_segments = 0;

         Pop[i]->xover_pts = (int *)malloc(Max_gen_len * sizeof(int));
         if (Pop[i]->xover_pts == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Pop[i]->xover_pts\n");
            return ERROR;
            }  /* if */
         Pop[i]->mut_pts = (int *)malloc(Max_gen_len * sizeof(int));
         if (Pop[i]->mut_pts == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Pop[i]->mut_pts\n");
            return ERROR;
            }  /* if */
         Pop[i]->kids = (int *)malloc(Pop_size * sizeof(int));
         if (Pop[i]->kids == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Pop[i]->kids\n");
            return ERROR;
            }  /* if */
         Pop[i]->seg_parent = (int *)malloc((Max_gen_len+1) * sizeof(int));
         if (Pop[i]->seg_parent == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Pop[i]->seg_parent\n");
            return ERROR;
            }  /* if */
         Pop[i]->seg_start = (int *)malloc((Max_gen_len+1) * sizeof(int));
         if (Pop[i]->seg_start == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Pop[i]->seg_start\n");
            return ERROR;
            }  /* if */
         Pop[i]->seg_len = (int *)malloc((Max_gen_len+1) * sizeof(int));
         if (Pop[i]->seg_len == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Pop[i]->seg_len\n");
            return ERROR;
            }  /* if */
         Kids[i]->xover_pts = (int *)malloc(Max_gen_len * sizeof(int));
         if (Kids[i]->xover_pts == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Kids[i]->xover_pts\n");
            return ERROR;
            }  /* if */
         Kids[i]->mut_pts = (int *)malloc(Max_gen_len * sizeof(int));
         if (Kids[i]->mut_pts == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Kids[i]->mut_pts\n");
            return ERROR;
            }  /* if */
         Kids[i]->kids = (int *)malloc(Pop_size * sizeof(int));
         if (Kids[i]->kids == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Kids[i]->kids\n");
            return ERROR;
            }  /* if */
         Kids[i]->seg_parent = (int *)malloc((Max_gen_len+1) * sizeof(int));
         if (Kids[i]->seg_parent == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Kids[i]->seg_parent\n");
            return ERROR;
            }  /* if */
         Kids[i]->seg_start = (int *)malloc((Max_gen_len+1) * sizeof(int));
         if (Kids[i]->seg_start == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Kids[i]->seg_start\n");
            return ERROR;
            }  /* if */
         Kids[i]->seg_len = (int *)malloc((Max_gen_len+1) * sizeof(int));
         if (Kids[i]->seg_len == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Kids[i]->seg_len\n");
            return ERROR;
            }  /* if */
         }  /* if Trace */

     /* CoG file */
      if (CoG)
         {
         Pop[i]->cog_count = (int *)malloc(Alphabet_size * sizeof(int));
         if (Pop[i]->cog_count == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Pop[i]->cog_count\n");
            return ERROR;
            }  /* if */
         Pop[i]->cog_avg = (double *)malloc(Alphabet_size * sizeof(double));
         if (Pop[i]->cog_avg == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Pop[i]->cog_avg\n");
            return ERROR;
            }  /* if */
         Pop[i]->cog_sd = (double *)malloc(Alphabet_size * sizeof(double));
         if (Pop[i]->cog_sd == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Pop[i]->cog_sd\n");
            return ERROR;
            }  /* if */
         Pop[i]->cog_low = (double *)malloc(Alphabet_size * sizeof(double));
         if (Pop[i]->cog_low == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Pop[i]->cog_low\n");
            return ERROR;
            }  /* if */
         Pop[i]->cog_high = (double *)malloc(Alphabet_size * sizeof(double));
         if (Pop[i]->cog_high == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Pop[i]->cog_high\n");
            return ERROR;
            }  /* if */
         Kids[i]->cog_count = (int *)malloc(Alphabet_size * sizeof(int));
         if (Kids[i]->cog_count == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Kids[i]->cog_count\n");
            return ERROR;
            }  /* if */
         Kids[i]->cog_avg = (double *)malloc(Alphabet_size * sizeof(double));
         if (Kids[i]->cog_avg == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Kids[i]->cog_avg\n");
            return ERROR;
            }  /* if */
         Kids[i]->cog_sd = (double *)malloc(Alphabet_size * sizeof(double));
         if (Kids[i]->cog_sd == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Kids[i]->cog_sd\n");
            return ERROR;
            }  /* if */
         Kids[i]->cog_low = (double *)malloc(Alphabet_size * sizeof(double));
         if (Kids[i]->cog_low == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Kids[i]->cog_low\n");
            return ERROR;
            }  /* if */
         Kids[i]->cog_high = (double *)malloc(Alphabet_size * sizeof(double));
         if (Kids[i]->cog_high == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Kids[i]->cog_high\n");
            return ERROR;
            }  /* if */
         }  /* if CoG */

     /* xter track file on */
      if (file_on("xtrack"))
         {
         Pop[i]->xtrack = (XTER_TRACK *)malloc(sizeof(XTER_TRACK));
         if (Pop[i]->xtrack == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Pop[i]->xtrack\n");
            return ERROR;
            }  /* if */
         Kids[i]->xtrack = (XTER_TRACK *)malloc(sizeof(XTER_TRACK));
         if (Kids[i]->xtrack == NULL)
            {
            printf(" Error(pop_allocate): cannot allocate space:");
            printf(" Kids[i]->xtrack\n");
            return ERROR;
            }  /* if */

        /* init */
         Pop[i]->xtrack->dist_sum = -1;
         Pop[i]->xtrack->dist_avg = -1;
         Pop[i]->xtrack->dist_sd = -1;
         Pop[i]->xtrack->dist_max = INT_MIN;
         Pop[i]->xtrack->dist_min = INT_MAX;
         Pop[i]->xtrack->num = -1;
         Pop[i]->xtrack->dist_sum2 = -1;

         Kids[i]->xtrack->dist_sum = -1;
         Kids[i]->xtrack->dist_avg = -1;
         Kids[i]->xtrack->dist_sd = -1;
         Kids[i]->xtrack->dist_max = INT_MIN;
         Kids[i]->xtrack->dist_min = INT_MAX;
         Kids[i]->xtrack->num = -1;
         Kids[i]->xtrack->dist_sum2 = -1;
         }  /* if xtrack file on */

      }  /* for i */

#ifdef DEBUG
   printf(" ---end pop_allocate()---\n");
#endif
   return OK;
   }  /* pop_allocate */

/********** indv_fill **********/
/* parameters:	genome		to fill
		len		of genome
		x		character to fill with
   called by:	
   actions:	fills a genome with a given character
*/
void indv_fill(int *genome, int len, int x)
   {
   int i;
   for (i=len-1; i>=0; i--)
      genome[i] = x;
   }  /* indv_fill */

/********** pop_fill **********/
/* parameters:  how             what to initialize individuals to:
                                -1 = don't generate initial members, just
                             _       set up structure
                            |   0 = initial indvs all 0's
                 Init_pop  _|   1 = initial indvs all 1's
                            |   2 = initial indvs random
                            |_  3 = read initial indvs from Init_pop_file

   called by:	init_pop(), pop.c
   actions:     get initial individuals for population.
*/
int pop_fill()
   {
#ifdef DEBUG
   printf(" ---in pop_fill()---\n");
#endif

   if (!strcmp(Base, "binary"))
      {
      switch(Init_pop)
         {
         case -1:	break;
         case 0:	pop_individuals_0();
   			break;
         case 1:	pop_individuals_1();
			break;
         case 2:	pop_individuals_random();
			break;
         case 3:	if (pop_individuals_readin() == ERROR)  return ERROR;
			break;
         default:	printf(" Error(pop_fill): Init_pop value invalid");
			printf(" for binary base: %d\n", Init_pop);
			return ERROR;
         }  /* switch */
      }  /* if */
   else if (!strcmp(Base, "alphabet"))
      {
      switch(Init_pop)
         {
         case -1:	break;
         case 3:	if (pop_individuals_readin() == ERROR)  return ERROR;
			break;
         case 4:	pop_individuals_random_alphabet();
   			break;
         default:	printf(" Error(pop_fill): Init_pop value invalid");
			printf(" for alphabet base: %d\n", Init_pop);
			return ERROR;
         }  /* switch */
      }  /* else if */
   else if (!strcmp(Base, "multichar"))
      {
      switch(Init_pop)
         {
         case -1:       break;
         case 3:        if (pop_individuals_readin() == ERROR)  return ERROR;
                        break;
         case 4:        pop_individuals_random_multichar();
                        break;
         default:       printf(" Error(pop_fill): Init_pop value invalid");
                        printf(" for multichar base: %d\n", Init_pop);
                        return ERROR;
         }  /* switch */
      }  /* else if */
   else if (!strcmp(Base, "integers"))
      {
      switch(Init_pop)
         {
         case -1:       break;
         case 3:        if (pop_individuals_readin() == ERROR)  return ERROR;
                        break;
         case 4:        pop_individuals_random_int();
                        break;
         default:       printf(" Error(pop_fill): Init_pop value invalid");
                        printf(" for integers base: %d\n", Init_pop);
                        return ERROR;
         }  /* switch */
      }  /* else if */

#ifdef DEBUG
   printf(" ---end pop_fill()---\n");
#endif
   return OK;
   }  /* pop_fill */

/********** pop_individuals_0 **********/
/* parameters:
   called by:   pop_fill(), pop.c
   actions:     fills population with individuals made up of all 0's.
*/
void pop_individuals_0()
   {
   int i, j;

#ifdef DEBUG
   printf(" ---in pop_individuals_0()---\n");
#endif

   for (i=Pop_size-1; i>=0; i--)
      {
      for (j=Pop[i]->length-1; j>=0; j--)
         {
         Pop[i]->genome[j] = '0';
         }  /* for j */
      }  /* for i */
 
#ifdef DEBUG
   printf(" ---end pop_individuals_0()---\n");
#endif
   }  /* pop_individuals_0 */

/********** pop_individuals_1 **********/
/* parameters:
   called by:   pop_fill(), pop.c
   actions:     fills population with individuals made up of all 1's.
*/
void pop_individuals_1()
   {
   int i, j;

#ifdef DEBUG
   printf(" ---in pop_individuals_1()---\n");
#endif
 
   for (i=Pop_size-1; i>=0; i--)
      {
      for (j=Pop[i]->length-1; j>=0; j--)
         {
         Pop[i]->genome[j] = '1';
         }  /* for j */
      }  /* for i */

#ifdef DEBUG
   printf(" ---end pop_individuals_1()---\n");
#endif
   }  /* pop_individuals_1 */

/********** pop_individuals_random **********/
/* parameters:
   called by:   pop_fill(), pop.c
   actions:     fills population with individuals of random 1's and 0's.
*/
void pop_individuals_random()
   {
   int i, j;

#ifdef DEBUG
   printf(" ---in pop_individuals_random()---\n");
#endif
 
   for (i=Pop_size-1; i>=0; i--)
      {
      for (j=Pop[i]->length-1; j>=0; j--)
         {
         if (uniform(2))  Pop[i]->genome[j] = '0';
         else  Pop[i]->genome[j] = '1';
         }  /* for j */
      }  /* for i */

#ifdef DEBUG
   printf(" ---end pop_individuals_random()---\n");
#endif
   }  /* pop_individuals_random */

/********** pop_individuals_readin **********/
/* parameters:
   called by:   pop_fill(), pop.c
   actions:     fills population with individuals read in from 
		Init_pop_file.
   09.08.98	If the population size read in here does not match Pop_size,
		return error instead of resetting.  Resetting could be a
		problem if not enough space has been allocated to Pop (done
		in pop_allocate() which is called before this routine).
*/
int pop_individuals_readin()
   {
   FILE *fp;
   int i, j;

#ifdef DEBUG
   printf(" ---in pop_individuals_readin()---\n");
#endif

   fp = fopen(Init_pop_file, "r");
   if (fp == NULL)
      {
      printf(" Error(pop_indvs_readin): cannot open file: %s\n",Init_pop_file);
      return ERROR;
      }  /* if */
   printf(" Reading initial population from file: %s\n", Init_pop_file);
 
  /* population size is first item in file -- in case of differences
     between the Pop_size given in this file and that given in the params
     file, this one wins. */
   fscanf(fp, "%d %d", &i, &j);
   if (i != Pop_size)
      {
      printf(" Error(pop_indvs_readin): Pop_size inconsistent: %d (new)", i);
      printf(" %d (params file)\n", Pop_size);
      }  /* if */
   if (j != Variable_gen_len)
      {
      printf(" Note(pop_indvs_readin): Variable_gen_len reset from %d to %d\n",
		Variable_gen_len, j);
      Variable_gen_len = j;
      }  /* if */

   for (i=0; i<Pop_size; i++)
      {
      fscanf(fp, "%d", &Pop[i]->length);
     /* read eoln xter */
      getc(fp);
      for (j=0; j<Pop[i]->length; j++)
         {
         Pop[i]->genome[j] = getc(fp);
         }  /* for j */
     /* read eoln xter */
      getc(fp);
      }  /* for i */

   fclose(fp);
#ifdef DEBUG
   printf(" ---end pop_individuals_readin()---\n");
#endif
   return OK;
   }  /* pop_individuals_readin */

/********** pop_individuals_random_alphabet **********/
/* parameters:
   called by:   pop_fill(), pop.c
   actions:     fills population with individuals of random letters.
*/
void pop_individuals_random_alphabet()
   {
   int i, j;

#ifdef DEBUG
   printf(" ---in pop_individuals_random_alphabet()---\n");
#endif
 
   for (i=Pop_size-1; i>=0; i--)
      {
      for (j=Pop[i]->length-1; j>=0; j--)
         {
         Pop[i]->genome[j] = uniform(26) + Ascii_offset;
         }  /* for j */
      }  /* for i */

#ifdef DEBUG
   printf(" ---end pop_individuals_random_alphabet()---\n");
#endif
   }  /* pop_individuals_random_alphabet */

/********** pop_individuals_random_multichar **********/
/* parameters:
   called by:   pop_fill(), pop.c
   actions:     fills population with individuals of random letters.
*/
void pop_individuals_random_multichar()
   {
   int i, j;

#ifdef DEBUG
   printf(" ---in pop_individuals_random_multichar()---\n");
#endif

   for (i=Pop_size-1; i>=0; i--)
      {
      for (j=Pop[i]->length-1; j>=0; j--)
         {
         Pop[i]->genome[j] = Xters[uniform(Alphabet_size)];
         }  /* for j */
      }  /* for i */

#ifdef DEBUG
   printf(" ---end pop_individuals_random_multichar()---\n");
#endif
   }  /* pop_individuals_random_multichar */
   
/********** pop_individuals_random_int **********/
/* parameters:
   called by:   pop_fill(), pop.c
   actions:     fills population with individuals of random integers,
                both positive and negative.
*/
void pop_individuals_random_int()
   {
   int i, j;

#ifdef DEBUG
   printf(" ---in pop_individuals_random_int()---\n");
#endif

   for (i=Pop_size-1; i>=0; i--)
      {	  
      for (j=Pop[i]->length-1; j>=0; j--)
         {
         Pop[i]->genome[j] = 1 + uniform(Alphabet_size);    
         if (funiform(1) > 0.5) 
		    {
	        Pop[i]->genome[j] = -1 * Pop[i]->genome[j];		 
		    }
         }  /* for j */
      }  /* for i */

#ifdef DEBUG
   printf(" ---end pop_individuals_random_int()---\n");
#endif
   }  /* pop_individuals_random_int */

/********** pop_eval **********/
/* parameters:
   called by:
   actions:	evaluates the individuals in the population.
		also sets chosen (which is the actual # of offspring when
		parent replacement is off) to zero.
		and set fitness and raw_fitness to zero.

   990201.AW	Added check for Flat_fitness in calculating indv fitness.
		Need to run eval_indv() no matter what for the bbb_exists
		data.  But if Flat_fitness > 0, override actual fitness
		with Flat_fitness.
*/
void pop_eval()
   {
   int i;

#ifdef DEBUG
   printf(" ---in pop_eval()---\n");
#endif

   for (i=Pop_size-1; i>=0; i--)
      {
      init_indv(Pop[i], Gen.index);
      eval_indv(Pop[i]);

      if (Flat_fitness > 0)  Pop[i]->fitness = Flat_fitness;
      }  /* for i */

#ifdef DEBUG
   printf(" ---end pop_eval()---\n");
#endif
   }  /* pop_eval */

/********** close_pop **********/
/* parameters:
   called by:	ga_end(), ga.c
   actions:	whatever needs to be done to free up population structure.
*/
void close_pop()
   {
   int i;
#ifdef DEBUG
   printf(" ---in close_pop()---\n");
#endif

  /* trace */
   if (Trace)
      {
      for (i=0; i<Pop_size; i++)
         {
         free(Pop[i]->xover_pts);
         free(Pop[i]->mut_pts);
         free(Pop[i]->kids);
         free(Pop[i]->seg_parent);
         free(Pop[i]->seg_start);
         free(Pop[i]->seg_len);
         free(Kids[i]->xover_pts);
         free(Kids[i]->mut_pts);
         free(Kids[i]->kids);
         free(Kids[i]->seg_parent);
         free(Kids[i]->seg_start);
         free(Kids[i]->seg_len);
         }  /* for i */
      }  /* if Trace */

  /* cog */
   if (CoG)
      {
      for (i=0; i<Pop_size; i++)
         {
         free(Pop[i]->cog_count);
         free(Pop[i]->cog_avg);
         free(Pop[i]->cog_sd);
         free(Pop[i]->cog_low);
         free(Pop[i]->cog_high);
         free(Kids[i]->cog_count);
         free(Kids[i]->cog_avg);
         free(Kids[i]->cog_sd);
         free(Kids[i]->cog_low);
         free(Kids[i]->cog_high);
         }  /* for i */
      }  /* if CoG */

  /* multichar */
   if (!strcmp(Base, "multichar"))
      {
      for (i=0; i<Pop_size; i++)
         {
         free(Pop[i]->xter_count);
         free(Kids[i]->xter_count);
         } /* for i */
      } /* if multichar */

  /* deallocate pop */
   for (i=0; i<Pop_size; i++)
      {
      free(Pop[i]->genome);
      free(Kids[i]->genome);

      free(Pop[i]);
      free(Kids[i]);
      }  /* for i */

   free(Pop);
   free(Kids);
   free(Parents);

#ifdef DEBUG
   printf(" ---end close_pop()---\n");
#endif
   }  /* close_pop */

