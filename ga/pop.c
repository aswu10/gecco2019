/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* pop.c
   07.03.98.AW	Created.
		Routines that have to do with populations.
		Such as setting up the data structure and initializing it.

   Routines:	init_pop()			07.03.98.AW
		pop_allocate()			07.19.17 RN
		pop_fill()			07.19.17 RN
		pop_individuals_1()		07.07.98.AW
		pop_individuals_random()	07.07.98.AW
		pop_eval()			07.08.98.AW
   		close_pop()			07.03.98.AW
*/

#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include "types.h"
#include "extern.h"
#include "pop.h"
#include "params.h"
#include "random.h"



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
         Pop[i]->length = Max_gen_len;
         Kids[i]->length = Max_gen_len;
         }  /* if */
      else if (Variable_gen_len > 1 && Variable_gen_len <= Max_gen_len)
         {
         Pop[i]->length = Variable_gen_len;
         Kids[i]->length = Variable_gen_len;
         }  /* else if */
      else if (Variable_gen_len == 1)
         {
         Pop[i]->length = Min_gen_len + uniform(Max_gen_len - Min_gen_len +1);
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
      
      if (Init_pop == 2)  
         {
         Pop[i]->floats_genome = (double *)malloc(Max_gen_len * sizeof(double));
        Kids[i]->floats_genome = (double *)malloc(Max_gen_len * sizeof(double));
         }
      }  /* for i */

#ifdef DEBUG
   printf(" ---end pop_allocate()---\n");
#endif
   return OK;
   }  /* pop_allocate */

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

   switch(Init_pop)
      {
      case -1:	break;
	  case 1:    pop_individuals_random();     
         break;
      case 2:	pop_individuals_rkeys();
	     break;
      default:	printf(" Error(pop_fill): Init_pop value invalid");
		 printf(" for int base: %d\n", Init_pop);
		 return ERROR;
      }  /* switch */

#ifdef DEBUG
   printf(" ---end pop_fill()---\n");
#endif
   return OK;
   }  /* pop_fill */

/********** pop_individuals_random **********/
/* parameters:
   called by:   pop_fill(), pop.c
   actions:     fills population with individuals of random 1's and 0's.
   07.23.17 RN: Modified for int & random keys
*/
void pop_individuals_random()
   {
   int i, j;

#ifdef DEBUG
   printf(" ---in pop_individuals_random()---\n");
#endif
   for (i=Pop_size-1; i>=0; i--)
      {
      for (j = 0; j < Pop[i]->length; j++)   
         {
         int random = uniform(Max_gen_len);
         if (j > 0) 
            {
            while (int_in_list(random, Pop[i]->genome, j) == 1) 
               {
               random = uniform(Max_gen_len);
               }
            }
         Pop[i]->genome[j] = random;
         }  /* for j */
      }  /* for i */
      
#ifdef DEBUG
   printf(" ---end pop_individuals_random()---\n");
#endif
   }  /* pop_individuals_random */

/********** pop_individuals_rkeys **********/
/* parameters:
   called by:   pop_fill(), pop.c
   actions:     fills population with individuals of random keys between 0,1
   07.23.17 RN: Added for Random Keys
*/
void pop_individuals_rkeys()
   {
   int i, j;

#ifdef DEBUG
   printf(" ---in pop_individuals_rkeys()---\n");
#endif

   for (i=Pop_size-1; i>=0; i--)
      {
      for (j=Pop[i]->length-1; j>=0; j--)
         {
         Pop[i]->floats_genome[j] = funiform(1);
		 decode(Pop[i]);
         }  /* for j */
      }  /* for i */

#ifdef DEBUG
   printf(" ---end pop_individuals_rkeys()---\n");
#endif
   }  /* pop_individuals_rkeys */

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
   07.19.17 RN: Commented out for integer conversion
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

  /* deallocate pop */
   for (i=0; i<Pop_size; i++)
      {
      free(Pop[i]->genome);
      free(Kids[i]->genome);

      if (Init_pop == 2) // 07.23.17 RN: for Random Keys
         {
         free(Pop[i]->floats_genome);
         free(Kids[i]->floats_genome);      
         }
      
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

