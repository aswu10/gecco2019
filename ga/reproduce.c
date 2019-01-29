/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* reproduce.c
   07.24.98.AW	Created.  Routines having to do with reproduction,
		selecting parents.

   Routines:	reproduce()			07.24.98.AW
		calc_num_offspring()		07.24.98.AW
		select_parents()		07.24.98.AW
		select_p_with_replacement()	07.24.98.AW
		select_p_without_replacement()	07.24.98.AW
*/

#include <stdio.h>
#include "types.h"
#include "extern.h"
#include "reproduce.h"
#include "util.h"
#include "random.h"
#include "genops.h"
#include <float.h>

/********** reproduce **********/
/* parameters:
   called by:	ga_loop(), ga.c
   actions:	top level routine for reproduction.  Currently
		running a generational GA.
		- calculate expected number of offspring
		- select parents
		- crossover
		- mutate
*/
int reproduce()
   {
   int i;
   int g;
#ifdef DEBUG
   printf(" ---in reproduce()---\n");
#endif

  /* calculate expected number of offspring */
   calc_num_offspring();

#ifdef SMALLSTEP
   printf(" gen %d, after calculating expected number of offspring\n",
		Gen.index);
   print_population(Pop, 0, Pop_size-1);
   printf(" --- Press any key to continue ---\n");  fgetc(stdin);
#endif

  /* select parents */
   if (select_parents() == ERROR)  return ERROR;
    
#ifdef SMALLSTEP
   printf(" gen %d, after selecting parents.  Selected parents:\n",
                Gen.index);
   print_population(Parents, 0, Pop_size-1);
   printf(" --- Press any key to continue ---\n");  fgetc(stdin);
#endif

  /* crossover */
   if (crossover_pop() == ERROR)  return ERROR;
    
#ifdef SMALLSTEP
   printf(" gen %d, after crossover.  Offspring:\n", Gen.index);
   print_population(Kids, 0, Pop_size-1);
   printf(" --- Press any key to continue ---\n");  fgetc(stdin);
#endif

  /* mutate */
   mutate_pop();
    
#ifdef SMALLSTEP
   printf(" gen %d, after mutation.  Offspring (gen %d):\n",
		Gen.index, Gen.index+1);
   print_population(Kids, 0, Pop_size-1);
   printf(" --- Press any key to continue ---\n");  fgetc(stdin);
#endif

  /* if elitism is on, copy the best indv from Pop to the first slot of Kids */
   if (Elite)
      {
      //printf(" Save elite in gen %d\n", Gen.index);
      copy_indv(Pop[Gen.best_indv_index], Kids[0]);
      }

  /* if Random_immigrants > 0, generate random immigrants starting at the
     bottom of the population, e.g. starting at slot Pop_size -1 */
   if (Gen.index % RI_interval == 0)
      {
      if (Random_immigrants > 0)
         {
         for (i=Pop_size-1; i>=Pop_size-Random_immigrants; i--)
            {
            for (g=Kids[i]->length-1; g>=0; g--)
               {
               Kids[i]->floats_genome[g] = funiform(1);
               }
            decode(Kids[i]);
            }
         }
      }

  /* if Mass_extinction > 0, implement mass extinction starting at the
     bottom of the population, e.g. starting at slot Pop_size -1 */
  /* ME overwrites the RI individuals */
   if (Gen.index % ME_interval == 0)
      {
      if (Mass_extinction > 0)
         {
         for (i=Pop_size-1; i>=Pop_size-Mass_extinction; i--)
            {
            for (g=Kids[i]->length-1; g>=0; g--)
               {
               Kids[i]->floats_genome[g] = funiform(1);
               }
            decode(Kids[i]);
            }
         }
      }

#ifdef SMALLSTEP
   printf(" gen %d, after elitism and random immmigrants.  Offspring (gen %d):\n",
		Gen.index, Gen.index+1);
   print_population(Kids, 0, Pop_size-1);
   printf(" --- Press any key to continue ---\n");  fgetc(stdin);
#endif

#ifdef DEBUG
   printf(" ---end reproduce()---\n");
#endif
   return OK;
   }  /* reproduce */

/********** calc_num_offspring **********/
/* parameters:
   called by:	reproduce(), reproduce.c
   actions:	calculates expected number of offspring for each
		individual based on its fitness.  If Sigma_scaling_on
		is true, expected number of offspring is directly
		proportional to fitness.  If Sigma_scaling_on is false,
		the expected number of offspring for each individual
		will range from Sigma_scale_min to Sigma_scale_max.
*/
void calc_num_offspring()
   {
   int i;

#ifdef DEBUG
   printf(" ---in calc_num_offspring()---\n");
#endif

   if (!Sigma_scaling_on)
      {
      for (i=Pop_size-1; i>=0; i--)
        {
         Pop[i]->calc_num_offspring = 1/(Pop[i]->fitness/
			Gen.fitness_sum*Pop_size);
         // printf("fitness: %.2f\n", Pop[i]->fitness);
         // printf("offspring: %.2f\n", Pop[i]->calc_num_offspring);
        }
        // printf(" --- Press any key to continue ---\n");  fgetc(stdin);
      }  /* if */
   else
      {
      if (Gen.std_dev < 0.0001)
         {
         for (i=Pop_size-1; i>=0; i--)
            Pop[i]->calc_num_offspring = 1.0;
         }  /* if */
      else
         {
         for (i=Pop_size-1; i>=0; i--)
            {
            Pop[i]->calc_num_offspring = 
                1/(fmin(Sigma_scale_max,
                    fmax(Sigma_scale_min,
                         (1.0 + (Pop[i]->fitness - Gen.avg_fitness) /
                              (2.0 * Gen.std_dev) 
                         )
                        )));
            }
         }  /* else */
      }  /* else */

#ifdef DEBUG
   printf(" ---end calc_num_offspring()---\n");
#endif
   }  /* calc_num_offspring */

/********** select_parents **********/
/* parameters:
   called by:   reproduce(), reproduce.c
   actions:     selects parents for next generation from current population
		and saves them as pointers in the array parents.
		Can select with or without replacement.
		If entire population is eligible to be parents, then selecting
		without replacement essentially simply mixes up the individuals.
*/
int select_parents()
   {
#ifdef DEBUG
   printf(" ---in select_parents()---\n");
#endif

   if (!strcmp(Parent_selection, "proportional"))
      {
      if (Parent_replacement_on)
         select_p_with_replacement();
      else
         select_p_without_replacement();
      }  /* if prop */
   else if (!strcmp(Parent_selection, "tournament"))
      {
      tournament_selection();
      }  /* else if tourn */
   else
      {
      printf(" Error(select_parents): unknown selection method: %s\n",
		Parent_selection);
      return ERROR;
      }  /* else */

#ifdef DEBUG
   printf(" ---end select_parents()---\n");
#endif
   return OK;
   }  /* select_parents */

/********** select_p_with_replacement **********/
/* parameters:
   called by:   select_parents(), reproduce.c
   actions:     selects parents from current population with replacement.
		An individual my be selected as a parent more than once,
		meaning that it can have multiple offspring.
*/
void select_p_with_replacement()
   {
   int c, p;	/* ptrs to current pop and chosen parents */
   double random_num;
   double sum;
   int selected;

#ifdef DEBUG
   printf(" ---in select_p_with_replacement()---\n");
#endif

   for (p=0; p<Pop_size; p++)
      {
     /* generate a random number between 0 and Pop_size-1 */
      random_num = funiform((double)Pop_size);
      sum = 0.0;
      selected = -1;
      for (c=0; c<Pop_size; c++)
         {
         sum += Pop[c]->calc_num_offspring;
         if (sum >= random_num)
            {
            selected = c;
            break;
            }  /* if */
         }  /* for c */
      if (selected == -1)  selected = Pop_size-1;
      Parents[p] = Pop[selected];
      }  /* for p */

#ifdef DEBUG
   printf(" ---end select_p_with_replacement()---\n");
#endif
   }  /* select_p_with_replacement */

/********** select_p_without_replacement **********/
/* parameters:
   called by:   select_parents(), reproduce.c
   actions:     selects parents from current population with replacement.
		An individual can only be selected to be a parent once.
		When Pct_breeding and Pct_bred are both 1.0, this routine
		simple mixes up the individuals in the current population
		so that they reproduce with someone other than their 
		current neighbor.
*/
void select_p_without_replacement()
   {
   int c, p;    /* ptrs to current pop and chosen parents */
   int selected;
   double random_num;
   double sum;

#ifdef DEBUG
   printf(" ---in select_p_without_replacement()---\n");
#endif

   p = 0;
   while (p < Pop_size)
      {
     /* generate a random number between 0 and Pop_size-1 */
      random_num = funiform((double)Pop_size);
      sum = 0.0;
      selected = -1;
      for (c=0; c<Pop_size; c++)
         {
         sum += Pop[c]->calc_num_offspring;
         if (sum >= random_num)
            {
            selected = c;
            break;
            }  /* if */
         }  /* for c */
      if (selected == -1)  selected = Pop_size-1;
      if (!Pop[selected]->chosen)
         {
         Parents[p] = Pop[selected];
         Pop[selected]->chosen = 1;
         p++;
         }  /* if */
      }  /* while */

#ifdef DEBUG
   printf(" ---end select_p_without_replacement()---\n");
#endif
   }  /* select_p_without_replacement */

/********** tournament_selection **********/
/* parameters:
   called by:   select_parents(), reproduce.c
   actions:     selects parents from current population with tournament
		selection.  Tournament size specified in params.
		For each tournament, best fitness wins.
*/
void tournament_selection()
   {
   int p;	/* ptrs to parent pop to fill */
   int t;	/* tournament ptrs */
   double best_fitness;
   int best_index;
   int rand_num;
 
#ifdef DEBUG
   printf(" ---in tournament_selection()---\n");
#endif

   for (p=Pop_size-1; p>=0; p--)
      {
      best_fitness = DBL_MAX;
      best_index = -1;

      for (t=Tournament_size-1; t>=0; t--)
         {
         rand_num = uniform(Pop_size);
         if (Pop[rand_num]->fitness < best_fitness)
            {
            best_index = rand_num;
            best_fitness = Pop[rand_num]->fitness;
            }  /* if found better fitness */
         }  /* for t */   
      Parents[p] = Pop[best_index];
      }  /* for p */

#ifdef DEBUG
   printf(" ---end tournament_selection()---\n");
#endif
   }  /* tournament_selection */

