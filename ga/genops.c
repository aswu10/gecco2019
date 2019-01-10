/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* genops.c
   07.24.98.AW	Created.  Routines having to do with genetic operators,
		such as crossover and mutation.

   Routines:	crossover_pop()		07.27.98.AW
		mutate_pop()		07.27.98.AW
		onept_crossover()	07.27.98.AW
		twopt_crossover()	07.27.98.AW
		uniform_crossover()	07.27.98.AW
		mutate()		07.19.17 RN
		poisson()		07.27.98.AW
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "types.h"
#include "extern.h"
#include "genops.h"
#include "random.h"
#include "params.h"
#include "util.h"

#define DEBUG_MUT	/* debug mutation */
#define HX		/* debug homologous crossover */
#undef DEBUG_MUT
#undef HX

/********** crossover_pop **********/
/* parameters:
   called by:	reproduce(), reproduce.c
   actions:	Perform xover on each pair of breeding individuals
		from selected parents.
*/
int crossover_pop()
   {
   int i, j;
   int num_pairs;

#ifdef DEBUG
   printf(" ---in crossover_pop()---\n");
#endif

   num_pairs = Num_breeding/2;

   if (!Variable_gen_len)
      {
      if (!strcmp(Xover_type, "one-point"))
         for (i=0, j=0; i<num_pairs; i++, j+=2)
            onept_crossover(Parents[j], Parents[j+1], Kids[j], Kids[j+1]);
      else if (!strcmp(Xover_type, "two-point"))
         for (i=0, j=0; i<num_pairs; i++, j+=2)
            twopt_crossover(Parents[j], Parents[j+1], Kids[j], Kids[j+1]);
      else if (!strcmp(Xover_type, "uniform")) {
         for (i=0, j=0; i<num_pairs; i++, j+=2)
            uniform_crossover(Parents[j], Parents[j+1], Kids[j], Kids[j+1]); }
      else if (!strcmp(Xover_type, "position"))
         for (i=0, j=0; i<num_pairs; i++, j+=2)
            position_based_crossover(Parents[j], Parents[j+1], Kids[j], Kids[j+1]);
      else if (!strcmp(Xover_type, "edge-recombination"))
         for (i=0, j=0; i<num_pairs; i++, j+=2) 
            edge_recombination_crossover(Parents[j], Parents[j+1], Kids[j], Kids[j+1]);
      else if (!strcmp(Xover_type, "partially-mapped"))
         for (i=0, j=0; i<num_pairs; i++, j+=2)
             
            partially_mapped_crossover(Parents[j], Parents[j+1], Kids[j], Kids[j+1]);
      else
         {
         printf(" Error(crossover_pop): invalid crossover for fixed len: %s\n",
		Xover_type);
         return ERROR;
         }  /* else */
      }  /* if not variable length */
   else  /* if variable length */
      {
      if (!strcmp(Xover_type, "homologous"))
         for (i=0, j=0; i<num_pairs; i++, j+=2)
            homologous_crossover(Parents[j], Parents[j+1], Kids[j], Kids[j+1]);
      else
         {
         printf(" Error(crossover_pop): invalid crossover for var len: %s\n",
		Xover_type);
         return ERROR;
         }  /* else */
      }  /* else */

#ifdef DEBUG
   printf(" ---end crossover_pop()---\n");
#endif
   return OK;
   }  /* crossover_pop */

/********** mutate_pop **********/
/* parameters:
   called by:	reproduce(), reproduce.c
   actions:
*/
void mutate_pop()
   {
   int i;
 
#ifdef DEBUG
   printf(" ---in mutate_pop()---\n");
#endif
 
   for (i=0; i<Num_breeding; i++)
      {
      if (!strcmp(Mut_type, "random")) random_mutate(Kids[i]);   
      else if (!strcmp(Mut_type, "gaussian")) gaussian_mutate(Kids[i]);
      else if (!strcmp(Mut_type, "displacement"))
         {
         int start = uniform(Kids[i]->length - 1);
         int length = uniform(Kids[i]->length - start);
         displacement_mutate(Kids[i], start, length, 0, 0);        
         } /* if */          
      else if (!strcmp(Mut_type, "insertion"))
         {
         // insertion is essential a special case of displacement
         // where the length == 1
         int location = uniform(Kids[i]->length - 1);
         displacement_mutate(Kids[i], location, 1, 0, 0);        
         } /* else if */       
      else if (!strcmp(Mut_type, "insertion_per_gene"))
         {
         // insertion is essential a special case of displacement
         // where the length == 1
         int location = uniform(Kids[i]->length - 1);
         displacement_mutate(Kids[i], location, 1, 0, 1);        
         } /* else if */  
      else if (!strcmp(Mut_type, "inversion"))
         {
         // insertion is essential a special case of displacement
         // where the the values are inverted before re-insertion
         int start = uniform(Kids[i]->length - 1);
         int length = uniform(Kids[i]->length - start);
         displacement_mutate(Kids[i], start, length, 1, 0);        
         } /* else if */       
      }  /* for i */
 
#ifdef DEBUG
   printf(" ---end mutate_pop()---\n");
#endif
   }  /* mutate_pop */

/********** onept_crossover **********/
/* parameters:
   called by:	crossover_pop(), genop.c
   actions:	Given ptrs to two parents and two kids.
		Randomly generate a number between 0 and 1.
		If value is less than Xover_rate,
		perform 1pt crossover on parents to create new kids.
		Else simply copy parents to kids structures.
		Includes setting relevant parent data in the kid structure.
		Randomly generates a crossover point within length
		of parent1 (assumes fixed length individuals) and 
		stores in cross_point.  Bits 0 to cross_point-1
		go to kid1, bits cross_point to length-1 go to kid2.
   restricxns:	for fixed length individuals only.
*/
void onept_crossover(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
			INDIVIDUAL *kid1, INDIVIDUAL *kid2)
   {
   double cross_prob;
   int cross_point;
   int i;
 
#ifdef DEBUG
   printf(" ---in onept_crossover()---\n");
#endif

   cross_prob = knuth_random();

   if (cross_prob < Xover_rate)
      {
     /* yes, crossover, do crossover and copy data */
      cross_point = uniform(parent1->length);
      for (i=0; i<cross_point; i++)
         {
         kid1->floats_genome[i] = parent1->floats_genome[i];
         kid2->floats_genome[i] = parent2->floats_genome[i];
         }  /* for i */
      for (i=cross_point; i<parent1->length; i++)
         {
         kid1->floats_genome[i] = parent2->floats_genome[i];
         kid2->floats_genome[i] = parent1->floats_genome[i];
         }  /* for i */
     /* set parent info, parent 1 contributes first portion */
      kid1->parent1_index = parent1->index;
      kid1->parent2_index = parent2->index;
      kid2->parent1_index = parent2->index;
      kid2->parent2_index = parent1->index;

      // decode(kid1);
      // decode(kid2);
      
      }  /* if */
   else
      {
     /* no, don't crossover, just copy parents to kid structures */
      for (i=0; i<parent1->length; i++)
         {
         kid1->floats_genome[i] = parent1->floats_genome[i];
         kid2->floats_genome[i] = parent2->floats_genome[i];
         
         kid1->genome[i] = parent1->genome[i];
         kid2->genome[i] = parent2->genome[i];
         }  /* for i */
     /* set parent info, set parent2 to -1 if clones (no xover) */
      kid1->parent1_index = parent1->index;
      kid1->parent2_index = -1;
      kid2->parent1_index = parent2->index;
      kid2->parent2_index = -1;

      }  /* else */

#ifdef SMALLERSTEP
   if (cross_prob < Xover_rate)
      printf(" ********** one point crossover at %d\n", cross_point);
   else
      printf(" ********** no crossover\n");
   printf(" %3d ", parent1->index);
   print_genome(parent1, 0);
   printf("  ");
   printf("%3d ", kid1->index);
   print_genome(kid1, 1);
   printf("     ");
   if (cross_prob < Xover_rate)
      {
      for (i=0; i<cross_point; i++)  printf(" ");
      printf("|");
      for (i=cross_point+1; i<parent1->length; i++)  printf(" ");
      printf("->    ");
      for (i=0; i<cross_point; i++)  printf(" ");
      printf("|\n");
      }  /* if */
   else
      {
      for (i=0; i<parent1->length; i++)  printf(" ");
      printf("->\n");
      }  /* if */
   printf(" %3d ", parent2->index);
   print_genome(parent2, 0);
   printf("  ");
   printf("%3d ", kid2->index);
   print_genome(kid2, 1);
   printf(" --- Press any key to continue ---\n");  fgetc(stdin);
#endif  /* SMALLERSTEP */

#ifdef DEBUG
   printf(" ---end onept_crossover()---\n");
#endif
   }  /* onept_crossover */

/********** twopt_crossover **********/
/* parameters:
   called by:	crossover_pop(), genop.c
   actions:	Given ptrs to two parents and two kids,
                perform 2pt crossover on parents to create
                new kids, if cross_prob is less than Xover_rate
		(see description in onept_crossover).
                Includes copying all of the parent data to the
                kid structure.
   restricxns:	for fixed length individuals only.
*/
void twopt_crossover(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
			INDIVIDUAL *kid1, INDIVIDUAL *kid2)
   {
   double cross_prob;
   int cross_point;
   int cross_point2;
   int i;
 
#ifdef DEBUG
   printf(" ---in twopt_crossover()---\n");
#endif
 
   cross_prob = knuth_random();
 
   if (cross_prob < Xover_rate)
      {
     /* yes, crossover, do crossover and copy data */
      cross_point = uniform(parent1->length);
      cross_point2 = uniform(parent1->length);
      if (cross_point > cross_point2)
         {
         i = cross_point;
         cross_point = cross_point2;
         cross_point2 = i;
         }  /* if */
      for (i=0; i<cross_point; i++)
         {
         kid1->floats_genome[i] = parent1->floats_genome[i];
         kid2->floats_genome[i] = parent2->floats_genome[i];
         }  /* for i */
      for (i=cross_point; i<cross_point2; i++)
         {
         kid1->floats_genome[i] = parent2->floats_genome[i];
         kid2->floats_genome[i] = parent1->floats_genome[i];
         }  /* for i */
      for (i=cross_point2; i<parent1->length; i++)
         {
         kid1->floats_genome[i] = parent1->floats_genome[i];
         kid2->floats_genome[i] = parent2->floats_genome[i];
         }  /* for i */
     /* set parent info, parent 1 contributes first portion */
      kid1->parent1_index = parent1->index;
      kid1->parent2_index = parent2->index;
      kid2->parent1_index = parent2->index;
      kid2->parent2_index = parent1->index;

	  // decode(kid1);
      // decode(kid2);
	  
      }  /* if */
   else
      {
     /* no, don't crossover, just copy parents to kid structures */
      for (i=0; i<parent1->length; i++)
         {
         kid1->floats_genome[i] = parent1->floats_genome[i];
         kid2->floats_genome[i] = parent2->floats_genome[i];
         
         kid1->genome[i] = parent1->genome[i];
         kid2->genome[i] = parent2->genome[i];
         }  /* for i */
     /* set parent info, set parent2 to -1 if clones (no xover) */
      kid1->parent1_index = parent1->index;
      kid1->parent2_index = -1;
      kid2->parent1_index = parent2->index;
      kid2->parent2_index = -1;

      }  /* else */

#ifdef SMALLERSTEP
   if (cross_prob < Xover_rate)
      printf(" ********** two point crossover at %d and %d\n",
		cross_point, cross_point2);
   else
      printf(" ********** no crossover\n");
   printf(" %3d ", parent1->index);
   print_genome(parent1->floats_genome, 0);
   printf("  ");
   printf("%3d ", kid1->index);
   print_genome(kid1->floats_genome, 1);
   printf("     ");
   if (cross_prob < Xover_rate)
      {
      for (i=0; i<cross_point; i++)  printf(" ");
      printf("|");
      for (i=cross_point+1; i<cross_point2; i++)  printf(" ");
      printf("|");
      for (i=cross_point2+1; i<parent1->length; i++)  printf(" ");
      printf("->    ");
      for (i=0; i<cross_point; i++)  printf(" ");
      printf("|");
      for (i=cross_point+1; i<cross_point2; i++)  printf(" ");
      printf("|\n");
      }  /* if */
   else
      {
      for (i=0; i<parent1->length; i++)  printf(" ");
      printf("->\n");
      }  /* if */
   printf(" %3d ", parent2->index);
   print_genome(parent2, 0);
   printf("  ");
   printf("%3d ", kid2->index);
   print_genome(kid2, 1);
   printf(" --- Press any key to continue ---\n");  fgetc(stdin);
#endif  /* SMALLERSTEP */

#ifdef DEBUG
   printf(" ---end twopt_crossover()---\n");
#endif
   }  /* twopt_crossover */

/********** uniform_crossover **********/
/* parameters:
   called by:   crossover_pop(), genop.c
   actions:     Given ptrs to two parents and two kids,
                perform uniform crossover on parents to create
                new kids.
                Includes copying all of the parent data to the
                kid structure.

		Child 1 gets from parent 1 if random number less than
		Uniform_x; child one gets from parent 2 if random number
		greater than Uniform_x.

   restricxns:  for fixed length individuals only.
*/
void uniform_crossover(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
                        INDIVIDUAL *kid1, INDIVIDUAL *kid2)
   {
   int i;
   double cross_prob;   /* xover_rate */
   double x_prob;       /* uniform_x rate */
   int inher;           /* stores previous inheritance pattern */
                        /* p1 -> o1 or p1 -> o2, for Trace purposes */
 
#ifdef DEBUG
   printf(" ---in uniform_crossover()---\n");
#endif
 
   cross_prob = knuth_random();
 
   if (cross_prob < Xover_rate)
      {
     /* yes, do crossover */
 
     /* copye nucleotides from parents to kids */
      for (i=parent1->length-1; i>=0; i--)
         {
         x_prob = knuth_random();
         if (x_prob < Uniform_x)
            {
           /* inher = 1; parent1 to kid1, parent2 to kid2 */
            kid1->floats_genome[i] = parent1->floats_genome[i];
            kid2->floats_genome[i] = parent2->floats_genome[i];
 
            }  /* if */
         else
            {
           /* inher = -1; parent1 to kid2, parent2 to kid1 */
            kid1->floats_genome[i] = parent2->floats_genome[i];
            kid2->floats_genome[i] = parent1->floats_genome[i];
 
            }  /* else */
         }  /* for i */
 
     /* set parent info, parent 1 contributes first portion */
      kid1->parent1_index = parent1->index;
      kid1->parent2_index = parent2->index;
      kid2->parent1_index = parent2->index;
      kid2->parent2_index = parent1->index;
 
	  // decode(kid1);
      // decode(kid2);
 
      }  /* if yes crossover */
   else
      {
     /* no, don't crossover, just copy parents to kid structures */
      for (i=0; i<parent1->length; i++)
         {
         kid1->floats_genome[i] = parent1->floats_genome[i];
         kid2->floats_genome[i] = parent2->floats_genome[i];
         
         kid1->genome[i] = parent1->genome[i];
         kid2->genome[i] = parent2->genome[i];
         }  /* for i */
     /* set parent info, set parent2 to -1 if clones (no xover) */
      kid1->parent1_index = parent1->index;
      kid1->parent2_index = -1;
      kid2->parent1_index = parent2->index;
      kid2->parent2_index = -1;
 
      }  /* else don't crossover */

#ifdef DEBUG
   printf(" ---end uniform_crossover()---\n");
#endif
   }  /* uniform_crossover */
   
   /********** position_based_crossover **********/
/* parameters:
   called by:   crossover_pop(), genop.c
   actions:     

		

   restricxns:  
*/
void position_based_crossover(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
                        INDIVIDUAL *kid1, INDIVIDUAL *kid2)
   {
   int i;
   double cross_prob;   /* xover_rate */
   double x_prob;       /* uniform_x rate */

#ifdef DEBUG
   printf(" ---in position_based_crossover()---\n");
#endif

   for (i=0; i<parent1->length; i++)
      {
      kid1->genome[i] = parent1->genome[i];
      kid2->genome[i] = parent2->genome[i];
      }  /* for i */
 
   cross_prob = knuth_random();
 
   if (cross_prob < Xover_rate)
      {
      
      for (i=parent1->length-1; i>=0; i--)
         {
         x_prob = knuth_random();
         /* this position selected for crossover */
         if (x_prob < Uniform_x)  
            {
            int value1 = kid1->genome[i];
            int parent2_value = parent2->genome[i];
            int idx_parent2_value = idx_in_list(parent2_value, kid1->genome, parent2->length);
            kid1->genome[i] = parent2_value;
            kid1->genome[idx_parent2_value] = value1;  
            
            int value2 = kid2->genome[i];
            int parent1_value = parent1->genome[i];
            int idx_parent1_value = idx_in_list(parent1_value, kid2->genome, parent1->length);
            kid2->genome[i] = parent1_value;
            kid2->genome[idx_parent1_value] = value2;
            }
         }  /* for i */
 
     /* set parent info, parent 1 contributes first portion */
      kid1->parent1_index = parent1->index;
      kid1->parent2_index = parent2->index;
      kid2->parent1_index = parent2->index;
      kid2->parent2_index = parent1->index;
      }  /* if yes crossover */
   else
      {
     /* no, don't crossover, just copy parents to kid structures */
     /* set parent info, set parent2 to -1 if clones (no xover) */
      kid1->parent1_index = parent1->index;
      kid1->parent2_index = -1;
      kid2->parent1_index = parent2->index;
      kid2->parent2_index = -1;
 
      }  /* else don't crossover */
      
#ifdef DEBUG
   printf(" ---end position_based_crossover()---\n");
#endif
   }  /* position_based_crossover */
   
void partially_mapped_crossover(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
                        INDIVIDUAL *kid1, INDIVIDUAL *kid2)
   {
   int i;
   // int j;
   double cross_prob;   /* xover_rate */

#ifdef DEBUG
   printf(" ---in partially_mapped_crossover()---\n");
#endif

   // for (i=0; i<parent1->length; i++)
      // {
      // kid1->genome[i] = parent1->genome[i];
      // kid2->genome[i] = parent2->genome[i];
      // }  /* for i */
 
   cross_prob = knuth_random();
 
   if (cross_prob < Xover_rate)
      {
      int start = uniform(Max_gen_len);
      int end = uniform(Max_gen_len - start) + start + 1;
  
      // copy over mapped section from parents to children
      for (i = start; i < end; i++)
         {
         kid1->genome[i] = parent2->genome[i];
         kid2->genome[i] = parent1->genome[i];         
         }
      
      // for rest of chromo, copy over values directly if not in the already mapped section
      // if already in mapped section, get what that value replaced in the parent's chromo
      for (i = 0; i < start; i++)
         {
         kid1->genome[i] = parent1->genome[i];
         while (int_in_list(kid1->genome[i], &kid1->genome[start], end - start) == 1)
            {
            int idx = start + idx_in_list(kid1->genome[i], &parent2->genome[start], end - start);
            kid1->genome[i] = parent1->genome[idx];
            }
         
         kid2->genome[i] = parent2->genome[i];        
         while (int_in_list(kid2->genome[i], &kid2->genome[start], end - start) == 1) 
            {
            int idx = start + idx_in_list(kid2->genome[i], &parent1->genome[start], end - start);
            kid2->genome[i] = parent2->genome[idx];
            }     
         }
        
      for (i = end; i < parent1->length; i++)
         {
         kid1->genome[i] = parent1->genome[i];
         while (int_in_list(kid1->genome[i], &kid1->genome[start], end - start) == 1)
            {
            int idx = start + idx_in_list(kid1->genome[i], &parent2->genome[start], end - start);
            kid1->genome[i] = parent1->genome[idx];
            }
         
         kid2->genome[i] = parent2->genome[i];        
         while (int_in_list(kid2->genome[i], &kid2->genome[start], end - start) == 1) 
            {
            int idx = start + idx_in_list(kid2->genome[i], &parent1->genome[start], end - start);
            kid2->genome[i] = parent2->genome[idx];
            }     
         }
 
     // /* set parent info, parent 1 contributes first portion */
      kid1->parent1_index = parent1->index;
      kid1->parent2_index = parent2->index;
      kid2->parent1_index = parent2->index;
      kid2->parent2_index = parent1->index;
      }  /* if yes crossover */
   else
      {
     /* no, don't crossover, just copy parents to kid structures */
     for (i=0; i<parent1->length; i++)
         {
         kid1->genome[i] = parent1->genome[i];
         kid2->genome[i] = parent2->genome[i];
         }  /* for i */
     /* set parent info, set parent2 to -1 if clones (no xover) */
      kid1->parent1_index = parent1->index;
      kid1->parent2_index = -1;
      kid2->parent1_index = parent2->index;
      kid2->parent2_index = -1;
      }  /* else don't crossover */
      
#ifdef DEBUG
   printf(" ---end partially_mapped_crossover()---\n");
#endif
   }  /* partially_mapped_crossover */
   
/********** edge_recombination_crossover **********/
/* parameters:
   called by:   crossover_pop(), genop.c
   actions:     

		

   restricxns: assumes constant length chromosomes
*/ 
void edge_recombination_crossover(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
                        INDIVIDUAL *kid1, INDIVIDUAL *kid2)
   {
    
#ifdef DEBUG
   printf(" ---in edge_recombination_crossover()---\n");
#endif
    
    int i;
    int j;
    double cross_prob;   /* xover_rate */
    EDGE_MAP *maps;
    EDGE_MAP *maps_copy;
    
    cross_prob = knuth_random();
 
    if (cross_prob < Xover_rate) 
       {
       
       // construct edge map
       maps = (EDGE_MAP *)malloc(parent1->length * sizeof(EDGE_MAP));
       maps_copy = (EDGE_MAP *)malloc(parent1-> length * sizeof(EDGE_MAP));
       for (i=0; i<parent1->length; i++)
          {
          maps[i].edges = (int *)malloc(4 * sizeof(int));
          maps[i].length = 0;
          maps_copy[i].edges = (int *)malloc(4 * sizeof(int));
          maps_copy[i].length = 0;
          }
          
       get_edges_from_parent(maps, parent1);
       get_edges_from_parent(maps, parent2);
          
       for (j=0; j<parent1->length; j++)
          {
          maps_copy[j].length = maps[j].length;
          int k;
          for (k=0; k<maps_copy[j].length; k++) 
             {
             maps_copy[j].edges[k] = maps[j].edges[k]; 
             }
          }
          
       // printf("parent1: \n");
       // for (j=0; j<parent1->length; j++) {
           // printf("%d ", parent1->genome[j]);
       // }
       // printf("\n");
       // printf("parent2: \n");
       // for (j=0; j<parent2->length; j++) {
           // printf("%d ", parent2->genome[j]);
       // }
       // printf("\n");
       // printf("edge-map: \n");
       // int k;
       // for (j=0; j<parent1->length; j++) {
           // printf("%d:  ",  j);
           // for (k=0; k < maps[j].length; k++) {
              // printf("%d ", maps[j].edges[k]);
           // }
           // printf("\n");
       // }
       // printf("edge-map copy: \n");
       // for (j=0; j<parent1->length; j++) {
           // printf("%d:  ",  j);
           // for (k=0; k < maps_copy[j].length; k++) {
               // printf("%d ", maps_copy[j].edges[k]);
           // }
           // printf("\n");
       // }
       // printf(" --- Press any key to continue ---\n");  fgetc(stdin);

    
       edge_recombination_to_child(kid1, maps, parent1, parent2);
       edge_recombination_to_child(kid2, maps_copy, parent1, parent2);
       
       // printf("kid1: \n");
       // for (j=0; j<kid1->length; j++) {
           // printf("%d ", kid1->genome[j]);
       // }
       // printf("\n");
       // printf("kid2: \n");
       // for (j=0; j<kid2->length; j++) {
           // printf("%d ", kid2->genome[j]);
       // }
       // printf("\n");
              // printf("edge-map: \n");
       // for (j=0; j<parent1->length; j++) {
           // printf("%d:  ",  j);
           // for (k=0; k < maps[j].length; k++) {
               // printf("%d ", maps[j].edges[k]);
           // }
           // printf("\n");
       // }
       // printf("edge-map copy: \n");
       // for (j=0; j<parent1->length; j++) {
           // printf("%d:  ",  j);
           // for (k=0; k < maps_copy[j].length; k++) {
               // printf("%d ", maps_copy[j].edges[k]);
           // }
           // printf("\n");
       // }
       // printf(" --- Press any key to continue ---\n");  fgetc(stdin);
       
       for (i=0; i<parent1->length; i++)
          {
          free(maps[i].edges);
          free(maps_copy[i].edges);
          }
       free(maps);
       free(maps_copy);
       
       // /* set parent info, parent 1 contributes first portion */
       kid1->parent1_index = parent1->index;
       kid1->parent2_index = parent2->index;
       kid2->parent1_index = parent2->index;
       kid2->parent2_index = parent1->index;
       }  /* if yes crossover */
    else
       {
       /* no, don't crossover, just copy parents to kid structures */
       for (i=0; i<parent1->length; i++)
          {
          kid1->genome[i] = parent1->genome[i];
          kid2->genome[i] = parent2->genome[i];
          }  /* for i */
       /* set parent info, set parent2 to -1 if clones (no xover) */
       kid1->parent1_index = parent1->index;
       kid1->parent2_index = -1;
       kid2->parent1_index = parent2->index;
       kid2->parent2_index = -1;
       }  /* else don't crossover */
       


    
#ifdef DEBUG
   printf(" ---end edge_recombination_crossover()---\n");
#endif      
   } /* edge_recombination_crossover */
   
void get_edges_from_parent(EDGE_MAP *maps, INDIVIDUAL *parent)
   {
    int i;
    int j;
    for (i=0; i<parent->length; i++)
       {
       
       int first_value = parent->genome[0];
       int last_value;
       
       for (j=0; j<parent->length; j++) 
          {
          int key = parent->genome[j];
          
          if (j > 0) 
             {         
             int value = parent->genome[j-1];
             if (!int_in_list(value, maps[key].edges, maps[key].length))
                {
                maps[key].edges[maps[key].length] = value;
                maps[key].length++;  
                } /* if */
             } /* if */
          
          if (j < parent->length-1) 
             {
             int value = parent->genome[j+1];
             if (!int_in_list(value, maps[key].edges, maps[key].length)) 
                {
                maps[key].edges[maps[key].length] = value;
                maps[key].length++;           
                } /* if */
             } /* if */
             
          if (j == parent->length-1 && !int_in_list(first_value, maps[key].edges, maps[key].length))
             {
             maps[key].edges[maps[key].length] = first_value;
             maps[key].length++; 
             
             last_value = key;
             maps[first_value].edges[maps[first_value].length] = last_value;
             maps[first_value].length++; 
             
             } /* if */
          } /* for */

       } /* for */    
   }
   
void edge_recombination_to_child(INDIVIDUAL *child, EDGE_MAP *maps, INDIVIDUAL *parent1, INDIVIDUAL *parent2)
   {
    int i;
    int j;
    
    // create list of unvisited cities
    int unvisited[Max_gen_len];
    int num_unvisited = Max_gen_len;
    for (i=0; i < Max_gen_len; i++)
       {
       unvisited[i] = i;
       }
        
    // 1
    // choose initial from random parent
    int current;
    double random_double = funiform(1);
    if (random_double < 0.5) 
       {
       current = parent1->genome[0];
       }
    else 
       {
       current = parent2->genome[0];    
       }
    child->genome[0] = current;
       
    // repeat until full tour complete 
    for (i=1; i<Max_gen_len; i++)
       {
       // 2
       // remove current from edge map
       remove_from_maps(current, maps, Max_gen_len, unvisited, num_unvisited);
       num_unvisited--;
       
       // 3
       // select next 
       if (maps[current].length > 0)
          {
          // 4
          // next current is neighbor w/ fewest edges
          int new_current;
          int num_edges = REALLY_BIG_NUMBER;
          int tied[4];
          int num_ties = 0;
          for (j=0; j<maps[current].length; j++)
             {
             int value = maps[current].edges[j];
             if (maps[value].length < num_edges)
                {
                num_edges = maps[value].length;
                new_current = value;
                } /* if */
             else if (maps[value].length == num_edges) // record ties
                {
                if (num_ties == 0)
                   {
                   tied[num_ties] = new_current;
                   num_ties++;
                   } /* if */
                tied[num_ties] = value;
                num_ties++;
                } /* else if */
             } /* for */ 
             if (num_ties > 0) //resolve ties
                {
                int random = uniform(num_ties);
                new_current = tied[random];
                } /* if */
             current = new_current;
          } /* if */
       
       // 5       
       // if no neighbors, choose randomly from edgeless    
       else
          {
          int random = uniform(num_unvisited);
          current = unvisited[random];
          } /* else */

       child->genome[i] = current;
       } /* for */
   }  

void remove_from_maps(int value, EDGE_MAP *maps, int length, int *unvisited, int num_unvisited)
   {
   int i;
   int j;
   
   // remove from right side of edge maps
   for (i=0; i<length; i++)
      {
      int item_found = 0;
      for (j=0; j<maps[i].length; j++)
         {
         if (maps[i].edges[j] == value)
            {
            item_found = 1;   
            } /* if */
         if (item_found && j < maps[i].length-1) 
            {
            maps[i].edges[j] = maps[i].edges[j+1];
            } /* if */
         } /* for */
      if (item_found)
         {
         maps[i].length--;
         } /* if */
      } /* for */
      
   // remove from list of unvisited items
   int item_found = 0;
   for (i=0; i<num_unvisited; i++)
      {
      if (unvisited[i] == value)
         {
         item_found = 1; 
         } /* if */
      if (item_found && i < num_unvisited-1) 
         {
         unvisited[i] = unvisited[i+1];
         } /* if */
      } /* for */
   } /* remove_from_maps */
   
/********** switch_crossover **********/
/* parameters:
   called by:   crossover_pop(), genop.c
   actions:	Given ptrs to two parents and two kids,
                perform uniform crossover on parents to create
                new kids.
                Includes copying all of the parent data to the
                kid structure.

		Each time the randomly generated number is below Uniform_x,
		switch parents.

   restricxns:	for fixed length individuals only.
*/

/********** random_mutate **********/
/* parameters:
   called by:	mutate_pop(), genop.c
   actions:	given an individual, mutate it.
		Use a poisson distribution to decide how many mutations
		should be performed on this individual.  Then chooses
		mutation locations randomly, with replacement.
   07.19.17 RN: Converted to int, mutates to a random int 0 - 9, 
   is possible to change to current int, effenctively making no change
   07.23.17 RN: Added random keys version
   08.02.17 RN: changed name to random_mutate 
*/
void random_mutate(INDIVIDUAL *indv)
   {
   int i;
   int num_mutations;
   int random_num;

#ifdef DEBUG
   printf(" ---in random_mutate()---\n");
#endif
   
   num_mutations = poisson(((double)Mut_rate * (double)indv->length));

   for (i=num_mutations-1; i>=0; i--)
      {
      random_num = uniform(indv->length);
      if (Init_pop != 2)
         {
         indv->genome[random_num] = uniform(Max_gen_len); 
         } 
      else // random keys
         {
         indv->floats_genome[random_num] = funiform(1); 
         // decode(indv);
         }
      }  /* for i */

#ifdef DEBUG
   printf(" ---end random_mutate()---\n");
#endif
   }  /* random_mutate */
   
void gaussian_mutate(INDIVIDUAL *indv)
   {
   int i;
   int num_mutations;
   int random_num;

#ifdef DEBUG
   printf(" ---in gaussian_mutate()---\n");
#endif

   num_mutations = poisson(((double)Mut_rate * (double)indv->length));

   for (i=num_mutations-1; i>=0; i--)
      {
      random_num = uniform(indv->length);
      if (Init_pop != 2)
         {
         indv->genome[random_num] = uniform(Max_gen_len); 
         } 
      else // random keys
         {
          // printf("origingal value: %.2lf", indv->floats_genome[random_num]);
         double value = gaussian(indv->floats_genome[random_num], Gaussian_sd);
         if (value > 1.0) 
            {
            value = 1.0;
            }
         else if (value < 0.0) 
            {
            value = 0.0;
            }
         
         indv->floats_genome[random_num] = value; 
         }
      }  /* for i */

#ifdef DEBUG
   printf(" ---end gaussian_mutate()---\n");
#endif
   }  /* gaussian_mutate */
   
void displacement_mutate(INDIVIDUAL *indv, int start, int length, int invert, int per_gene)
   {
#ifdef DEBUG
   printf(" ---in displacement_mutate()---\n");
#endif

   double mut_prob;
   if (per_gene == 1)
      {
      mut_prob = knuth_random();       
      }
   else 
      {
      mut_prob = 0; // percent is per gene, not per individual
      }
   if (mut_prob < Mut_rate)
      {
          
      // if per gene, get amount of muts, if per indv there is only 1 mut.
      int num_mutations;        
      if (per_gene == 1)
         {
         num_mutations = poisson(((double)Mut_rate * (double)indv->length));
         }          
      else 
         {
         num_mutations = 1;
         }

      int k;
      for (k=num_mutations-1; k>=0; k--)
         {
         int end = start + length;
         int i;
         int j;
         int displaced_items[length];

         // remove items & shift
         for (i = end - 1; i >= start; i--)
            {
            displaced_items[i - start] = indv->genome[i];
            for (j = i; j < indv->length - 1; j++)
               {
               indv->genome[j] = indv->genome[j + 1];   
               }
            }       
         
         // re-add items at new location
         int new_start = uniform(indv->length - length);
         int new_end = new_start + length;
         int k = 0; // used for inversion
         for (i = new_end - 1; i >= new_start; i--)
            {
         
            // get value to insert, inverting order if necessary
            int value;
            if (invert == 1) 
               {
               value = displaced_items[k];
               } /* if */
            else 
               {
               value = displaced_items[i - new_start];
               } /* else */
         
            // insert value & shift       
            for (j = indv->length - 1; j >= new_start; j--)
               {
               if (j > new_start) 
                  {
                  indv->genome[j] = indv->genome[j - 1];
                  } /* if */
               else 
                  {
                  indv->genome[j] = value;
                  } /* else */
               } /* for j */
            if (invert) k++;
            } /* for i */
         } /* for */
      } /* if */

#ifdef DEBUG
   printf(" ---end displacement_mutate()---\n");
#endif
   }  /* displacement_mutate */

/********** poisson **********/
/* parameters:  lambda
   called by:   mutate()
   actions:     given a mutation probability, returns the number of
                mutations that will occur on a chromosome using the
                Poisson distribution.
                This algorithm was copied directly from the Lisp
                royal road code.
*/
int poisson(double lambda)
   {
   double lambda_term;
   double p;
   double sum;
   double unif_rand;
   int i;
 
#ifdef DEBUG
   printf("     ---in poisson---\n");
#endif
#ifdef DEBUG_MUT
printf(" lambda = %lf\n", lambda);
#endif
 
   if (lambda < 0.0)
      {
      printf(" Error(poisson): bad (neg.) lambda value.\n");
      return -1;
      }  /* if */
 
   unif_rand = knuth_random();
   lambda_term = exp(-lambda);
   sum = p = lambda_term;
#ifdef DEBUG_MUT
printf(" lambda = %lf, unif_rand = %lf\n", lambda, unif_rand);
printf(" sum = %e, lambda_term = %e, p = %e\n", sum, lambda_term, p);
#endif
 
   i = 0;
   while ((sum <= unif_rand)/* && (p > 0.0000005)*/)
      {
      i++;
      p = p * lambda / i;
      sum += p;
#ifdef DEBUG_MUT
printf(" i = %d, sum = %e, lambda_term = %e, p = %e\n",i,sum,lambda_term,p);
#endif
      }  /* while */
 
#ifdef DEBUG_MUT
printf(" i = %d\n", i);
#endif

#ifdef DEBUG
   printf("     ---end poisson---\n");
#endif
 
   return i;
   }  /* poisson */

/******************** begin homologous crossover routines ********************/

/********** homologous_crossover **********/
/* parameters:
   called by:   crossover_pop(), genop.c
   actions:     Given ptrs to two parents and two kids, perform
		one-point homologous crossover on parents to create
                new kids.
                Includes copying all of the parent data to the
                kid structure.
		If crossover results in a kid that is longer than
		Max_gen_len, the the crossover just does not occur.
		Returns 1 if xover happened, 0 if it did not happen.
   restricxns:  for variable length individuals only.

   note:	In all of the homologous crossover code, to simplify
		implementation, I treat the homology window as the Hx_window
		bits to the right of the current bit.  So the current bit
		is the leftmost bit of the homology window and not the
		center bit of the homology window.  The end result is the
		same either way.
*/
int homologous_crossover(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
                        INDIVIDUAL *kid1, INDIVIDUAL *kid2)
   {
   double cross_prob;
   int parent1_xpt, parent2_xpt;
   int random_in_window;
   int i, j;

#ifdef DEBUG
   printf(" ---in homologous_crossover()---\n");
#endif

   cross_prob = knuth_random();

#ifdef HX
   printf(" ********** potential homologous crossover \n");
   printf("     crossover rate: %lf, cross_prob: %lf\n",
	Xover_rate, cross_prob);
   printf("     parent1: ");
   print_genome(parent1->genome, 1);
   printf("     parent2: ");
   print_genome(parent2->genome, 1);
#endif

   if (cross_prob < Xover_rate)
      {
     /* attempt homologous crossover */
     /* within the homology window, this is where xover pt will be */
      random_in_window = uniform(Hx_window);
      parent1_xpt = uniform(parent1->length - Hx_window);
      parent2_xpt = get_parent2_xpt(parent1, parent2, parent1_xpt) +
			random_in_window;
      parent1_xpt += random_in_window;
#ifdef HX
      printf("     actual xover pts shifted %d - p1 at %d, p2 at %d\n",
		random_in_window, parent1_xpt, parent2_xpt);
#endif
      if (valid_crossover(parent1, parent2, parent1_xpt, parent2_xpt))
         {
         for (i=0; i<parent1_xpt; i++)
            kid1->genome[i] = parent1->genome[i];
         for (j=parent2_xpt, i=parent1_xpt; j<parent2->length; j++, i++)
            kid1->genome[i] = parent2->genome[j];
         for (j=0; j<parent2_xpt; j++)
            kid2->genome[j] = parent2->genome[j];
         for (i=parent1_xpt, j=parent2_xpt; i<parent1->length; j++, i++)
            kid2->genome[j] = parent1->genome[i];
        /* get new lengths */
         kid1->length = parent1_xpt + parent2->length - parent2_xpt;
         kid2->length = parent2_xpt + parent1->length - parent1_xpt;
        /* set parent info, parent 1 contributes first portion */
         kid1->parent1_index = parent1->index;
         kid1->parent2_index = parent2->index;
         kid2->parent1_index = parent2->index;
         kid2->parent2_index = parent1->index;
         }  /* if valid xover */
      else
         {
        /* don't cross over */
        /* no, don't crossover, just copy parents to kid structures */
         for (i=0; i<parent1->length; i++)
            kid1->genome[i] = parent1->genome[i];
         for (i=0; i<parent2->length; i++)
            kid2->genome[i] = parent2->genome[i];
        /* get new lengths */
         kid1->length = parent1->length;
         kid2->length = parent2->length;
        /* set parent info, set parent2 to -1 if clones (no xover) */
         kid1->parent1_index = parent1->index;
         kid1->parent2_index = -1;
         kid2->parent1_index = parent2->index;
         kid2->parent2_index = -1;
         }  /* else not a valid xover*/
      }  /* if cross_prob < Xover_rate */
   else
      {
     /* don't attempt homologous crossover, just clone parents */
     /* no, don't crossover, just copy parents to kid structures */
      for (i=0; i<parent1->length; i++)
         kid1->genome[i] = parent1->genome[i];
      for (i=0; i<parent2->length; i++)
         kid2->genome[i] = parent2->genome[i];
     /* get new lengths */
      kid1->length = parent1->length;
      kid2->length = parent2->length;
     /* set parent info, set parent2 to -1 if clones (no xover) */
      kid1->parent1_index = parent1->index;
      kid1->parent2_index = -1;
      kid2->parent1_index = parent2->index;
      kid2->parent2_index = -1;
      }  /* else cross_prob >= Xover_rate */

#ifdef HX
   getchar();
#endif

#ifdef DEBUG
   printf(" ---end homologous_crossover()---\n");
#endif
   }  /* homologous_crossover */

/********** get_parent2_xpt **********/
/* parameters:	parent1, parent2, parent1_xpt
   called by:	homologous_crossover()
   actions:	give a crossover point on parent 1, slide parent 2 along
		and look at the homology window (size dtmd by Hx_window)
		for the best match and return that point.  If there is
		only one best match, return that point; if there are 
		multiple best matches, randomly pick one of them.

   note:	In all of the homologous crossover code, to simplify
		implementation, I treat the homology window as the Hx_window
		bits to the right of the current bit.  So the current bit
		is the leftmost bit of the homology window and not the
		center bit of the homology window.  The end result is the
		same either way.
*/
int get_parent2_xpt(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
                        int parent1_xpt)
   {
   int i, j, k;
   int best_match;		/* number of matching bits of best match */
   int number_best_match;	/* number of locations with best_match */
   int p2_range;		/* slide parent2 from bit 0 to p2_range 
				   along parent 1 for comparison */
   int *loc;			/* keep track of number of matching bits
				   at each relative position of parent1 and
				   parent 2 */
   int l;			/* save loc value to return so loc can
				   be deallocated */
   int sum;

#ifdef DEBUG
   printf(" ---in get_parent2_xpt()---\n");
#endif

  /* initialize values */
   best_match = -1;
   number_best_match = 0;
   p2_range = parent2->length - Hx_window + 1;
   loc = (int *)malloc(p2_range * sizeof(int));

#ifdef HX
   printf("     parent1_xpt = %d\n", parent1_xpt);
   printf("     parent2_len = %d, Hx_window = %d\n",parent2->length,Hx_window);
   printf("     parents range from 0 to %d\n", p2_range-1);
#endif

  /* slide parent2 along parent1 one bit at a time and count the number
     of matching bits in the homology window.  Keep track of the highest
     match values and the number of locations with that highest match
     value and save the locations in loc[]. */
   for (i=p2_range-1; i>=0; i--)
      {
      sum = 0;
      for (j=Hx_window-1; j>=0; j--)
         {
         if (parent1->genome[parent1_xpt+j] == parent2->genome[i+j])
            sum++;
         }  /* for j */
      if (sum == best_match)
         {
         loc[number_best_match] = i;
         number_best_match++;
         }  /* if */
      else if (sum > best_match)
         {
         best_match = sum;
         loc[0] = i;
         number_best_match = 1;
         }  /* else if */

#ifdef HX
      printf("        p1(%d) ", parent1_xpt);
      print_genome(&parent1, 0);
      printf(" p2(%d) ", i);
      print_genome(&parent2, 0);
      printf(" match:%d best:%d numbest:%d\n",
		sum, best_match, number_best_match);
#endif
      }  /* for i */

  /* at this point should know what best match is and how many and
     which points have that match.  Now decide which location to return. */
   if (number_best_match == 1)
      {
#ifdef HX
   printf("     parent2_xpt = %d\n", loc[0]);
#endif
      l = loc[0];
      free(loc);
      return l;
      }  /* if */
   else  /* there was more than one best match */
      {
      j = uniform(number_best_match);
#ifdef HX
   printf("     parent2_xpt = %d\n", loc[j]);
#endif
      l = loc[j];
      free(loc);
      return l;
      }  /* else */

#ifdef DEBUG
   printf(" ---end get_parent2_xpt()---\n");
#endif
   }  /* get_parent2_xpt */

/********** valid_crossover **********/
/* parameters:	parent1, parent2, parent1_xpt, parent2_xpt
   called by:	homologous_crossover()
   actions:	Given parents and selected crossover points on each,
		checks to make sure that the offspring produced from
		such a crossover will be within allowed length retrictions.
		If both offspring do not exceed Max_gen_len, return 1.
		If either or both offspring exceed Max_gen_len, return 0.
*/
int valid_crossover(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
			int parent1_xpt, int parent2_xpt)
   {
   int x, y;

   x = parent1_xpt + parent2->length - parent2_xpt;
   y = parent2_xpt + parent1->length - parent1_xpt;

   if (x >= Max_gen_len || x < Min_gen_len)
      return 0;
   else if (y >= Max_gen_len || y < Min_gen_len)
      return 0;
   else  return 1;
   }  /* valid_crossover */

/******************** end homologous crossover routines ********************/
