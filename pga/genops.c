/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* genops.c
   07.24.98.AW	Created.  Routines having to do with genetic operators,
		such as crossover and mutation.

   Routines:	crossover_pop()		07.27.98.AW
		mutate_pop()		07.27.98.AW
		onept_crossover()	07.27.98.AW
		twopt_crossover()	07.27.98.AW
		uniform_crossover()	07.27.98.AW
		switch_crossover()	02.01.99.AW
		mutate()		07.27.98.AW
		poisson()		07.27.98.AW
		mutate_pop_multichar()	03.08.01.AW
		mutate_multichar()	03.08.01.AW
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "types.h"
#include <stdlib.h>
#include "extern.h"
#include "genops.h"
#include "random.h"
#include "params.h"

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
      else if (!strcmp(Xover_type, "uniform"))
         for (i=0, j=0; i<num_pairs; i++, j+=2)
            uniform_crossover(Parents[j], Parents[j+1], Kids[j], Kids[j+1]);
      else if (!strcmp(Xover_type, "switch"))
         for (i=0, j=0; i<num_pairs; i++, j+=2)
            switch_crossover(Parents[j], Parents[j+1], Kids[j], Kids[j+1]);
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
   if (!strcmp(Base, "binary"))
      {
      for (i=0; i<Num_breeding; i++) mutate(Kids[i]);
	  }
   else if (!strcmp(Base, "alphabet"))
	  {
      for (i=0; i<Num_breeding; i++) mutate_alpha(Kids[i]);
	  }
   else if (!strcmp(Base, "multichar"))
	  {
      for (i=0; i<Num_breeding; i++) mutate_multichar(Kids[i]);
	  }
   else if (!strcmp(Base, "integers"))
	  {
      for (i=0; i<Num_breeding; i++) mutate_int(Kids[i]);
	  }
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
         kid1->genome[i] = parent1->genome[i];
         kid2->genome[i] = parent2->genome[i];
         }  /* for i */
      for (i=cross_point; i<parent1->length; i++)
         {
         kid1->genome[i] = parent2->genome[i];
         kid2->genome[i] = parent1->genome[i];
         }  /* for i */
     /* set parent info, parent 1 contributes first portion */
      kid1->parent1_index = parent1->index;
      kid1->parent2_index = parent2->index;
      kid2->parent1_index = parent2->index;
      kid2->parent2_index = parent1->index;

      if (Trace)
         {
         Gen.num_x++;
        /* set pointers from kids to parents */
         kid1->num_xover = 1;	kid1->xover_pts[0] = cross_point;
         kid2->num_xover = 1;	kid2->xover_pts[0] = cross_point;
        /* set pointers from parents to kids */
         parent1->kids[parent1->num_kids] = kid1->index;
         parent1->kids[parent1->num_kids + 1] = kid2->index;
         parent1->num_kids += 2;
         parent2->kids[parent2->num_kids] = kid1->index;
         parent2->kids[parent2->num_kids + 1] = kid2->index;
         parent2->num_kids += 2;
        /* 000728AW: set more pointers from kids to parents */
         kid1->num_segments = kid2->num_segments = 2;
         kid1->seg_parent[0] = kid2->seg_parent[1] = parent1->index;
         kid1->seg_parent[1] = kid2->seg_parent[0] = parent2->index;
         kid1->seg_start[0] = kid2->seg_start[0] = 0;
         kid1->seg_start[1] = kid2->seg_start[1] = cross_point;
         kid1->seg_len[0] = kid2->seg_len[0] = cross_point;
         kid1->seg_len[1] = kid2->seg_len[1] = parent1->length-cross_point;
         }  /* if Trace */

      }  /* if */
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

      if (Trace)
         {
         kid1->num_xover = 0;
         kid2->num_xover = 0;
        /* set pointers from parents to kids */
         parent1->kids[parent1->num_kids] = kid1->index;
         parent1->num_kids++;
         parent2->kids[parent2->num_kids] = kid2->index;
         parent2->num_kids++;
        /* 000728AW: set more pointers from kids to parents */
         kid1->num_segments = kid2->num_segments = 1;
         kid1->seg_parent[0] = parent1->index;
         kid2->seg_parent[0] = parent2->index;
         kid1->seg_start[0] = kid2->seg_start[0] = 0;
         kid1->seg_len[0] = kid2->seg_len[0] = parent1->length;
         }  /* if Trace */

      }  /* else */

#ifdef SMALLERSTEP
   if (cross_prob < Xover_rate)
      printf(" ********** one point crossover at %d\n", cross_point);
   else
      printf(" ********** no crossover\n");
   printf(" %3d ", parent1->index);
   print_genome(parent1->genome, parent1->length, 0);
   printf("  ");
   printf("%3d ", kid1->index);
   print_genome(kid1->genome, kid1->length, 1);
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
   print_genome(parent2->genome, parent2->length, 0);
   printf("  ");
   printf("%3d ", kid2->index);
   print_genome(kid2->genome, kid2->length, 1);
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
			 
			 
         kid1->genome[i] = parent1->genome[i];
         kid2->genome[i] = parent2->genome[i];
         }  /* for i */
      for (i=cross_point; i<cross_point2; i++)
         {
         kid1->genome[i] = parent2->genome[i];
         kid2->genome[i] = parent1->genome[i];
         }  /* for i */
      for (i=cross_point2; i<parent1->length; i++)
         {
         kid1->genome[i] = parent1->genome[i];
         kid2->genome[i] = parent2->genome[i];
         }  /* for i */
     /* set parent info, parent 1 contributes first portion */
      kid1->parent1_index = parent1->index;
      kid1->parent2_index = parent2->index;
      kid2->parent1_index = parent2->index;
      kid2->parent2_index = parent1->index;

      if (Trace)
         {
         Gen.num_x++;
         kid1->num_xover = 2;
         kid1->xover_pts[0] = cross_point;
         kid1->xover_pts[1] = cross_point2;
         kid2->num_xover = 2;
         kid2->xover_pts[0] = cross_point;
         kid2->xover_pts[1] = cross_point2;
        /* set pointers from parents to kids */
         parent1->kids[parent1->num_kids] = kid1->index;
         parent1->kids[parent1->num_kids + 1] = kid2->index;
         parent1->num_kids += 2;
         parent2->kids[parent2->num_kids] = kid1->index;
         parent2->kids[parent2->num_kids + 1] = kid2->index;
         parent2->num_kids += 2;
        /* 000728AW: set more pointers from kids to parents */
         kid1->num_segments = 3;
         kid2->num_segments = 3;
         kid1->seg_parent[0] = kid2->seg_parent[1] = kid1->seg_parent[2]
		= parent1->index;
         kid2->seg_parent[0] = kid1->seg_parent[1] = kid2->seg_parent[2]
		= parent2->index;
         kid1->seg_start[0] = kid2->seg_start[0] = 0;
         kid1->seg_start[1] = kid2->seg_start[1] = cross_point;
         kid1->seg_start[2] = kid2->seg_start[2] = cross_point2;
         kid1->seg_len[0] = kid2->seg_len[0] = cross_point;
         kid1->seg_len[1] = kid2->seg_len[1] = cross_point2-cross_point;
         kid1->seg_len[2] = kid2->seg_len[2] = parent1->length-cross_point2;
         }  /* if Trace */

      }  /* if */
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

      if (Trace)
         {
         kid1->num_xover = 0;
         kid2->num_xover = 0;
        /* set pointers from parents to kids */
         parent1->kids[parent1->num_kids] = kid1->index;
         parent1->num_kids++;
         parent2->kids[parent2->num_kids] = kid2->index;
         parent2->num_kids++;
        /* 000728AW: set more pointers from kids to parents */
         kid1->num_segments = kid2->num_segments = 1;
         kid1->seg_parent[0] = parent1->index;
         kid2->seg_parent[0] = parent2->index;
         kid1->seg_start[0] = kid2->seg_start[0] = 0;
         kid1->seg_len[0] = kid2->seg_len[0] = parent1->length;
         }  /* if Trace */

      }  /* else */

#ifdef SMALLERSTEP
   if (cross_prob < Xover_rate)
      printf(" ********** two point crossover at %d and %d\n",
		cross_point, cross_point2);
   else
      printf(" ********** no crossover\n");
   printf(" %3d ", parent1->index);
   print_genome(parent1->genome, parent1->length, 0);
   printf("  ");
   printf("%3d ", kid1->index);
   print_genome(kid1->genome, kid1->length, 1);
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
   print_genome(parent2->genome, parent2->length, 0);
   printf("  ");
   printf("%3d ", kid2->index);
   print_genome(kid2->genome, kid2->length, 1);
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
 
   if (Trace)
      {
      kid1->num_xover = 0;
      kid2->num_xover = 0;
      inher = 0;
      }  /* if Trace */
 
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
            kid1->genome[i] = parent1->genome[i];
            kid2->genome[i] = parent2->genome[i];
 
            if (Trace)
               {
               if (inher == -1)  /* xover occurred */
                  {
                  kid1->xover_pts[kid1->num_xover] = i;
                  kid2->xover_pts[kid2->num_xover] = i;
                  kid1->num_xover++;
                  kid2->num_xover++;
                  }
               inher = 1;
               }  /* if Trace */
 
            }  /* if */
         else
            {
           /* inher = -1; parent1 to kid2, parent2 to kid1 */
            kid1->genome[i] = parent2->genome[i];
            kid2->genome[i] = parent1->genome[i];
 
            if (Trace)
               {
               if (inher == 1)  /* xover occurred */
                  {
                  kid1->xover_pts[kid1->num_xover] = i;
                  kid2->xover_pts[kid2->num_xover] = i;
                  kid1->num_xover++;
                  kid2->num_xover++;
                  }
               inher = -1;
               }  /* if Trace */
 
            }  /* else */
         }  /* for i */
 
     /* set parent info, parent 1 contributes first portion */
      kid1->parent1_index = parent1->index;
      kid1->parent2_index = parent2->index;
      kid2->parent1_index = parent2->index;
      kid2->parent2_index = parent1->index;
 
      if (Trace)
         {
         Gen.num_x++;
        /* set pointers from parents to kids */
         parent1->kids[parent1->num_kids] = kid1->index;
         parent1->kids[parent1->num_kids + 1] = kid2->index;
         parent1->num_kids += 2;
         parent2->kids[parent2->num_kids] = kid1->index;
         parent2->kids[parent2->num_kids + 1] = kid2->index;
         parent2->num_kids += 2;
         }  /* if Trace */
 
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
 
      if (Trace)
         {
        /* set pointers from parents to kids */
         parent1->kids[parent1->num_kids] = kid1->index;
         parent1->num_kids++;
         parent2->kids[parent2->num_kids] = kid2->index;
         parent2->num_kids++;
         }  /* if Trace */
 
      }  /* else don't crossover */
 
#ifdef DEBUG
   printf(" ---end uniform_crossover()---\n");
#endif
   }  /* uniform_crossover */

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
void switch_crossover(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
			INDIVIDUAL *kid1, INDIVIDUAL *kid2)
   {
   int i;
   double cross_prob;	/* xover_rate */
   double x_prob;	/* uniform_x rate */
   int inher;		/* stores previous inheritance pattern */
			/* p1 -> o1 or p1 -> o2, for Trace purposes */
 
#ifdef DEBUG
   printf(" ---in switch_crossover()---\n");
#endif

   if (Trace)
      {
      kid1->num_xover = 0;
      kid2->num_xover = 0;
      }  /* if Trace */
 
   cross_prob = knuth_random();

   if (cross_prob < Xover_rate)
      {
     /* yes, do crossover */

     /* copye nucleotides from parents to kids */
     /* if inher = 1, p1-c1; if inher = -1 p1-c2 */
     /* first bit special case */
      kid1->genome[0] = parent1->genome[0];
      kid2->genome[0] = parent2->genome[0];

     /* rest of bits */
      inher = 1;
      for (i=1; i<parent1->length; i++)
         {
         x_prob = knuth_random();
         if (x_prob < Uniform_x)
            {
            inher = inher * -1;

            if (Trace)
               {
               kid1->xover_pts[kid1->num_xover] = i;
               kid2->xover_pts[kid2->num_xover] = i;
               kid1->num_xover++;
               kid2->num_xover++;
               }  /* if Trace */

            }  /* if */

         if (inher == 1)
            {
           /* inher = 1; parent1 to kid1, parent2 to kid2 */
            kid1->genome[i] = parent1->genome[i];
            kid2->genome[i] = parent2->genome[i];
            }  /* if */
         else if (inher == -1)
            {
           /* inher = -1; parent1 to kid2, parent2 to kid1 */
            kid1->genome[i] = parent2->genome[i];
            kid2->genome[i] = parent1->genome[i];
            }  /* else if */

         }  /* for i */

     /* set parent info, parent 1 contributes first portion */
      kid1->parent1_index = parent1->index;
      kid1->parent2_index = parent2->index;
      kid2->parent1_index = parent2->index;
      kid2->parent2_index = parent1->index;
 
      if (Trace)
         {
         Gen.num_x++;
        /* set pointers from parents to kids */
         parent1->kids[parent1->num_kids] = kid1->index;
         parent1->kids[parent1->num_kids + 1] = kid2->index;
         parent1->num_kids += 2;
         parent2->kids[parent2->num_kids] = kid1->index;
         parent2->kids[parent2->num_kids + 1] = kid2->index;
         parent2->num_kids += 2;
         }  /* if Trace */

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
 
      if (Trace)
         {
        /* set pointers from parents to kids */
         parent1->kids[parent1->num_kids] = kid1->index;
         parent1->num_kids++;
         parent2->kids[parent2->num_kids] = kid2->index;
         parent2->num_kids++;
         }  /* if Trace */

      }  /* else don't crossover */

#ifdef DEBUG
   printf(" ---end switch_crossover()---\n");
#endif
   }  /* switch_crossover */

/********** mutate **********/
/* parameters:
   called by:	mutate_pop(), genop.c
   actions:	given an individual, mutate it.
		Use a poisson distribution to decide how many mutations
		should be performed on this individual.  Then chooses
		mutation locations randomly, with replacement.
*/
void mutate(INDIVIDUAL *indv)
   {
   int i;
   int num_mutations;
   int random_num;
#ifdef SMALLERSTEP
   int *oldgenome;
   char *mask;
   int *muts;
#endif

#ifdef DEBUG
   printf(" ---in mutate()---\n");
#endif

   num_mutations = poisson(((double)Mut_rate * (double)indv->length));

#ifdef SMALLERSTEP
   muts = (int *)malloc(num_mutations * sizeof(int));
   oldgenome = (int *)malloc(indv->length * sizeof(int));
   mask = (char *)malloc(indv->length * sizeof(char));
   for (i=0; i<indv->length; i++)
      {
      oldgenome[i] = indv->genome[i];
      mask[i] = ' ';
      }  /* for */
#endif

   if (Trace)  indv->num_mut = num_mutations;

   for (i=num_mutations-1; i>=0; i--)
      {
      random_num = uniform(indv->length);
      if (indv->genome[random_num] == '0')  indv->genome[random_num] = '1';
      else  indv->genome[random_num] = '0';

      if (Trace)  indv->mut_pts[i] = random_num;

#ifdef SMALLERSTEP
      muts[i] = random_num;
      mask[random_num] = '|';
#endif
      }  /* for i */

#ifdef SMALLERSTEP
   if (num_mutations == 0)
      {
      printf(" ***** expect %6.4lf * %d = %5.2lf ***** no mutations\n", 
		Mut_rate, indv->length,
		(Mut_rate * (double)indv->length) );
      printf(" %3d ", indv->index);
      print_genome(indv->genome, indv->length, 1);
      }  /* if */
   else
      {
      printf(" ***** expect %6.4lf * %d = %5.2lf ***** %d mutation(s) at:",
		Mut_rate, indv->length,
		(Mut_rate * (double)indv->length) , num_mutations);
      for (i=0; i<num_mutations; i++)  printf(" %d", muts[i]);
      printf("\n");
      printf(" %3d ", indv->index);
      print_genome(oldgenome, indv->length, 1);
      printf("     ");
      print_genome(mask, indv->length, 1);
      printf("     ");
      print_genome(indv->genome, indv->length, 1);
      }  /* else */
   printf(" --- Press any key to continue ---\n");  fgetc(stdin);
   free(muts);
   free(oldgenome);
   free(mask);
#endif  /* SMALLERSTEP */

#ifdef DEBUG
   printf(" ---end mutate()---\n");
#endif
   }  /* mutate */

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
   print_genome(parent1->genome, parent1->length, 1);
   printf("     parent2: ");
   print_genome(parent2->genome, parent2->length, 1);
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

         if (Trace)
            {
            Gen.num_x++;
            kid1->num_xover = 1;   kid1->xover_pts[0] = parent1_xpt;
            kid2->num_xover = 1;   kid2->xover_pts[0] = parent2_xpt;
           /* set pointers from parents to kids */
            parent1->kids[parent1->num_kids] = kid1->index;
            parent1->kids[parent1->num_kids + 1] = kid2->index;
            parent1->num_kids += 2;
            parent2->kids[parent2->num_kids] = kid1->index;
            parent2->kids[parent2->num_kids + 1] = kid2->index;
            parent2->num_kids += 2;
           /* 000728AW: set more pointers from kids to parents */
            kid1->num_segments = kid2->num_segments = 2;
            kid1->seg_parent[0] = kid2->seg_parent[1] = parent1->index;
            kid1->seg_parent[1] = kid2->seg_parent[0] = parent2->index;
            kid1->seg_start[0] = kid2->seg_start[0] = 0;
            kid1->seg_len[0] = parent1_xpt;
            kid2->seg_len[0] = parent2_xpt;
            kid1->seg_start[1] = parent2_xpt;
            kid1->seg_len[1] = parent2->length - parent2_xpt;
            kid2->seg_start[1] = parent1_xpt;
            kid2->seg_len[1] = parent1->length - parent1_xpt;
            }  /* if Trace */

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

         if (Trace)
            {
            kid1->num_xover = 0;
            kid2->num_xover = 0;
           /* set pointers from parents to kids */
            parent1->kids[parent1->num_kids] = kid1->index;
            parent1->num_kids++;
            parent2->kids[parent2->num_kids] = kid2->index;
            parent2->num_kids++;
           /* 000728AW: set more pointers from kids to parents */
            kid1->num_segments = kid2->num_segments = 1;
            kid1->seg_parent[0] = parent1->index;
            kid2->seg_parent[0] = parent2->index;
            kid1->seg_start[0] = kid2->seg_start[0] = 0;
            kid1->seg_len[0] = parent1->length;
            kid2->seg_len[0] = parent2->length;
            }  /* if Trace */
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

      if (Trace)
         {
         kid1->num_xover = 0;
         kid2->num_xover = 0;
        /* set pointers from parents to kids */
         parent1->kids[parent1->num_kids] = kid1->index;
         parent1->num_kids++;
         parent2->kids[parent2->num_kids] = kid2->index;
         parent2->num_kids++;
        /* 000728AW: set more pointers from kids to parents */
         kid1->num_segments = kid2->num_segments = 1;
         kid1->seg_parent[0] = parent1->index;
         kid2->seg_parent[0] = parent2->index;
         kid1->seg_start[0] = kid2->seg_start[0] = 0;
         kid1->seg_len[0] = parent1->length;
         kid2->seg_len[0] = parent2->length;
         }  /* if Trace */

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
      print_genome(&parent1->genome[parent1_xpt], Hx_window, 0);
      printf(" p2(%d) ", i);
      print_genome(&parent2->genome[i], Hx_window, 0);
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

/******************** start alphabet routines ********************/

// /********** mutate_pop_alpha **********/
// /* parameters:
   // called by:	reproduce(), reproduce.c
   // actions:
// */
// void mutate_pop_alpha()
   // {
   // int i;
 
// #ifdef DEBUG
   // printf(" ---in mutate_pop_alpha()---\n");
// #endif
 
   // for (i=0; i<Num_breeding; i++)
      // {
      // mutate_alpha(Kids[i]);
      // }  /* for i */
 
// #ifdef DEBUG
   // printf(" ---end mutate_pop_alpha()---\n");
// #endif
   // }  /* mutate_pop_alpha */

/********** mutate_alpha **********/
/* parameters:
   called by:	mutate_pop(), genop.c
   actions:	given an individual, mutate it.
		Use a poisson distribution to decide how many mutations
		should be performed on this individual.  Then chooses
		mutation locations randomly, with replacement.
*/
void mutate_alpha(INDIVIDUAL *indv)
   {
   int i;
   int num_mutations;
   int random_num;
   int random_num2;
#ifdef SMALLERSTEP
   char *oldgenome;
   char *mask;
   int *muts;
#endif

#ifdef DEBUG
   printf(" ---in mutate_alpha()---\n");
#endif

   num_mutations = poisson((double)(Mut_rate * indv->length));

#ifdef SMALLERSTEP
   muts = (int *)malloc(num_mutations * sizeof(int));
   oldgenome = (char *)malloc(indv->length * sizeof(char));
   mask = (char *)malloc(indv->length * sizeof(char));
   for (i=0; i<indv->length; i++)
      {
      oldgenome[i] = indv->genome[i];
      mask[i] = ' ';
      }  /* for */
#endif

   if (Trace)  indv->num_mut = num_mutations;

   for (i=num_mutations-1; i>=0; i--)
      {
      random_num = uniform(indv->length);
      random_num2 = uniform(Alphabet_size);
/*
printf(" old value: (%d,", indv->genome[random_num]);
putchar(indv->genome[random_num]);
*/
      indv->genome[random_num] = 
		((indv->genome[random_num] -97 + random_num2) % Alphabet_size)
		+ 97;
/*
printf(") - %d - new: (%d,", random_num2, indv->genome[random_num]);
putchar(indv->genome[random_num]);
printf(")\n");
getchar();
*/

      if (Trace)  indv->mut_pts[i] = random_num;

#ifdef SMALLERSTEP
      muts[i] = random_num;
      mask[random_num] = '|';
#endif
      }  /* for i */

#ifdef SMALLERSTEP
   if (num_mutations == 0)
      {
      printf(" ********** no mutations\n");
      printf(" %3d ", indv->index);
      print_genome(indv->genome, indv->length, 1);
      }  /* if */
   else
      {
      printf(" ********** %d mutation(s) at:", num_mutations);
      for (i=0; i<num_mutations; i++)  printf(" %d", muts[i]);
      printf("\n");
      printf(" %3d ", indv->index);
      print_genome(oldgenome, indv->length, 1);
      printf("     ");
      print_genome(mask, indv->length, 1);
      printf("     ");
      print_genome(indv->genome, indv->length, 1);
      }  /* else */
   printf(" --- Press any key to continue ---\n");  fgetc(stdin);
   free(muts);
   free(oldgenome);
   free(mask);
#endif  /* SMALLERSTEP */

#ifdef DEBUG
   printf(" ---end mutate_alpha()---\n");
#endif
   }  /* mutate_alpha */

/******************** end alphabet routines ********************/

/******************** start multichar routines ********************/

// /********** mutate_pop_multichar **********/
// /* parameters:
   // called by:   reproduce(), reproduce.c
   // actions:
// */
// void mutate_pop_multichar()
   // {
   // int i;

// #ifdef DEBUG
   // printf(" ---in mutate_pop_multichar()---\n");
// #endif

   // for (i=0; i<Num_breeding; i++)
      // {
      // mutate_multichar(Kids[i]);
      // }  /* for i */

// #ifdef DEBUG
   // printf(" ---end mutate_pop_multichar()---\n");
// #endif
   // }  /* mutate_pop_multichar */

/********** mutate_multichar **********/
/* parameters:
   called by:   mutate_pop(), genop.c
   actions:     given an individual, mutate it.
                Use a poisson distribution to decide how many mutations
                should be performed on this individual.  Then chooses
                mutation locations randomly, with replacement.
*/
void mutate_multichar(INDIVIDUAL *indv)
   {
   int i;
   int num_mutations;
   int random_num;
   int random_num2;
#ifdef SMALLERSTEP
   char *oldgenome;
   char *mask;
   int *muts;
#endif

#ifdef DEBUG
   printf(" ---in mutate_multichar()---\n");
#endif

   num_mutations = poisson((double)(Mut_rate * indv->length));

#ifdef SMALLERSTEP
   muts = (int *)malloc(num_mutations * sizeof(int));
   oldgenome = (char *)malloc(indv->length * sizeof(char));
   mask = (char *)malloc(indv->length * sizeof(char));
   for (i=0; i<indv->length; i++)
      {
      oldgenome[i] = indv->genome[i];
      mask[i] = ' ';
      }  /* for */
#endif

   if (Trace)  indv->num_mut = num_mutations;

   for (i=num_mutations-1; i>=0; i--)
      {
      random_num = uniform(indv->length);
      random_num2 = uniform(Alphabet_size);
/*
printf(" old value: (%d,", indv->genome[random_num]);
putchar(indv->genome[random_num]);
*/
      indv->genome[random_num] = Xters[random_num2];
/*
printf(") - %d - new: (%d,", random_num2, indv->genome[random_num]);
putchar(indv->genome[random_num]);
printf(")\n");
getchar();
*/

      if (Trace)  indv->mut_pts[i] = random_num;

#ifdef SMALLERSTEP
      muts[i] = random_num;
      mask[random_num] = '|';
#endif
      }  /* for i */

#ifdef SMALLERSTEP
   if (num_mutations == 0)
      {
      printf(" ********** no mutations\n");
      printf(" %3d ", indv->index);
      print_genome(indv->genome, indv->length, 1);
      }  /* if */
   else
      {
      printf(" ********** %d mutation(s) at:", num_mutations);
      for (i=0; i<num_mutations; i++)  printf(" %d", muts[i]);
      printf("\n");
      printf(" %3d ", indv->index);
      print_genome(oldgenome, indv->length, 1);
      printf("     ");
      print_genome(mask, indv->length, 1);
      printf("     ");
      print_genome(indv->genome, indv->length, 1);
      }  /* else */
   printf(" --- Press any key to continue ---\n");  fgetc(stdin);
   free(muts);
   free(oldgenome);
   free(mask);
#endif  /* SMALLERSTEP */

#ifdef DEBUG
   printf(" ---end mutate_multichar()---\n");
#endif
   }  /* mutate_multichar */

/******************** end multichar routines ********************/

/******************** start int routines ********************/

void mutate_int(INDIVIDUAL *indv)
   {
   int i;
   int num_mutations;
   int random_num;
#ifdef SMALLERSTEP
   int *oldgenome;
   char *mask;
   int *muts;
#endif

#ifdef DEBUG
   printf(" ---in mutate()---\n");
#endif


   num_mutations = poisson(((double)Mut_rate * (double)indv->length));

   // printf("Mut_rate: %lf\n", Mut_rate);
   // printf("indv->length: %d\n", indv->length);
   // printf("num_mutations: %d\n", num_mutations);
   // printf(" --- Press any key to continue ---\n");  fgetc(stdin);
   
#ifdef SMALLERSTEP
   muts = (int *)malloc(num_mutations * sizeof(int));
   oldgenome = (int *)malloc(indv->length * sizeof(int));
   mask = (char *)malloc(indv->length * sizeof(char));
   for (i=0; i<indv->length; i++)
      {
      oldgenome[i] = indv->genome[i];
      mask[i] = ' ';
      }  /* for */
#endif

   if (Trace)  indv->num_mut = num_mutations;

   for (i=num_mutations-1; i>=0; i--)
      {
      random_num = uniform(indv->length);
      int new_value = 1 + uniform(Alphabet_size);
      
	  if (funiform(1) > 0.5) 
		    {
	        indv->genome[random_num] = new_value;		 
		    }
	  else
            {
			indv->genome[random_num] = -1 * new_value;		
			}

      if (Trace)  indv->mut_pts[i] = random_num;

#ifdef SMALLERSTEP
      muts[i] = random_num;
      mask[random_num] = '|';
#endif
      }  /* for i */

#ifdef SMALLERSTEP
   if (num_mutations == 0)
      {
      printf(" ***** expect %6.4lf * %d = %5.2lf ***** no mutations\n", 
		Mut_rate, indv->length,
		(Mut_rate * (double)indv->length) );
      printf(" %3d ", indv->index);
      print_genome(indv->genome, indv->length, 1);
      }  /* if */
   else
      {
      printf(" ***** expect %6.4lf * %d = %5.2lf ***** %d mutation(s) at:",
		Mut_rate, indv->length,
		(Mut_rate * (double)indv->length) , num_mutations);
      for (i=0; i<num_mutations; i++)  printf(" %d", muts[i]);
      printf("\n");
      printf(" %3d ", indv->index);
      print_genome(oldgenome, indv->length, 1);
      printf("     ");
      print_genome(mask, indv->length, 1);
      printf("     ");
      print_genome(indv->genome, indv->length, 1);
      }  /* else */
   printf(" --- Press any key to continue ---\n");  fgetc(stdin);
   free(muts);
   free(oldgenome);
   free(mask);
#endif  /* SMALLERSTEP */

#ifdef DEBUG
   printf(" ---end mutate()---\n");
#endif
   }  /* mutate */
