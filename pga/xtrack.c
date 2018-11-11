/**** Copyright 2003, Annie S. Wu, All Rights Reserved ****/

/* cog.c
   031208AW	Created.  Routines that track xter locations and count.

   031208.AW	First attempt to calculate distance between identical
		characters.  Assume circular chromosome.  Found that
		the average distance will always be equal to the number
		of unique characters on an individual.  This conclusion
		confirmed mathematically.  Save this code as xtrack.c.v1.

		Instead of circular chromosome, assume linear chromosome,
		read left to right.  Characters that do not have an identical
		to the right will have distance zero.  Sum two sets of values:
		sum_all includes the zero distances and
		sum_not_zero includes only non-zero distances.
*/

#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "types.h"
#include "extern.h"
#include "xtrack.h"

#define DEBUG_XTRACK 1
#define DEBUG_XTRACK_DETAIL 1
#define DEBUG_GEN 1
#undef DEBUG_XTRACK
#undef DEBUG_XTRACK_DETAIL
#undef DEBUG_GEN

/********** xtrack_dist_indv **********/
/* parameters:  genome, start_bit, end_bit, xtrack
   created:     031208AW
   called by:
   actions:     Calculate avg and stdev of distances between every two
		"adjacent" identical xters on an individual.
		Adjacent means there are none of the same xter in between.
*/
void xtrack_dist_indv(char *genome, int start_bit, int end_bit,
			XTER_TRACK *xtrack)
   {
   int count;		/* counts all characters on genome */
   int next; 		/* search for next identical xter */
   int dist;

  /* count all bits */
   int sum_all;		/* sum of distances */
   int sum_all2;	/* sum of distances squared */
   int length;		/* length of segment to search */

  /* count only bits with non zero dist to next bit */
   int sum_not_zero;
   int sum_not_zero2;
   int length_not_zero;  /* count number of bits with non-zero dist */

#ifdef DEBUG
   printf(" ---in xtrack_dist_indv---\n");
#endif

#ifdef DEBUG_XTRACK
printf(" enter xtrack_dist_indv: start %d, end %d\n",
       start_bit, end_bit);
 printf(" %s\n", genome);
#endif

   sum_all = sum_all2 = 0;
   length = end_bit - start_bit + 1;

   sum_not_zero = sum_not_zero2 = 0;
   length_not_zero = 0;

   for (count=start_bit; count<end_bit; count++)
      {
      dist = 0;
      for (next=count+1; next<=end_bit; next++)
         {
#ifdef DEBUG_XTRACK_DETAIL
printf("        ");
printf(" compare genome[%d] and genome[%d] distance %d\n", count+start_bit,
((count+next)%length)+start_bit, next);
#endif
         if (genome[count+start_bit] == genome[((count+next)%length)+start_bit])
            {
            dist = next-count;
#ifdef DEBUG_XTRACK_DETAIL
printf("        ");
printf(" next = %d\n", next);
#endif
            break;
            }  /* if */
         }  /* for next */

      if (dist > xtrack->dist_max)  xtrack->dist_max = dist;
      if (dist < xtrack->dist_min)  xtrack->dist_min = dist;

      sum_all += dist;
      sum_all2 += dist * dist;

      if (dist > 0)
         {
         length_not_zero++;
         sum_not_zero += dist;
         sum_not_zero2 += dist * dist;
         }  /* if */

      }  /* for count */

   xtrack->num = length_not_zero;
   xtrack->dist_sum = sum_not_zero;
   xtrack->dist_sum2 = sum_not_zero2;
   xtrack->dist_avg = (double)sum_not_zero/(double)length;
   xtrack->dist_sd = sqrt(fabs(sum_not_zero2 -
			sum_not_zero*sum_not_zero/(double)length)/
			(double)(length - 1) );
/*
   xtrack->num = length;
   xtrack->dist_sum = sum_all;
   xtrack->dist_sum2 = sum_all2;
   xtrack->dist_avg = (double)sum_all/(double)length;
   xtrack->dist_sd = sqrt(fabs(sum_all2 - sum_all*sum_all/(double)length)/
			(double)(length - 1) );
*/
 
#ifdef DEBUG_XTRACK
printf(" leave xtrack_dist_indv: start %d, end %d, sum %d, avg %lf, sd %lf max %d min %d sum2 %d num %d\n",
       start_bit, end_bit, xtrack->dist_sum, xtrack->dist_avg, xtrack->dist_sd,
	xtrack->dist_max, xtrack->dist_min, xtrack->dist_sum2, xtrack->num);
// printf(" %s\n", genome);
#endif
 
#ifdef DEBUG
   printf(" ---end xtrack_dist_indv---\n");
#endif
   }  /* xtrack_dist_indv */

/********** xtrack_gen_start **********/
/* parameters:  
   created:     031227AW
   called by:   gen_start(), stats.c
   actions:     Initialize xtrack data for current generation.
*/              
void xtrack_gen_start()
   {
   Gen.avg_xtrack->dist_sum = 0;
   Gen.avg_xtrack->dist_avg = -1;
   Gen.avg_xtrack->dist_sd = -1;
   Gen.avg_xtrack->dist_max = INT_MIN;
   Gen.avg_xtrack->dist_min = INT_MAX;
   Gen.best_xtrack->dist_sum = 0;
   Gen.best_xtrack->dist_avg = -1;
   Gen.best_xtrack->dist_sd = -1;
   Gen.best_xtrack->dist_max = INT_MIN;
   Gen.best_xtrack->dist_min = INT_MAX;
   }  /* xtrack_gen_start */

/********** xtrack_gen_end **********/
/* parameters:  
   created:     031208AW
   called by:   gen_end(), stats.c
   actions:     Calculate avg of xtrack data from entire generation.
*/              
void xtrack_gen_end()
   {                    
   int i;
   int sum;             /* sum of distances */
   int sum2;            /* sum of distances */
   int num;		/* num of pairs, should = # bits in pop */
   
#ifdef DEBUG
   printf(" ---in xtrack_dist_indv---\n");
#endif

  /* get data for entire population */

   num = sum = sum2 = 0;

   for (i=Pop_size-1; i>=0; i--)
      {
      sum += Pop[i]->xtrack->dist_sum;
      sum2 += Pop[i]->xtrack->dist_sum2;
      num += Pop[i]->xtrack->num;

      if (Pop[i]->xtrack->dist_max > Gen.avg_xtrack->dist_max)
         Gen.avg_xtrack->dist_max = Pop[i]->xtrack->dist_max;
      if (Pop[i]->xtrack->dist_min < Gen.avg_xtrack->dist_min)
         Gen.avg_xtrack->dist_min = Pop[i]->xtrack->dist_min;
      }  /* for i */

   Gen.avg_xtrack->num = num;
   Gen.avg_xtrack->dist_sum = sum;
   Gen.avg_xtrack->dist_sum2 = sum2;
   Gen.avg_xtrack->dist_avg = (double)sum/(double)num;
   Gen.avg_xtrack->dist_sd = sqrt(fabs(sum2 - sum*sum/(double)num)/
				(double)(num - 1) );

#ifdef DEBUG_GEN
printf(" ----- Generation %d\n", Gen.index);
printf("       num %d sum %d sum2 %d avg %lf sd %lf\n",
		Gen.avg_xtrack->num,
		Gen.avg_xtrack->dist_sum,
		Gen.avg_xtrack->dist_sum2,
		Gen.avg_xtrack->dist_avg,
		Gen.avg_xtrack->dist_sd);
#endif

  /* get data for just best individual in population */
   Gen.best_xtrack->dist_avg = Pop[Gen.best_indv_index]->xtrack->dist_avg;
   Gen.best_xtrack->dist_sd = Pop[Gen.best_indv_index]->xtrack->dist_sd;
   Gen.best_xtrack->dist_max = Pop[Gen.best_indv_index]->xtrack->dist_max;
   Gen.best_xtrack->dist_min = Pop[Gen.best_indv_index]->xtrack->dist_min;

#ifdef DEBUG
   printf(" ---end xtrack_dist_indv---\n");
#endif
   }  /* xtrack_gen_end */
