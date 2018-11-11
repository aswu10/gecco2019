/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* Utility functions

   09/19/93 AW  Created.
*/
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#include "types.h"
#include "util.h"

/********** fmin **********/
/* parameters:	a, b	two double numbers
   actions:	return the lesser of the two values.
*/
/***** 031208AW Comment out and use same fxn in math.h
double fmin(double a, double b)
   {
   return(a<=b ? a : b);
   }  /* fmin */

/********** fmax **********/
/* parameters:	a, b	two double numbers
   actions:	return the greater of the two values.
*/
/***** 031208AW Comment out and use same fxn in math.h
double fmax(double a, double b)
   {
   return(a>=b ? a : b);
   }  /* fmax */

/********** get_date_time **********/
/* parameters:	none
   actions:	returns a string with the date and time.
		*** This routine doesn't work yet *** 9/29/93
*/
char *get_date_time()
   {
   char time_string[100];
   time_t *tstart;

   tstart = (time_t *)malloc(sizeof(time_t));
   time(tstart);
   strftime(time_string, 99, " %m/%d/%y  %H:%M:%S ", localtime(tstart));
   return (time_string);
   }  /* get_date_time */

/********** int_in_list **********/
/* parameters:	val		value we are checking
		alist		list of integers
		len		number of elements to check in alist
				(check elements 0 to len-1)
   created:	990112AW
   called by:	check_mut_fx(), fxfrr.c   and related routines
   actions:	Given a value and an array, check to see if any other
		element of the array has the same value.  Returns 1
		if the integer value was found in the list, returns 0
		if not.
		Routine used for checking for multiple same mutations in
		list of mutations.
*/
int int_in_list(int val, int *alist, int len)
   {
   int i;

   for (i=len-1; i>=0; i--)
      {
      if (alist[i] == val)  return 1;
      }  /* for i */
   return 0;
   }  /* int_in_list */

/********** which_parent **********/
/* parameters:  indv, bit
   created:     990113AW
   called by:
   actions:     given a bit and an individual in which crossover has
                occurred, return whether that bit came from that
                individual's parent1 or parent2.
 
                Right now will only work with one or two crossovers.
*/
int which_parent(INDIVIDUAL *indv, int bit)
   {
   if (indv->num_xover == 1)
      {
      if (bit < indv->xover_pts[0])  return 1;
      else  return 2;
      }
   else if (indv->num_xover == 2)
      {
      if (bit < indv->xover_pts[0] || bit >= indv->xover_pts[1])
         return 1;
      else  return 2;
      }
   }  /* which_parent */

/********** bin_to_int_pos **********/
/* parameters:  bstring         character string
                left_bit        left most bit of substring
                right_bit       right most bit of substring
   note:        left_bit and right_bit count from zero.
   actions:     evaluates a substring of bstring which is
                expected to be binary and returns the integer
                value.  The length of the substring to be
                evaluated is right_bit - left_bit + 1
                and if the entire bstring is to be evaluated,
                then enter 0 for left_bit and length-1 for
                right_bit.

                Converts to positive numbers only.  No negative numbers.
		Zero is acceptable.
*/
int bin_to_int_pos(char *bstring, int left_bit, int right_bit)
   {
   int i, j, sum;
 
   sum = 0;

   j = 0;
   for (i=right_bit; i>=left_bit; i--)
      {
      sum += ( (bstring[i] == '0') ? 0 : ipow(2,j) );
      j++;
      }  /* for */
   return sum;

   }  /* bin_to_int_pos */
 
/********** ipow **********/
int ipow(int x, int y)
   {
   int i, j;
   j = 1;
   for (i=0; i<y; i++)
      j = j * x;
   return j;
   }  /* ipow */

/********** bin_to_int **********/
/* parameters:  bstring         character string
                left_bit        left most bit of substring
                right_bit       right most bit of substring
   note:        left_bit and right_bit count from zero.
   actions:     evaluates a substring of bstring which is
                expected to be binary and returns the integer
                value.  The length of the substring to be
                evaluated is right_bit - left_bit + 1
                and if the entire bstring is to be evaluated,
                then enter 0 for left_bit and length-1 for
                right_bit.

   030801AW     Calculates negative numbers using 1's complement.
		Simply flip bits.
		2's complement requires: subtract 1, then flip bits.
		Did not feel like figuring out subtract.
*/
int bin_to_int(char *bstring, int left_bit, int right_bit)
   {
   int i, j;
   int length;
   char *tempstring;

   length = right_bit - left_bit + 1;
   tempstring = (char *)malloc( (length + 1) * sizeof(char) );

   if (bstring[left_bit] == '0')
      return ( bin_to_int_pos(bstring, left_bit, right_bit) );
   else if (bstring[left_bit] == '1')
      {
      for (j=0, i=left_bit; i<=right_bit; i++, j++)
         {
         if (bstring[i] == '1')        tempstring[j] = '0';
         else if (bstring[i] == '0')   tempstring[j] = '1';
         else
            printf(" Error(bin_to_int):  String not binary.\n");
         }  /* fir */
      return ( -1 * bin_to_int_pos(tempstring, 0, length-1) );
      }  /* else if */
   else
      printf(" Error(bin_to_int):  String not binary.\n");
   }  /* bin_to_int */
