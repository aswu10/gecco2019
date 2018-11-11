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
/*
double fmin(double a, double b)
   {
   return(a<=b ? a : b);
   }  /* fmin */

/********** fmax **********/
/* parameters:	a, b	two double numbers
   actions:	return the greater of the two values.
*/
/*
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


   // printf("int_in_list\n");
   // for (i=0; i < len; i++) {
       // printf("%d ", alist[i]);
   // }
   // printf("\n");
   // printf(" --- Press any key to continue ---\n");  fgetc(stdin);
   for (i=len-1; i>=0; i--)
      {
          
      if (alist[i] == val)  return 1;
      }  /* for i */

      return 0;
   }  /* int_in_list */
   
/********** idx_in_list **********/
/* parameters:	val		value we are checking
		alist		list of integers
		len		number of elements to check in alist
				(check elements 0 to len-1)
   created:	08.02.17RN
   called by:	check_mut_fx(), fxfrr.c   and related routines
   actions:	Given a value and an array, find the index of the value. 
        Returns -1 if value not found
*/
int idx_in_list(int val, int *alist, int len)
   {
   int i;

   for (i=len-1; i>=0; i--)
      {
      if (alist[i] == val)  return i;
      }  /* for i */
   return -1;
   }  /* idx_in_list */

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

