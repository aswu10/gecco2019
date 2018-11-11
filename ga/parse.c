/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* This file contains the routines to parse a text string using
   the state tree created by build_goto and build_fail.
   This is algorithm 1 of Aho and Corasick.

   08/06/94 AW  Created.
   08/06/94 AW  parse works!
   08/23/94 AW  moved to ~aswu/rr directory and modified to work
		with our GA program.  Added run_stats data collection.
		Seems to be parsing the individuals correctly.
   06/01/95 AW  Changed printout in #ifdef MARK to print numbers underneath
		all of the bits of a building block, not just the first bit.
		Haven't checked to see if it works yet.
   06/02/95 AW	Randomizes all non-coding bits if RANDOMIZE is defined
		in types.h.
   06/06/95 AW  Instead of randomizing all of the non-coding bits,
		sets them all to one.
   09.18.98.AW	Modified.
   09.24.98.AW	Change parse() so that there is an input parameter that
		is a string of xters the same length as the individuals
		to be parsed.  Expects this string to be initialized to
		a neutral xter, then parse() will print the (one xter) name
		of any string in the location of the start of the string.
        	Right now the name of each character is a number.  May have
		to switch to xters later, but leave numbers for now.  If
		two different lengthed bbs start at the same location, I
		think the longer one will overwrite the shorter one.
   10.09.98.AW	See comment in build.c.  I am now using the variable
		Ascii_offset instead of hardcoding the value 48 to convert
		from character to integer value.  This is because the
		conversion from alphabetic xters to integer values (assuming
		that a is 1, b is 2, c is 3, etc.) is 97, so using a variable
		allows this to be set differently for binary bases and for
		alphabetic bases.  Ascii_offset is set in init_run(), init.c.
		***** Important note: Ascii_offset refers to conversion of
		bases to integers, not the conversion of bb names to integers.
*/

#include <stdio.h>
#include <stdlib.h>
#include "types.h"
#include "extern.h"
#include "parse.h"
#include "random.h"

#define PARSE
#undef  PARSE

/********** parse **********/
/* parameters:	tree		goto and fail functions
		textstring	string to parse
		textlen		length of textstring
		marker		where to mark found strings
   called by:	eval_indv()
   actions:	Searches for predefined patterns or building blocks in a
		text string. 
		Keeps track of where building blocks are found by printing
		to marker where the start of each existing building block is.
		Expects marker to be allocated and filled with a neutral xter.
*/
void parse(struct state_type *tree, char *textstring, int textlen, char *marker)
   {
   struct state_type *state;
   int i, j;

#ifdef DEBUG
   printf(" ---in parse---\n");
#endif
#ifdef PARSE
   printf(" textstring = %s, textlen = %d\n", textstring, textlen);
#endif

   state = tree;
#ifdef PARSE
   printf("    AT state %d\n", state->number);
#endif

   for (i=0; i<textlen; i++)
      { 
#ifdef PARSE
      printf(" Next xter = %s\n", &textstring[i]);
#endif
      while (state->g[textstring[i]-Ascii_offset] == NULL)
         {
#ifdef PARSE
         printf("     state %d -> [f] -> state %d\n", state->number,
                state->f->number);
#endif
         state = state->f;
         }

#ifdef PARSE
      printf("     state %d -> [%d] -> state %d\n", state->number,
		textstring[i]-Ascii_offset,
		state->g[textstring[i]-Ascii_offset]->number);
#endif
      state = state->g[textstring[i]-Ascii_offset];

      if (state->output[0] >= 0)
         {
        /***** print individual with the building block number underneath
               the first (left-most) bit of the building block only *****/
        /* didn't use Ascii_offset here because we are converting bb
           names to integers here, not bases to integers - this
           means we can't have more than ten bb, may switch to xters later.
           This use of 48 requires that I also hard code 48 in eval_indv()
           in fxfrr.c where I refer to the name of the bb */
        /* change 48 to 97 -- use small letters to name the bbs --
           now can have 26 bb.  Also change in eval_indv() in fxfrr.c */
         marker[i-(state->length[0])+1] = state->output[0]+97;
#ifdef PARSE
         printf(" Found string %d starting at bit %d \n", state->output[0], 
                 i-state->length[0]+1);
#endif
         }
      }  /* for i */

#ifdef DEBUG
   printf(" ---end parse---\n");
#endif
   }  /* parse */
