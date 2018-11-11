/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* This file contains the code for building a finite state
   machine that will recognize a specific set of keywords.
   It is basically algorithms 2 and 3 from (Aho & Corasick, 1975).

   07/28/94 AW  Created.
   07/29/94 AW  build_goto works -- builds a tree with outputs
		labelled at the appropriate nodes.  (Algorithm 2)
  		Added new link called next.  Now all the states
		are linked in chronological order by this link.
   07/31/94 AW	build_fail works -- builds the failure function
		into the tree from build_goto.  (Algorithm 3)
		The failure and output field still needs work.
   08/06/94 AW  all of build_fail works.
   08/23/94 AW  moved to ~aswu/rr directory and modified to work
		with our GA program.  Works ok on one test file.
   08/26/94 AW  Distinguish between keyword_len and len_w_tags.
		The first is the actual length of the basic bb.
		The second is the length of the bb plus all the
	 	associated tags.  This is the number of xters that
		are read in from the fbb file.  We need to make
		this distinction because the exponential coefficient
		runs use the length of the bb to calculate fitness.
   09.18.98.AW	Modified to work generally on list of patterns.  Patterns
		are xter strings.  There is no longer a keyword_len or
		len_w_tags, length has replaced both of those.
   10.09.98.AW	Added global variable called Ascii_offset.  When this was
		only working on binary bases, I could hard code the
		conversion -- the ascii value for numbers is the value
		of the number plus 48 -- '0' = 48, '1' = 49, etc.
		For small letters of the alphabet, the offset is 97, that
		is 'a' = 97, 'b' = 98, 'c' = 99, etc.  So in the cases
		when we have alphabet bases and an array of pointers where
		the element 0 refers to a, element 1 refers to b, etc.,
		I will need to do a different conversion.  So now the 
		conversion factor is set in Ascii_offset in init_run() in
		init.c, and Ascii_offset is used as the conversion variable
		instead of the hard coded value of 48, which would only
		work for decimal numbers.
		***** I make the assumption then that any arrays used with
		alphabet bases puts a in element 0, b in element 1, and so on.
		***** important note -- Ascii_offset only tells how to
		convert bases to number, it doesn't have anything to do
		with the building block names.
*/

#include <stdio.h>
#include <stdlib.h>
#include "types.h"
#include "extern.h"
#include "build.h"

#define BUILD 1
#define ENTER 1
#define FAIL  1

#undef BUILD
#undef ENTER
#undef FAIL

/********** build_goto **********/
/* parameters:	fp		pointer to file -- where to read 
				in words that are to be built
				into a finite state machine(FSM)
 		num_keywords	# of words to be recognized
		context
   called by:
   actions:	reads the words from the input file and calls
		enter_word to add each word to the FSM tree.
*/
struct state_type *build_goto(int num_keywords, GA_WORD *list)
   {
   struct state_type *state;
   struct state_type *ordered;	/* pts to last state added */
   int state_ctr;		/* keeps track of current # of states */
   int bb_name;			/* which bb current word is */
   int i;
   
#ifdef DEBUG
   printf(" ---in build_goto---\n");
#endif

  /* initialize */
   state_ctr = -1;
   state = make_new_state(&state_ctr, num_keywords);
   ordered = state;
#ifdef BUILD
   printf(" state->number = %d state->num_outputs = %d\n", 
          state->number, state->num_outputs);
   if (state->f == NULL) printf(" state->f = NULL\n");
   else                  printf(" state->f = nonNULL\n");

   for (i=0; i<Alphabet_size; i++)
      {
      if (state->g[i] == NULL) printf(" state->g[%d] = NULL\n", i);
      else                     printf(" state->g[%d] = nonNULL\n", i);
      }  /* for i */
   printf(" Output: %d ", state->output[0]);
   printf(" state_ctr = %d\n", state_ctr);
#endif

  /* read in each word.  Format: length string bbname */
  /* enter word into tree */
   for (i=0; i<num_keywords; i++)
      {
      enter_word(list[i].length, list[i].word, list[i].name, &state_ctr, 
                 state, num_keywords, &ordered);
      }  /* for i */

  /* from state 0, point all letters that lead to failure back
     to state 0 */
   for (i=0; i<Alphabet_size; i++)
      {
      if (state->g[i] == NULL)
         {
         state->g[i] = state;
         }  /* if */
      }  /* for i */

#ifdef DEBUG
   printf(" ---end build_goto---\n");
#endif
   return(state);
   }  /* build_goto */

/********** enter_word **********/
/* parameters:	length		length of keyword 
 		keyword		keyword
		bb_name		name of building block to be added
		state_ctr	keeps count of # states so far
				(will be changed)
		state		ptr to state 0
  		num_keywords	total # that need to be recognized
		ordered		ptr to last state added
   called by:	build_goto()
   actions:	enters one new keyword into the FSM.  Starting from
		state 0 goes as far into the tree as possible before
		adding new states.
*/
void enter_word(int length, char *keyword, int bb_name, int *state_ctr,
                struct state_type *state, int num_keywords, 
                struct state_type **ordered)
   {
   struct state_type *new_state;
   struct state_type *current_state;	/* ptr to current state of FSM */
   int xter_ctr;			/* pts to current xter to insert into FSM */
   int i;

#ifdef DEBUG
   printf(" ---in enter_word---\n");
#endif
#ifdef BUILD
   printf(" ordered->number = %d\n", (*ordered)->number);
#endif

  /* initialize */
   current_state = state;
   xter_ctr = 0;

#ifdef JUNK
printf(" current_state->number = %d\n", current_state->number);
printf(" keyword = %s\n", keyword);
printf(" keyword[%d] = %s\n", xter_ctr, &keyword[xter_ctr]);
printf("  (integer)  = %d\n", atoi(&keyword[xter_ctr]));
printf("  (ascii integer)  = %d\n", keyword[xter_ctr]);
#endif

  /* Since our alphabet is simply the set {'0', '1'}, to index
     into the g and f arrays, we can just subtract the ascii 
     value for '0' which is 48.  This will give us the correct
     integer value, where g[0] is the goto function for '0' and
     g[1] is the goto function for '1'.  
     It would probably run even faster if we changed our alphabet
     to the integers {0, 1} but leaving it as characters makes
     the entire GA program more flexible.  */
  /* Ascii_offset is set in init_run(), init.c 10.09.98.AW */

  /* traverse as much of the existing tree as possible with this word */
   while((current_state->g[keyword[xter_ctr]-Ascii_offset] != NULL)
         && (xter_ctr < length))
      {
#ifdef ENTER
printf(" move to state: state %d -> [%d] -> state %d\n", 
       current_state->number, keyword[xter_ctr]-Ascii_offset, 
       current_state->g[keyword[xter_ctr]-Ascii_offset]->number);
#endif
      current_state = current_state->g[keyword[xter_ctr]-Ascii_offset];
      xter_ctr++;
      }  /* while */

  /* add new states to tree to finish off rest of word */
   for (i=xter_ctr; i<length; i++)
      {
      new_state = make_new_state(state_ctr, num_keywords);
      current_state->g[keyword[i]-Ascii_offset] = new_state;
      (*ordered)->next = new_state;
#ifdef ENTER
      printf(" add new state: state %d -> [%d] -> state %d\n",
             current_state->number, keyword[i]-Ascii_offset, new_state->number);
      printf(" ordered: %d -> %d\n", (*ordered)->number, 
             (*ordered)->next->number);
#endif
      current_state = new_state;
      *ordered = (*ordered)->next;
      }  /* for i */
   current_state->output[0] = bb_name;
   current_state->length[0] = length;

#ifdef DEBUG
   printf(" ---end enter_word---\n");
#endif
   }  /* enter_word */

/********** make_new_state **********/
/* parameters:	state_ctr	number of last state created
		num_keywords	in vocabulary
   called by:	enter_word()
   actions:	allocates space for a new state and initializes
		the values.  Increments state_ctr.
*/
struct state_type *make_new_state(int *state_ctr, int num_keywords)
   {
   struct state_type *state;
   int i;

#ifdef DEBUG1
   printf(" ---in make_new_state---\n");
#endif

   state = (struct state_type *)malloc(sizeof(struct state_type));
   (*state_ctr)++;
   state->number = *state_ctr;
   state->next = NULL;
   state->fnext = NULL;
   state->f = NULL;
   state->num_outputs = 1;
   state->output = (int *)malloc(sizeof(int));
   state->length = (int *)malloc(sizeof(int));
   state->output[0] = -1;
   state->length[0] = -1;
   state->g = (struct state_type **)
              malloc(Alphabet_size * sizeof(struct state_type));
   for (i=0; i<Alphabet_size; i++)
      {
      state->g[i] = NULL;
      }  /* for i */

#ifdef DEBUG1
   printf(" ---end make_new_state---\n");
#endif
   return(state);
   }  /* make_new_state */

/********** build_fail **********/
/* parameters:	state		tree of states from build_goto
   called by:
   actions:	adds the failure pointers to the states
		in the tree built by build_goto.  These
		pointers dtm where to go when a xter causes
		a failure condition.
*/
void build_fail(struct state_type *state)
   {
   struct state_type *queue_first;
   struct state_type *queue_last;
   struct state_type *current;		/* r from paper */
   struct state_type *fstate;		/* state from paper */
   int i, j;
   int temp;

#ifdef DEBUG
   printf(" ---in build_fail---\n");
#endif

#ifdef FAIL
  printf(" putting first level states into queue\n");
#endif

  /* set failure function for state 0 to state 0 */
   state->f = state;   

  /* initially the queue is empty, first element must be 
     added specially */
   j=0;
   while (state->g[j] == state)  j++;
   queue_first = state->g[j];
   queue_last  = state->g[j];
   state->g[j]->f = state;

  /* start at state 0 */
   for (i=j+1; i<Alphabet_size; i++)
      {
      if (state->g[i] != state)
         {
        /* add to tail of queue, state->g[i] is s from paper */
         queue_last->fnext = state->g[i];
         queue_last = queue_last->fnext;
         state->g[i]->f = state;
         }  /* if */
      }  /* for i */
#ifdef FAIL
   printf(" Queue:\n");
   print_queue(queue_first);
#endif

  /* recursively process rest of the states, one level at a time */
   while (queue_first != NULL)
      {
     /* take next state off head of queue */
      current = queue_first;
     /* put this at the end of the while loop:
        Leaving it here doesn't work.  Adding a new state to an
        empty list requires special code, as seen up above where
        we add the first state.  If we delete the current state
        before adding any of its children, we'll need to check 
        to see when the queue is empty and when it's not before
        adding the children.  If we add the children first, and
        delete the current state at the end of the while loop, 
        the algorithm is still the same, and we won't need to
        check for and empty queue.  The queue will never be empty
        because, at the very least, the current state will be in it. 
        When the queue does empty out, we should be done. */

/*      queue_first = queue_first->fnext; */

      for (i=0; i<Alphabet_size; i++)
         {
        /* current->g[i] is s from paper */
         if (current->g[i] != NULL)
            {
           /* add to tail of queue */
            queue_last->fnext = current->g[i]; 
            queue_last = queue_last->fnext;

            fstate = current->f;
            while (fstate->g[i] == NULL)
               {
               fstate = fstate->f;
               }  
            current->g[i]->f = fstate->g[i];

           /* copy outputs of current->g[i]->f to 
            current->g[i]->output */

            }  /* if */
         }  /* for i */

      queue_first = queue_first->fnext;
      }  /* while */

#ifdef DEBUG
   printf(" ---end build_fail---\n");
#endif
   }  /* build_fail */
