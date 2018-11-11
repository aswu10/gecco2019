/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* pop.h
   07.03.98.AW	Created.
   07.19.17 RN: Commented out for integer conversion
   07.23.17 RN: added pop_individuals_rkeys()
*/

/* prototypes */
int init_pop();
int pop_allocate();
int pop_fill();
void pop_individuals_random();
void pop_individuals_rkeys();
void pop_eval();
void close_pop();
