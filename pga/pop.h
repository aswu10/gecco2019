/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* pop.h
   07.03.98.AW	Created.
*/

/* prototypes */
int init_pop();
int pop_allocate();
void indv_fill(int *genome, int len, int x);
int pop_fill();
void pop_individuals_0();
void pop_individuals_1();
void pop_individuals_random();
int pop_individuals_readin();
void pop_individuals_random_alphabet();
void pop_individuals_random_multichar();
void pop_individuals_random_int();
void pop_eval();
void close_pop();
