/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* 07.28.98.AW Created.
*/

/* prototypes */
void copy_indv(INDIVIDUAL *indv1, INDIVIDUAL *indv2);
int init_run_best_indv();
void init_indv(INDIVIDUAL *indv, int gen);
void print_indv(INDIVIDUAL *indv);
void fprint_indv(FILE *fp, INDIVIDUAL *indv);
void free_indv(INDIVIDUAL *indv);
