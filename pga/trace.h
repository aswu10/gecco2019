/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* trace.h */

int init_trace();
int end_trace();
int trace_setup_output_files();
int trace_gen_start();
int trace_gen_end();
void trace_fprint_gen(FILE *fp, int median_fitness);
void trace_fprint_indv(FILE *fp, INDIVIDUAL *indv);
int get_median();
/*
void trace_fprint_genes(FILE *fp, INDIVIDUAL *indv);
*/
