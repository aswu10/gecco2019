/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* 06.21.98 AW  Created.
   07.23.17 RN: Modified print_genome & fprint_genome to handle
                random keys representation
*/

/* prototypes */
void print_params(FILE *fp);
void print_opfiles(FILE *fp);
void print_genome(INDIVIDUAL *indv, int endofline);
void fprint_genome(FILE *fp, INDIVIDUAL *indv, int endofline);
void print_population(POPULATION pop, int first, int last);
void print_pop(POPULATION pop, int first, int last);
void fprint_population(FILE *fp, POPULATION pop, int first, int last);
void fprint_individual(FILE *fp, INDIVIDUAL *indv, int gen);
void print_gen_best();
void fprint_gen_best(FILE *fp);
void gen_output();
void run_output();
