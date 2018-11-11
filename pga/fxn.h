/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* fxn.h
   07.03.98.AW	Created.
   030528.AW	Added fxn_gen_start() and fxn_gen_end().
*/

/* prototypes */
int read_fxn_file(char *fxn_file);
int init_function();
void end_function();
void eval_indv(INDIVIDUAL *indv);
int found_solution();
void fprint_fxn(FILE *fp);
void fprint_genes(FILE *fp, INDIVIDUAL *indv);
void fxn_fprint_gen_indv(FILE *fp, int indv);
void fxn_gen_output();
int fxn_gen_start();
int fxn_gen_end();
