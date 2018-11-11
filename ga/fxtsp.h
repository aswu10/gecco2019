/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* fxn2.h
   07.03.98.AW	Created.
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
void decode(INDIVIDUAL *indv); // 07.23.17: RN random keys

/* for RR functions only */
void trace_bb_data();
