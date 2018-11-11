/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* extern.h
   06.18.98.AW	Created.
*/

/********** values read in from params files **********/
extern int Rerun;
extern char *Run_num_file;
extern char *Output_path;
extern char *Base;
extern int Max_num_gen;
extern int Loop_until_find_opt;
extern double Min_pct_opt;
extern int Pop_size;
extern int Variable_gen_len;
extern int Hx_window;
extern int Parsimony_pressure;
extern int Max_gen_len;
extern int Min_gen_len;
extern char *Xover_type;
extern double Xover_rate;
extern char *Mut_type;
extern double Mut_rate;
extern double Uniform_x;
extern double Pct_breeding;
extern double Pct_bred;
extern char *Parent_selection;
extern int Parent_replacement_on;
extern int Sigma_scaling_on;
extern double Sigma_scale_min;
extern double Sigma_scale_max;
extern int Tournament_size;
extern double Flat_fitness;
extern int Init_pop;
extern char *Init_pop_file;
extern int Print_params;
extern int Print_function;
extern int Print_pop;
extern int Print_best;
extern int Print_stats;
extern int Print_fxn_best;
extern int Scientific_notation;

/********** values read in from  elsewhere **********/
extern int Run_num;
extern int Max_num_output_files;
extern long Seed;

/********** values calculated **********/
extern int Effective_mut_rate;
extern int Num_breeding;
extern int Num_bred;
extern int Alphabet_size;
extern int Ascii_offset;

/********** other necessary run variables **********/
extern POPULATION Pop;
extern POPULATION Kids;
extern POPULATION Parents;
extern GENERATION Gen;
extern GENERATION Prev_gen;
extern INDIVIDUAL *Run_best_indv;

/********** array of output files **********/
extern OUTPUT_FILE *Output_file;

/********** parse tree for floating representations **********/
extern struct state_type *Tree;

/********** general function stuff **********/
extern char *Function_name;
extern int Found_optimum;

/********** specific function stuff **********/
extern RR rr;
extern TSP tsp;
extern TENSOR tensor; 
