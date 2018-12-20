/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* global.h
   06.18.98.AW	Created.
		File of global variables.  Included in ga.c.
		All other .c files should include extern.h instead.
   08.2.17 RN: added Mut_type
*/

/********** run parameters read in **********/
int Rerun;
char *Run_num_file;
char *Output_path;
char *Base;			/* nucleotides: binary or letter */
int Max_num_gen;
int Loop_until_find_opt;
double Min_pct_opt;

int Metric;     /* 0 = distance, 1 = time */
int Pop_size;
int Variable_gen_len;		/* 0 = fixed genome length = Max_gen_len */
				/* 1 = variable, init pop indv sizes random */
				/* >1 = variable = init pop indv gen len */
				/* e.g. if this variable were set to 52, */
				/* that means variable length and initial */
				/* pop of individuals all of length 52. */
int Hx_window;			/* if Variable_gen_len, this is the size */
				/* of the homologous crossover window */
int Parsimony_pressure;		/* if 0 no pressure, else this is length bias */
int Max_gen_len;
int Min_gen_len;
char *Xover_type;		/* fixed len: one-point, two-point, uniform*/
				/* var len: homologous, random (both 1pt) */
double Xover_rate;		/* prob that each pair of parents will xover*/
char *Mut_type;         /* random, displacement, insertion, inversion */
double Mut_rate;
double Uniform_x;		/* uniform xover rate, if Xover_type=uniform */
double Pct_breeding;		/* % pop that become parents */
double Pct_bred;		/* % pop not clones from prev population */
				/* or % pop that were created via genops */
char *Parent_selection;		/* proportional or tournament */
int Parent_replacement_on;	/* indvs can be parents multiple times */
int Sigma_scaling_on;
double Sigma_scale_min;		/* min expected offspring if Sigma_scaling on */
double Sigma_scale_max;		/* max expected offspring if Sigma_scaling on */
int Tournament_size;		/* if tournament selection */

double Gaussian_sd;         /* standard deviation for Gaussian mutation */

double Flat_fitness;		/* if = 0, eval indvs as usual 
				   if > 0, assign all indvs that fitness val */

int Init_pop;			/* what is initial pop?  0 = all 0's, */
				/* 1 = all 1's, 2 = random, 3 = from file */
char *Init_pop_file;		/* if Init_pop=3, read from here */

				/* about printing to screen during run */
int Print_params;		/* if 1, print params at start of run */
int Print_function;		/* if 1, print function at start of run */
int Print_pop;			/* if 1, print population each generation */
int Print_best;			/* if 1, print best indv each generation */
int Print_stats;		/* if 1, print stats for each generation */
int Print_fxn_best;		/* if 1, overrides Print_best with fxn format */
				/* Print_best prints genome, Print_fxn_best
				   prints genome in fxn specific format */
int Scientific_notation; /* use scientific notation when printing fitness */

/********** parameters read in from elsewhere **********/
int Run_num;			/* read in from Run_num_file */
int Max_num_output_files;	/* read in from opfiles.default */
long Seed;			/* random seed, read in or generated */

/********** parameters calculated **********/
int Num_breeding;		/* Pop_size * Pct_breeding, must be even */
int Num_bred;			/* Pop_size * Pct_bred, must be even */
int Alphabet_size;		/* depends on Base */
int Ascii_offset;		/* depends on Base */
				/* tells how to change base to # */
				/* for integers: '0' = 48 -> '0' = 0+48 */
				/* for small letters: 'a' = 97 */

/********** other necessary run variables **********/
POPULATION Pop;
POPULATION Kids;
POPULATION Parents;
GENERATION Gen;
GENERATION Prev_gen;
INDIVIDUAL *Run_best_indv;

/********** array of output files **********/
OUTPUT_FILE *Output_file;

/********** parse tree for floating representations **********/
struct state_type *Tree;

/********** general function stuff **********/
char *Function_name;
int Found_optimum;		/* if fxn has optimum, set to 1 when found */

/********** specific function stuff **********/
RR rr;
TSP tsp;
TENSOR tensor;

