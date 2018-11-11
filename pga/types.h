/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* File containing descriptions of types, typedefs, and structures
   used in the program.
 
   To see what order the program calls subroutines, define DEBUG.  The program
   will now print a message as it enters and leaves each subroutine.

   09.05.93 AW  Created.  Uses a number of declarations from Terry
                Jones' program.
   06.17.98.AW	Recreated.
*/
 
#define DEBUG      1	/* notes when entering and exitting routines */
#undef DEBUG

#define CMD		/* print out system commands before executing */
#undef CMD

#define STEP		/* allows user to step thru GA one gen at a time */
#define SMALLSTEP	/* pauses at parent and kid generations too */
#define SMALLERSTEP	/* pauses after each reproduction event */
#undef STEP
#undef SMALLSTEP
#undef SMALLERSTEP

#define NUMBER_POP 1	/* numbers the population when printing */
#undef NUMBER_POP

/********** error codes **********/
#define OK	 0
#define ERROR	-1

/********** other defines **********/
#define ENDOFFILE 59
#define INPUT_LINE_LEN 100
#define LONG_LINE_LEN 200
#define REALLY_BIG_NUMBER 99999
#define REALLY_SMALL_NUMBER -99999

/********** default values for GA **********/
#define SIGMA_SCALE_MIN	0.0	/* min expected offspring of any indv */
#define SIGMA_SCALE_MAX	1.5	/* max expected offspring of any indv */
/*#define PARAM_FILE "params"	/* file containing default values */
/*#define OUTPUT_FILE "outputs"	/* default settings of which output
				   files are to be printed for a run */

/********** typedefs **********/

/* this is the type of each gene bit */
/* typedef char DNA; */

/* structure for storing output file information.  Includes
   standard output files and additional ones that may be added
   at the start of a run.
*/
typedef struct
   {
   int on;		/* if on=1, print; if on=0, don't print */
   char *extension;
   char *filename;
   FILE *fp;
   }  OUTPUT_FILE;

/* structure for tracking xters and their locations */
/* holds data for a single individual. */
typedef struct
   {
   int dist_sum;	/* sum of distances between pairs of id xters */
   double dist_avg;	/* avg distance between pairs of id xters */
   double dist_sd;	/* avg distance between pairs of id xters */
   int dist_min;	/* shortest distance between pairs of id xters */
   int dist_max;	/* longest distance between pairs of id xters */
  /* need to keep the next two values in order to calc avg over entire pop */
   int num;		/* number of pairs calculated for this indv */
   int dist_sum2;	/* sum of distances squared */
   }  XTER_TRACK;

/* structure to store data on each individual. */
typedef struct
   {
   int index;		/* to help identify individuals after selection */
   int gen;		/* only used when storing best indv of run */
   int length;
   int *genome;
   double fitness;
   double raw_fitness;	/* if parsimony pressure, fitness without pp */
   double calc_num_offspring;
   int parent1_index;
   int parent2_index;
   int chosen;	/* for use when selecting without replacement */
   int can_mate;	/* can the indv be chosen to be a parent, x&m */
  /* 000401AW:  count # of each xter */
   int *xter_count;
  /* only active when trace file is turned on */
   int num_xover;
   int *xover_pts;
   int num_mut;
   int *mut_pts;
   int num_kids;
   int *kids;
   int num_segments;    /* number of segments = num_xover pts + 1 */
   int *seg_parent;     /* which parent contributed segment */
   int *seg_start;      /* first bit of segment on parent, count from 0 */
   int *seg_len;        /* length of contributed segment */
  /* for floating RR function only for now 98.12.02 */
  /* allocated and deallodated in fxfrr.c for now */
  /* 990102AW: for any functions bbs, alloc in init_function,
     dealloc in end_function(), if used in eval_indv() */
   int *bbb_exists;	/* keep track of which bb exist on an individual */
			/* for evaluating ubb and trace printing */
   char *marker;	/* mark location of floating bb */
   int *ubb_exists;	/* keep track of which ubb exist on an individual */
  /* 990123AW:  init to same len as bbb_exists or marker, 0 = no bb,
     1 = newly cxn bb, >1 = inherited bb */
   int *dai;		/* how many gen a bb has been alive */
  /* 010412AW: only active when center of gravity (CoG) is on */
   int *cog_count;
   double *cog_avg;
   double *cog_sd;
   double *cog_low;
   double *cog_high;
  /* 031208AW: for tracking statistics about xters and locations */
  /* only allocated if Xter_track is on */
   XTER_TRACK *xtrack;
   
   
   
   double avg_gene_val; // avg of integers in genome
   double avg_val1_loc;
   double avg_val2_loc;
   double avg_val3_loc;
   }  INDIVIDUAL;

/* a population is a list of pointers that point to individuals. */
typedef INDIVIDUAL **POPULATION;

/* structure to store data for a generation */
typedef struct
   {
   int index;
   double fitness_sum;	/* used in reproduce.c ln 101 to do sigma scaling */
   double avg_fitness;
   double std_dev;
   int best_indv_index;
   double best_fitness;
   int worst_indv_index;
   double worst_fitness;
  /* added 10.09.98.AW for length statistics */
   double len_avg;
   double len_std;
   int longest;
   int longest_index;
   int shortest;
   int shortest_index;
  /* number of pairs of parents that crossover in gen */
   int num_x;
  /* 000401AW:  count # of each xter divided by genome length */
   double *avg_xter_proportion;
  /* 031227AW:  count # of each xter */
   double *avg_xter_count;
  /* 031208AW: for tracking statistics about xters and locations */
  /* only allocated if Xter_track is on */
   XTER_TRACK *avg_xtrack;
   XTER_TRACK *best_xtrack;
   }  GENERATION;

/********** declarations for floating representations **********/
/* for building parse tree for finding tags, from Aho & Corasick paper */

/* #define ALPHABET_SIZE 2   /* total number of xters in alphabet */
 
/* structure for information on each state:
   number - the number or name of the state
   num_outputs - the number of outputs at this state.  In other
            words, the number of strings that end at this state.
   output - this should be an array of length num_keywords since
            that is the maximum number of outputs there could be
            at a state.  Each value stored is the name of a bb
            that ends at this state.
   length - one for each output -- this is the length of the
            string for each output.
   level  - one for each output -- the level of the string for
            for each output.
   g      - the goto function for this state.  There should
            be one array element for each xter in the alphabet.
   f      - the failure function.
            pointer of where to go if a state fails.
   next   - pointer to the next state that was created.  This forms
            a one dimensional list of the states in order by number,
            used for printing list information.
   fnext  - pointer used for building the failure function.
*/
struct state_type
   {
   int number;
   int num_outputs;
   int *output;
   int *length;
   struct state_type **g;
   struct state_type *f;
   struct state_type *next;
   struct state_type *fnext;
   };  /* struct state_type */

/* structure for holding list of keywords to be passed into
   build to build the parse tree */
typedef struct
   {
   int length;
   char *word;
   int name;		/* name of this keyword */
   }  GA_WORD;

/********** function related structures and definitions **********/
#include "fxtypes.h"
