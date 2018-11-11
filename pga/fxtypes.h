/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* This file contains structures that have to do with specific
   functions or problems that the GA is trying to solve.
   For now, all function structures and declarations will be included
   in here and the entire file will be included in the types.h file.
   Since the GA will only be compiled to solve one problem at a time,
   structures and declarations for other problems will simply be
   inactive and not allocated.

   09.15.98.AW	Created.
		Royal Road functions.
   09.16.98.AW	Micro air vehicle sensor suite problem.
   010308AW	NASA planetary simulator.
*/

/********** Royal Road functions **********/

typedef struct
   {
   int index;
   int start;
   int length;
   char *optimum;
   double fitness;
   }  BBB;

typedef struct
   {
   int index;
   int num_components;
   int *component;
   double fitness;
   }  UBB;

typedef struct
   {
   char optimum_type[INPUT_LINE_LEN];	/* how to make the optimum indv */
   int num_basic_bb;		/* # basic building blocks */
   BBB *basic_bb;
   int num_upper_bb;		/* # upper level bb */
   UBB *upper_bb;
   double optimum_fitness;
/* moved to INDIVIDUAL */
/*   int *bbb_exists;		/* for evaluating ubb and trace printing */
/*   char *marker;		/* for floating bb only */
   int trace_bb_data;		/* if 1, trace cxn/dxn of bb; if 0, don't */
  /* function specific output files */
   int num_opfiles;
   char **opfile;
   }  RR;

/* for tracing the construction and disruption of bb by mutation and xover */
/* each structure holds one generations worth of data */
typedef struct
   {
   int sum;
   int cxn;
   int dxn;
   int bb_cxn;
   int bb_dxn;
   }  STUPID;

/********** NASA planetary simulator **********/

typedef struct
   {
   double min;			/* range of this input value */
   double max;
   char positive[2];		/* xter that represents */
   char negative[2];		/* xter that represents */
   char name[INPUT_LINE_LEN];
   char abbreviation[INPUT_LINE_LEN];
   int pos_count;
   int neg_count;
   double proportion;
   }  INPARAM;

typedef struct
   {
   int num_inputs;		/* to represent for sim */
   INPARAM *inparam;		/* info about inputs */
   OUTPUT_FILE opfile_simstats;	/* print basic simulator stats */
   }  SIM;

/********** matching floating point numbers *********
6.8.01.IG not only poss and neg chars, now use any base
		up to decimal.
*/

typedef struct
   {
   double val;
   double range_low;
   double range_high;
   double range;
   double pct_in_range;		/* ptr of relative position of value in range */
				/* e.g. range=0 to 5, pct for 1 is 0.2 */
				/* pct for 4.5 is 0.9, pct for 3 is 0.6 */
  /* PGA representation */
   char negative_xter[2];
   char positive_xter[2];
  /* binary representation */
   int start_bit;
   int num_bits;
   }  VALUE;
 
typedef struct
   {
   char name[INPUT_LINE_LEN];	/* the name of this target */
   int num_val;			/* how many values to match */
   VALUE *val;			/* the values */
   int duration;		/* how many generations this target is valid */
   }  TARGET;

typedef struct
   {
   int num_targets;		/* how many targets */
   TARGET *target;
   OUTPUT_FILE opfile_values;	/* print generated parameter values */
   int current;
   int gencount;
   }  MATCH;

/********** resource allocation *********/
typedef struct
   {
   double actual;
   double percent;
  /* for proportional GA */
   char xter[2];
  /* for binary GA */
   int start_bit;
   int length;
   }  COMPONENT;

typedef struct
   {
   int num_comp;
   COMPONENT *comp;
   double sum;
   OUTPUT_FILE opfile_values;
   int random_values;
   }  RA;

   /********** TSP *********/ 
typedef struct
   {
   int x;
   int y;
   } COORDS;

typedef struct 
   {
   int num_cities;
   COORDS **city_coordinates;
   } TSP;
   
   /********** tensor contraction *********/
typedef struct
   {
   double key;
   int value;
   } KEY_MAP;
   
typedef struct
   {
   int id;
   int x;
   int y;
   int z;
   } NODE;
   
typedef struct
   {
   // int id;
   int weight;
   int node_a;
   int node_b;
   int active;
   int redundant;
   // int *redundant_edges;
   } EDGE;

typedef struct
   {
   int length;
   int num_edges;
   EDGE **edges;
   } TENSOR;