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
   07.27.17.RN	Added TSP structs.
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
 
// Traveling Salesman   
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

// tensor contraction
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
   } EDGE;

typedef struct
   {
   int length;
   int num_edges;
   EDGE **edges;
   } TENSOR;
   
// for edge-recombination crossover   
typedef struct
   {
   int *edges;
   int length;
   }  EDGE_MAP;
   