/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* trace.c
   09.22.98.AW	Created
		Routines that have to do with saving the trace
		data for a run.

   Routines:	init_trace()			09.22.98.AW
		end_trace()			09.22.98.AW
		trace_setup_output_files()	09.23.98.AW
		trace_gen_start()		09.23.98.AW
		trace_gen_end()			09.23.98.AW
		trace_fprint_gen()		09.23.98.AW
		trace_fprint_indv()		09.23.98.AW
		trace_fprint_genes()		09.23.98.AW

   12.04.98.AW	All routines in here seem to have to do with printing
		data to the trace files.  Should only be called if Trace = 1
		and not if Trace = 2.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "types.h"
#include "extern.h"
#include "trace.h"
#include "output.h"

/********** init_trace **********/
/* parameters:
   called by:	ga_init(), ga.c
   actions:	initialize stuff for tracing data.
		e.g. make the trace directory and open the file
		pointers
*/
int init_trace()
   {
   char tempstring[LONG_LINE_LEN];
   int error;

#ifdef DEBUG
   printf(" ---in init_trace---\n");
#endif

  /* make the trace directory */
   sprintf(tempstring, "mkdir %s/run.%d/trace", Output_path, Run_num);

#ifdef CMD
   printf(" CMD: %s\n", tempstring);
#endif
   error = system(tempstring);
   if (error != 0)
      printf(" Error(init_trace): with system command: %s\n", tempstring);

  /* setup output files that are going to be written to once per
     generation */
   if (trace_setup_output_files() == ERROR)  return ERROR;

#ifdef DEBUG
   printf(" ---end init_trace---\n");
#endif
   return OK;
   }  /* init_trace */

/********** end_trace **********/
/* parameters:
   called by:
   actions:     finish up stuff for tracing data.
*/
int end_trace()
   {
   int i;
   char tempstring[INPUT_LINE_LEN];
   FILE *fp;

#ifdef DEBUG
   printf(" ---in end_trace---\n");
#endif

  /* print the vis_params file including start and final generation */
   sprintf(tempstring, "%s/run.%d/trace/vis_params", Output_path, Run_num);
   fp = fopen(tempstring, "w");
   fprintf(fp, "0\n%d\n%d\n", Gen.index, Pop_size);

   if (!strcmp(Base,"binary"))  fprintf(fp, "01\n");
   else if (!strcmp(Base,"alphabet"))
      fprintf(fp, "abcdefghijklmnopqrstuvwxyz\n");
   fclose(fp);

  /* deallocate Trace_file space (file should already be closed
     because we are appending and closing for each write */
   for (i=0; i<2; i++)
      {
      free(Trace_file[i].filename);
      }  /* fir */
   free(Trace_file);

#ifdef DEBUG
   printf(" ---end end_trace---\n");
#endif
   return OK;
   }  /* end_trace */

/********** trace_setup_output_files **********/
/* parameters:
   called by:	init_trace()
   actions:	Initialize output files that will be written to
		once per generation so that they can be appended
		to each generation.
   09.23.98.AW	vis_best and vis_median file pointers
*/
int trace_setup_output_files()
   {
   int i;

#ifdef DEBUG
   printf(" ---in trace_setup_output_files---\n");
#endif

  /* for now there are only two files that are written to once
     a generation, the best and median files.  The rest can be
     created on the fly as they are printed. */
  /* allocate space */
   Trace_file = (OUTPUT_FILE *)malloc(2 * sizeof(OUTPUT_FILE));
   if (Trace_file == NULL)
      {
      printf(" Error(trace_setup_output_files): cannot allocate space:");
      printf(" Trace_file\n");
      return ERROR;
      }

  /* allocate space for filenames */
   for (i=0; i<2; i++)
      {
      Trace_file[i].filename = (char *)malloc(INPUT_LINE_LEN * sizeof(char));
      if (Trace_file[i].filename == NULL)
         {
         printf(" Error(trace_setup_output_files): cannot allocate space:");
         printf(" Trace_file[%d].filename\n", i);
         return ERROR;
         }  /* if */
      }  /* for i */

  /* name the files */
   sprintf(Trace_file[0].filename, "%s/run.%d/trace/vis_best",
		Output_path, Run_num);
   sprintf(Trace_file[1].filename, "%s/run.%d/trace/vis_median",
		Output_path, Run_num);

  /* initialize files */
   for (i=0; i<2; i++)
      {
      Trace_file[i].fp = fopen(Trace_file[i].filename, "w");
      if (Trace_file[i].fp == NULL)
         {
         printf(" Error(trace_setup_output_files): cannot open file: %s",
		Trace_file[i].filename);
         return ERROR;
         }  /* if */
      fclose(Trace_file[i].fp);
      }  /* for */

#ifdef DEBUG
   printf(" ---end trace_setup_output_files---\n");
#endif
   return OK;
   }  /* trace_setup_output_files */

/********** trace_gen_start **********/
/* parameters:
   called by:	run_start(), gen_start(), stats.c
   actions:	outputs stuff at end of a generation
*/
int trace_gen_start()
   {
   char tempstring[INPUT_LINE_LEN];
   int error;

#ifdef DEBUG
   printf(" ---in trace_gen_start---\n");
#endif

  /* init to keep track of number of xovers in generation */
  /* count number of xovers -- max = pop_size/2 if 2 parents */
   Gen.num_x = 0;

  /* make directory for this generation */
   sprintf(tempstring, "mkdir %s/run.%d/trace/gen%d/",
                Output_path, Run_num, Gen.index);
#ifdef CMD
   printf(" CMD: %s\n", tempstring);
#endif
   error = system(tempstring);
   if (error != 0)
      {
      printf(" Error(trace_gen_start): with system command: %s\n", tempstring);
      return ERROR;
      }

#ifdef DEBUG
   printf(" ---end trace_gen_start---\n");
#endif
   return OK;
   }  /* trace_gen_start */

/********** trace_gen_end **********/
/* parameters:
   called by:	gen_end(), stats.c
   actions:	outputs stuff at end of a generation
*/
int trace_gen_end()
   {
   FILE *fp;
   char tempstring[INPUT_LINE_LEN];
   int i;
   int median_index;

#ifdef DEBUG
   printf(" ---in trace_gen_end---\n");
#endif
  /* find median */
   median_index = get_median();

  /* print to file genx/vis.genx */
   sprintf(tempstring, "%s/run.%d/trace/gen%d/vis.gen%d",
		Output_path, Run_num, Gen.index, Gen.index);
   fp = fopen(tempstring, "w");
   trace_fprint_gen(fp, median_index);
   fclose(fp);

  /* print best individual to vis.best */
   sprintf(tempstring, "%s/run.%d/trace/vis_best", Output_path, Run_num);
   fp = fopen(tempstring, "w");
   trace_fprint_indv(fp, Pop[Gen.best_indv_index]);
   fclose(fp);

  /* print median individual to vis.median */
   sprintf(tempstring, "%s/run.%d/trace/vis_median", Output_path, Run_num);
   fp = fopen(tempstring, "w");
   trace_fprint_indv(fp, Pop[median_index]);
   fclose(fp);

#ifdef DEBUG
   printf(" ---end trace_gen_end---\n");
#endif
   return OK;
   }  /* trace_gen_end */

/********** trace_fprint_gen **********/
/* parameters:	fp		file to print to
		median_index	median indv in population
   called by:
   actions:     print generation data to file
*/
void trace_fprint_gen(FILE *fp, int median_index)
   {
   fprintf(fp, "generation\n");
   fprintf(fp, "%d\n", Gen.index);
   fprintf(fp, "%lf\n", Gen.best_fitness);
   fprintf(fp, "%d\n", Gen.best_indv_index);
   fprintf(fp, "%lf\n", Gen.avg_fitness);
   fprintf(fp, "%lf\n", Gen.std_dev);
  /* median fitness individual */
   fprintf(fp, "%lf\n", Pop[median_index]->fitness);
   fprintf(fp, "%d\n", median_index);
  /* xover info */
   fprintf(fp, "%d\n", Gen.num_x);
  /* length info */
   fprintf(fp, "%d\n", Gen.longest);
   fprintf(fp, "%d\n", Gen.shortest);
   fprintf(fp, "%lf\n", Gen.len_avg);
   fprintf(fp, "%lf\n", Gen.len_std);
   }  /* trace_fprint_gen */

/********** trace_fprint_indv **********/
/* parameters:
   called by:	trace_gen_end()
		eval_indv(), function files
   actions:	print one individual to file
		all except for genes which will be printed by
		a function specific routine.
*/
void trace_fprint_indv(FILE *fp, INDIVIDUAL *indv)
   {
   int i;

   fprintf(fp, "individual\n");
   fprintf(fp, "%d\n", indv->index);
   fprintf(fp, "%d\n", indv->length);
   fprintf(fp, "%lf\n", indv->fitness);
   fprint_genome(fp, indv->genome, indv->length, 1);
  /* print parents */
   if (indv->parent2_index < 0 && indv->parent1_index < 0)
      {  /* if 0 parents -- should only be true if gen 0 */
      fprintf(fp, "0\n");
      }  /* if */
   else if (indv->parent2_index < 0)
      {  /* if 1 parent -- when no xover */
      fprintf(fp, "1\n");
      fprintf(fp, "%d\n", indv->parent1_index);
      }  /* else */
   else
      {  /* if 2 parents */
      fprintf(fp, "2\n");
      fprintf(fp, "%d\n", indv->parent1_index);
      fprintf(fp, "%d\n", indv->parent2_index);
      }  /* if */
  /* print xovers (locations from parents) */
   fprintf(fp, "%d\n", indv->num_segments);
   for (i=0; i<indv->num_segments; i++)
      {
      fprintf(fp, "%d %d %d\n", indv->seg_parent[i],
		indv->seg_start[i], indv->seg_len[i]);
      }  /* for i */
  /* print mutation locations */
   fprintf(fp, "%d\n", indv->num_mut);
   for (i=0; i<indv->num_mut; i++)
      fprintf(fp, "%d\n", indv->mut_pts[i]);
   }  /* trace_fprint_indv */

#ifdef LEFTOUT
/********** trace_fprint_genes **********/
/* parameters:
   called by:	trace_fprint_indv()
   actions:     print gene data to file
		This may differ from one function to the next.  If so,
		this may become one of those function specific routines.	
   09.23.98.AW	Royal Road function: print number of genes (total number
		of optimum building blocks found including duplicates.
		For each bb print bb name, bb first bit, bb last bit.
*/
void trace_fprint_genes(FILE *fp, INDIVIDUAL *indv)
   {
   }  /* trace_fprint_genes */
#endif  /* LEFTOUT */

/********** get_median **********/
/* parameters:
   called by:	trace_gen_end()
   actions:	returns index to individual in population with median
		fitness value.
   000613AW	Presently a null routine that simply returns the value 0
*/
int get_median()
   {
   return 0;
   }  /* get_median */
