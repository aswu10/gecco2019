/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* output.c
   06.21.98.AW	Created.  Contains routines that output to screen or file.

   Routines:	print_params()		06.21.98.AW
		print_opfiles()		07.08.98.AW
		print_genome()		07.08.98.AW
		fprint_genome()		07.08.98.AW
		print_population()	07.08.98.AW
		print_pop()		07.27.98.AW
		print_gen_best()	07.29.98.AW
		fprint_gen_best()	07.29.98.AW
		fprint_individual()	07.29.98.AW
		gen_output()		07.29.98.AW
		run_output()		07.29.98.AW
*/

#include <stdio.h>
#include <string.h>
#include "types.h"
#include "extern.h"
#include "output.h"
#include "params.h"
#include "indv.h"
#include "fxn.h"

/********** print_params **********/
/* parameters:	fp	where to print, includes stdout (screen)
			Assume that fp has already been fopened.
   called by:	read_params(), params.c
		ga_init(), ga.c
   actions:	prints out the values of the parameters for a particular
		run.  Values to be printed (currently) must be defined in
		global.h.
*/
void print_params(FILE *fp)
   {
  /* things that were read in */
   fprintf(fp, " Run_num = %ld\n", Run_num);
   fprintf(fp, " Seed = %ld\n", Seed);
   fprintf(fp, " Function_name = %s\n", Function_name);
   fprintf(fp, " Rerun = %d\n", Rerun);
   fprintf(fp, " Run_num_file = %s\n", Run_num_file);
   fprintf(fp, " Output_path = %s\n", Output_path);
   fprintf(fp, " Base = %s\n", Base);
   fprintf(fp, " Xters = %s\n", Xters);
   fprintf(fp, " Max_num_gen = %d\n", Max_num_gen);
   fprintf(fp, " Loop_until_find_opt = %d\n", Loop_until_find_opt);
   fprintf(fp, " Min_pct_opt = %lf\n", Min_pct_opt);
   fprintf(fp, " Pop_size = %d\n", Pop_size);
   fprintf(fp, " Variable_gen_len = %d\n", Variable_gen_len);
   fprintf(fp, " Hx_window = %d\n", Hx_window);
   fprintf(fp, " Parsimony_pressure = %d\n", Parsimony_pressure);
   fprintf(fp, " Max_gen_len = %d\n", Max_gen_len);
   fprintf(fp, " Min_gen_len = %d\n", Min_gen_len);
   fprintf(fp, " Init_max_gen_len = %d\n", Init_max_gen_len);
   fprintf(fp, " Init_min_gen_len = %d\n", Init_min_gen_len);
   fprintf(fp, " Xover_type = %s\n", Xover_type);
   fprintf(fp, " Xover_rate = %lf\n", Xover_rate);
   fprintf(fp, " Mut_rate = %lf\n", Mut_rate);
   fprintf(fp, " Uniform_x = %lf\n", Uniform_x);
   fprintf(fp, " Pct_breeding = %lf\n", Pct_breeding);
   fprintf(fp, " Pct_bred = %lf\n", Pct_bred);
   fprintf(fp, " Parent_selection = %s\n", Parent_selection);
   fprintf(fp, " Parent_replacement_on = %d\n", Parent_replacement_on);
   fprintf(fp, " Sigma_scaling_on = %d\n", Sigma_scaling_on);
   fprintf(fp, " Sigma_scale_min = %lf\n", Sigma_scale_min);
   fprintf(fp, " Sigma_scale_max = %lf\n", Sigma_scale_max);
   fprintf(fp, " Tournament_size = %d\n", Tournament_size);
   fprintf(fp, " Flat_fitness = %lf\n", Flat_fitness);
   fprintf(fp, " Init_pop = %d\n", Init_pop);
   fprintf(fp, " Init_pop_file = %s\n", Init_pop_file);
   fprintf(fp, " Print_params = %d\n", Print_params);
   fprintf(fp, " Print_function = %d\n", Print_function);
   fprintf(fp, " Print_pop = %d\n", Print_pop);
   fprintf(fp, " Print_best = %d\n", Print_best);
   fprintf(fp, " Print_stats = %d\n", Print_stats);
   fprintf(fp, " Print_fxn_best = %d\n", Print_fxn_best);
   fprintf(fp, " Print_targets = %d\n", Print_targets);
   fprintf(fp, " Trace = %d\n", Trace);
   fprintf(fp, " CoG = %d\n", CoG);
  /* calculated values */
   fprintf(fp, " Num_breeding = %d\n", Num_breeding);
   fprintf(fp, " Num_bred = %d\n", Num_bred);
   fprintf(fp, " Alphabet_size = %d\n", Alphabet_size);
   fprintf(fp, " Ascii_offset = %d\n", Ascii_offset);
   }  /* print_params */

/********** print_opfiles **********/
/* parameters:	fp	where to print, includes stdout (screen)
			Assume that fp has already been fopened.
   called by:	read_default_opfiles, read_opfiles, params.c
   actions:	print out status of output files -- which ones are
		to be printed and which ones are not.
*/
void print_opfiles(FILE *fp)
   {
   int i;

   fprintf(fp, " Max_num_output_files = %d\n", Max_num_output_files);

   for (i=0; i<Max_num_output_files; i++)
      {
      fprintf(fp, "     %d: %s %d", i, Output_file[i].extension,
			Output_file[i].on);

/*
      if (Output_file[i].on)  fprintf(fp, " %s\n", Output_file[i].filename);
      else  fprintf(fp, "\n");
*/
      fprintf(fp, "\n");
      }  /* for i */
   }  /* print_opfiles */

/********** print_genome **********/
/* parameters:	genome		to print
		length		of genome
		endofline	if = 1, print eoln xter, if = 0 don't
   called by:	
   actions:	
*/
void print_genome(int *genome, int length, int endofline)
   {
   int i;

   for (i=0; i<length; i++)
	  if (!strcmp(Base, "integers"))
	     {
		 printf("%d ", genome[i]);
		 }
      else
	     {
	     putchar(genome[i]); 
		 }
   if (endofline)  putchar('\n');
   }  /* print_genome */

/********** fprint_genome **********/
/* parameters:  fp		where to print
		genome          to print
                length          of genome
                endofline       if = 1, print eoln xter, if = 0 don't
   called by:   
   actions:     
*/
void fprint_genome(FILE *fp, int *genome, int length, int endofline)
   {
   int i;

   for (i=0; i<length; i++) {
      if (!strcmp(Base, "integers"))
	     {
		 fprintf(fp, "%d ", genome[i]);
		 }
      else
	     {
	     putc(genome[i], fp); 
		 }
   }
   if (endofline)  putc('\n', fp);
   }  /* fprint_genome */

/********** print_population **********/
/* parameters:  pop		to print
		first		indv to print
		last		indv to print
   called by:
   actions:	print individuals from a population.
		from and including <first> to <last> individuals.
*/
void print_population(POPULATION pop, int first, int last)
   {
   int i;

   for (i=first; i<=last; i++)
      {
      printf(" %3d %3d %l6.2f ", pop[i]->index,pop[i]->length,pop[i]->fitness);
      printf("%6.2lf ", pop[i]->calc_num_offspring);
      print_genome(pop[i]->genome, pop[i]->length, 1);
      }  /* for i */
   }  /* print_population */

void print_pop(POPULATION pop, int first, int last)
   {
   int i;
   int j;
 
   for (i=first; i<=last; i++)
      {
      printf(" %3d %3d %l6.2f ", pop[i]->index,pop[i]->length,pop[i]->fitness);
/*
      printf("%6.2lf ", pop[i]->calc_num_offspring);
*/
      printf("(%2d %2d) ", pop[i]->parent1_index, pop[i]->parent2_index);
      print_genome(pop[i]->genome, pop[i]->length, 0);
      if (Trace)
         {
         printf(" x%d ", pop[i]->num_xover);
         for (j=0;j<pop[i]->num_xover;j++) printf("%d ", pop[i]->xover_pts[j]);
         printf(" m%d ", pop[i]->num_mut);
         for (j=0;j<pop[i]->num_mut;j++) printf("%d ", pop[i]->mut_pts[j]);
         }  /* if Trace */
      printf("\n");
      }  /* for i */
   }  /* print_pop*/

/********** fprint_population **********/
/* parameters:  fp		ptr to file to print to
		pop		to print
		first		indv to print
		last		indv to print
   called by:
   actions:	print individuals from a population to a file.
		from and including <first> to <last> individuals.
		Prints length on a line, then prints genotype
		of individual on next line.
*/
void fprint_population(FILE *fp, POPULATION pop, int first, int last)
   {
   int i;

   for (i=first; i<=last; i++)
      {
      fprintf(fp, " %d\n", pop[i]->length);
      fprint_genome(fp, pop[i]->genome, pop[i]->length, 1);
      }  /* for */
   }  /* fprint_population */

/********** fprint_individual **********/
/* parameters:  fp
		indv
   called by:
   actions:     print one line to file about given individual.
*/
void fprint_individual(FILE *fp, INDIVIDUAL *indv, int gen)
   {
   fprintf(fp, " G %3d I %3d L %3d F %8.3lf O %8.3lf (%3d, %3d) ",
		gen, indv->index, indv->length,
		indv->fitness, indv->calc_num_offspring,
		indv->parent1_index, indv->parent2_index);
   fprint_genome(fp, indv->genome, indv->length, 1);
   }  /* fprint_individual */

/********** print_gen_best **********/
/* parameters:
   called by:
   actions:	print on one line: gen, avg fitness, stdev, best fitness
		and best individual.
*/
void print_gen_best()
   {
   printf(" Gen %3d %8.3lf %8.3lf %8.3lf ", Gen.index, Gen.avg_fitness,
		Gen.std_dev, Gen.best_fitness);
   print_genome(Pop[Gen.best_indv_index]->genome,
		Pop[Gen.best_indv_index]->length, 1);
   }  /* print_gen_best */

/********** fprint_gen_best **********/
/* parameters:	fp
   called by:
   actions:     print one line to file: gen, avg fitness, stdev, best fitness
                and best individual.
*/
void fprint_gen_best(FILE *fp)
   {
   fprintf(fp, " Gen %3d %8.3lf %8.3lf %8.3lf ", Gen.index, Gen.avg_fitness,
                Gen.std_dev, Gen.best_fitness);
   fprint_genome(fp, Pop[Gen.best_indv_index]->genome,
                Pop[Gen.best_indv_index]->length, 1);
   }  /* fprint_gen_best */

/********** gen_output **********/
/* parameters:
   called by:   gen_end(), stats.c
   actions:	prints everything that needs to be printed at the end
		of a generation.
*/
void gen_output()
   {
   int array_ptr;
   int i, j;

#ifdef DEBUG
   printf(" ---in gen_output---\n");
#endif

  /* output to screen */
   if (Print_stats)
      printf(" Gen %3d avg %lf stdev %lf best %lf (%d) worst %lf (%d)\n",
		Gen.index, Gen.avg_fitness, Gen.std_dev, Gen.best_fitness,
        	Gen.best_indv_index, Gen.worst_fitness, Gen.worst_indv_index);

   if (Print_fxn_best)
      fxn_fprint_gen_indv(stdout, Gen.best_indv_index);
   else if (Print_best)
      print_gen_best();

  // /* output to files */
   if (file_on("genbest"))
      {
      array_ptr = get_file_pointer("genbest");
      Output_file[array_ptr].fp = fopen(Output_file[array_ptr].filename, "a");
      fprint_gen_best(Output_file[array_ptr].fp);
      fclose(Output_file[array_ptr].fp);
      }  /* if genbest */
   if (file_on("genstats"))
      { 
      array_ptr = get_file_pointer("genstats");
      Output_file[array_ptr].fp = fopen(Output_file[array_ptr].filename, "a");
        fprintf(Output_file[array_ptr].fp,
	      " %4d av %6.3lf sd %6.3lf best %6.3lf %6.3lf %3d worst %6.3lf %6.3lf %3d\n",
		Gen.index, Gen.avg_fitness, Gen.std_dev,
		Gen.best_fitness, Pop[Gen.best_indv_index]->raw_fitness,
		Gen.best_indv_index,
		Gen.worst_fitness, Pop[Gen.worst_indv_index]->raw_fitness,
		Gen.worst_indv_index);
      fclose(Output_file[array_ptr].fp);
      }  /* if genstats */
   if (file_on("lenstats"))
      {
      array_ptr = get_file_pointer("lenstats");
      Output_file[array_ptr].fp = fopen(Output_file[array_ptr].filename, "a");
      fprintf(Output_file[array_ptr].fp,
		" %4d avg %7.2lf stdev %6.2lf long %4d %4d %5.2lf %5.2lf short %4d %4d %5.2lf %5.2lf\n",
		Gen.index, Gen.len_avg, Gen.len_std,
		Gen.longest, Gen.longest_index,
		Pop[Gen.longest_index]->fitness,
		Pop[Gen.longest_index]->raw_fitness,
		Gen.shortest, Gen.shortest_index,
		Pop[Gen.shortest_index]->fitness,
		Pop[Gen.shortest_index]->raw_fitness);
      fclose(Output_file[array_ptr].fp);
      }  /* if lenstats */
   if (file_on("popstats"))
      {
      array_ptr = get_file_pointer("popstats");
      Output_file[array_ptr].fp = fopen(Output_file[array_ptr].filename, "a");
      fprintf(Output_file[array_ptr].fp, "Gen %d\n", Gen.index);
      for (i = 0; i < Pop_size; i++)
         {
         fprintf(Output_file[array_ptr].fp, "%d: ", i+1);
         fprint_genome(Output_file[array_ptr].fp, Pop[i]->genome, Pop[i]->length, 1);
         }
      fclose(Output_file[array_ptr].fp);
      } /* if popstats */
// /***** 030802.AW  Added condition that Base is multichar because
       // this only is used when using the PGA representation.  Can
       // make it usable for all Bases by going into pop.c and making
       // sure that the indv_xter_count field is allocated for all
       // Bases.  Field is attached to Pop and Kids. *****/
   // if (file_on("repstats") && (!strcmp(Base, "multichar")))
      // {
      // array_ptr = get_file_pointer("repstats");
      // Output_file[array_ptr].fp = fopen(Output_file[array_ptr].filename, "a");
      // fprintf(Output_file[array_ptr].fp, " %4d", Gen.index);
      // for (i=0; i<Alphabet_size; i++)
         // {
         // fprintf(Output_file[array_ptr].fp, " ");
         // putc(Xters[i], Output_file[array_ptr].fp);
         // fprintf(Output_file[array_ptr].fp, " %lf",
                 // Gen.avg_xter_proportion[i]);
         // }  /* for i */
      // fprintf(Output_file[array_ptr].fp, "\n");
      // fclose(Output_file[array_ptr].fp);
      // }  /* if repstats */

  // /***** 031208AW xtrack file prints avg and sd distance between "adjacent"
   // pairs of identical xters for best indv and avgd over population *****/
   // if (file_on("xtrack") && (!strcmp(Base, "multichar")))
      // {
      // array_ptr = get_file_pointer("xtrack");
      // Output_file[array_ptr].fp = fopen(Output_file[array_ptr].filename, "a");
      // fprintf(Output_file[array_ptr].fp,
             // " %4d avg avg %7.2lf stdev %6.2lf max %d min %d best avg %7.2lf stdev %6.2lf max %d min %d \n",
                // Gen.index,
		// Gen.avg_xtrack->dist_avg,
		// Gen.avg_xtrack->dist_sd,
		// Gen.avg_xtrack->dist_max,
		// Gen.avg_xtrack->dist_min,
		// Gen.best_xtrack->dist_avg,
		// Gen.best_xtrack->dist_sd,
		// Gen.best_xtrack->dist_max,
		// Gen.best_xtrack->dist_min
              // );
      // fclose(Output_file[array_ptr].fp);
      // }  /* if xtrack */

  // /***** 031227AW xcount file prints avg and sd of the number of each xter
   // on an individual in the population.  Also prints number of each xter
   // on best indv of each pop *****/
   // if (file_on("xcount") && (!strcmp(Base, "multichar")))
      // {
      // array_ptr = get_file_pointer("xcount");
      // Output_file[array_ptr].fp = fopen(Output_file[array_ptr].filename, "a");
      // fprintf(Output_file[array_ptr].fp, " %4d ", Gen.index);

      // fprintf(Output_file[array_ptr].fp, "  avg  ");

      // for (i=0; i<Alphabet_size; i++)
         // {
         // putc(Xters[i], match.opfile_values.fp);
         // fprintf(match.opfile_values.fp, " %lf  ", Gen.avg_xter_count[i]);
         // }  /* for i */

      // fprintf(Output_file[array_ptr].fp, "  best  ");

      // for (i=0; i<Alphabet_size; i++)
         // {
         // putc(Xters[i], match.opfile_values.fp);
         // fprintf(match.opfile_values.fp, " %4d  ", 
		// Pop[Gen.best_indv_index]->xter_count[i]);
         // }  /* for i */

      // fprintf(Output_file[array_ptr].fp, "\n");
      // fclose(Output_file[array_ptr].fp);
      // }  /* if xcount */

  // /***** 031227AW xprop file prints avg and sd of the proportion of each xter
   // on an individual in the population.  Also prints the proportion of each xter
   // on best indv of each pop *****/
   // if (file_on("xprop") && (!strcmp(Base, "multichar")))
      // {
      // array_ptr = get_file_pointer("xprop");
      // Output_file[array_ptr].fp = fopen(Output_file[array_ptr].filename, "a");
      // fprintf(Output_file[array_ptr].fp, " %4d ", Gen.index);
         
      // fprintf(Output_file[array_ptr].fp, "  avg  ");
         
      // for (i=0; i<Alphabet_size; i++)
         // {
         // putc(Xters[i], match.opfile_values.fp);
         // fprintf(match.opfile_values.fp, " %lf  ", Gen.avg_xter_proportion[i]);
         // }  /* for i */
   
      // fprintf(Output_file[array_ptr].fp, "  best  ");
      
      // for (i=0; i<Alphabet_size; i++)
         // {
         // putc(Xters[i], match.opfile_values.fp);
         // fprintf(match.opfile_values.fp, " %lf  ",
            // (double)Pop[Gen.best_indv_index]->xter_count[i]/
            // (double)Pop[Gen.best_indv_index]->length
            // );
         // }  /* for i */
                
      // fprintf(Output_file[array_ptr].fp, "\n");
      // fclose(Output_file[array_ptr].fp);
      // }  /* if xprop */



  /* function specific output to files */
   fxn_gen_output();

#ifdef DEBUG
   printf(" ---end gen_output---\n");
#endif
   }  /* gen_output */

/********** run_output **********/
/* parameters:
   called by:   gen_end(), stats.c
   actions:     prints everything that needs to be printed at the end
                of a run.
*/
void run_output()
   {
#ifdef DEBUG
   printf(" ---in run_output---\n");
#endif
   if (file_on("runbest"))
      {
      Output_file[get_file_pointer("runbest")].fp = fopen(
		Output_file[get_file_pointer("runbest")].filename, "a");
      fprint_indv(Output_file[get_file_pointer("runbest")].fp, Run_best_indv);
      fclose(Output_file[get_file_pointer("runbest")].fp);
      }  /* if runbest */
   if (file_on("finalbest"))
      {
      Output_file[get_file_pointer("finalbest")].fp = fopen(
		Output_file[get_file_pointer("finalbest")].filename, "a");
      fprint_indv(Output_file[get_file_pointer("finalbest")].fp, 
			Pop[Gen.best_indv_index]);
      fclose(Output_file[get_file_pointer("finalbest")].fp);
      }  /* if finalbest */
   if (file_on("beststats"))
      {
      Output_file[get_file_pointer("beststats")].fp = fopen(
		Output_file[get_file_pointer("beststats")].filename, "a");
      fprintf(Output_file[get_file_pointer("beststats")].fp, 
             " run %3ld indv %3d gen %3d len %4d fit %lf rawfit %lf\n",
             Run_num,
             Run_best_indv->index,
             Run_best_indv->gen,
             Run_best_indv->length,
             Run_best_indv->fitness,
             Run_best_indv->raw_fitness
             );
      fclose(Output_file[get_file_pointer("beststats")].fp);
      }  /* if beststats */
#ifdef DEBUG
   printf(" ---end run_output---\n");
#endif
   }  /* run_output */
