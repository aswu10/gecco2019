/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* ga.c
   06.18.98.AW	Created.  Will contain routines for starting ga
		and main ga loop.
   Routines:
	ga_start()	06.18.98.AW
	ga_init()	06.18.98.AW
	ga_loop()	06.18.98.AW
	ga_continue()	06.18.98.AW
        ga_end()	06.19.98.AW
*/

#include <stdio.h>
#include <stdlib.h>
#include "types.h"
#include "global.h"
#include "ga.h"
#include "init.h"
#include "pop.h"
#include "stats.h"
#include "reproduce.h"

/********** ga_start **********/
/* parameters:  params_file     list of run parameters
		opfiles_file	list of output files
		fxn_file	file containing function specs.
   called by:   main(), main.c
   actions:     this is the ga program.
*/
int ga_start(char *params_file, char *opfiles_file, char *fxn_file, double xover_rate, double mut_rate, double uniform_xover_rate, char *xover_op, char *mut_op)
   {
   int error;

#ifdef DEBUG
   printf(" ---in ga_start()---\n");
#endif

   error = ga_init(params_file, opfiles_file, fxn_file, xover_rate, mut_rate, uniform_xover_rate, xover_op, mut_op);
   if (error == ERROR)
      {
      printf(" Error(ga): ga_init ends on error.\n");
      return ERROR;
      }

   error = ga_loop();
   if (error == ERROR)
      {
      printf(" Error(ga): ga_loop ends on error.\n");
      return ERROR;
      }

   error = ga_end();
   if (error == ERROR)
      {
      printf(" Error(ga): ga_end ends on error.\n");
      return ERROR;
      }

#ifdef DEBUG
   printf(" ---end ga_start()---\n");
#endif
   return OK;
   }  /* ga_start */

/********** ga_init **********/
/* parameters:	params_file     list of run parameters
                opfiles_file    list of output files
                fxn_file        file containing function specs.
   called by:   ga_start(), ga.c
   actions:     Initialize a GA run.
*/
int ga_init(char *params_file, char *opfiles_file, char *fxn_file, double xover_rate, double mut_rate, double uniform_xover_rate, char *xover_op, char *mut_op)
   {
#ifdef DEBUG
   printf(" ---in ga_init()---\n");
#endif

  /* allocate space */
   Output_path = (char *)malloc(INPUT_LINE_LEN * sizeof(char));
   Run_num_file = (char *)malloc(INPUT_LINE_LEN * sizeof(char));
   Base = (char *)malloc(INPUT_LINE_LEN * sizeof(char));
   Init_pop_file = (char *)malloc(INPUT_LINE_LEN * sizeof(char));
   Xover_type = (char *)malloc(INPUT_LINE_LEN * sizeof(char));
   Mut_type = (char *)malloc(INPUT_LINE_LEN * sizeof(char));
   Parent_selection = (char *)malloc(INPUT_LINE_LEN * sizeof(char));
   Function_name = (char *)malloc(INPUT_LINE_LEN * sizeof(char));

  /* read default parameter values and desired outputfiles */
   if (read_defaults() == ERROR)  return ERROR;

  /* read overriding parameter values and outputfiles */
   if (read_params(params_file) == ERROR)  return ERROR;
   if (read_opfiles(opfiles_file) == ERROR)  return ERROR;
   
   Xover_type = xover_op;
   Xover_rate = xover_rate;
   Mut_type = mut_op;
   Mut_rate = mut_rate;
   Uniform_x = uniform_xover_rate;
   
   // printf(" --- Press any key to continue ---\n");  fgetc(stdin);

  /* initialize run stuff */
   if (init_run() == ERROR)  return ERROR;
	/* include calculating global var in init_run? */

  /* set up output files */
   if (setup_output_files() == ERROR)  return ERROR;

  /* initialize fitness function -- function specific routine */
   if (read_fxn_file(fxn_file) == ERROR)  return ERROR;
   if (init_function() == ERROR)  return ERROR;

  /* initialize population */
   if (init_pop() == ERROR)  return ERROR;
   
   if (Print_params)
      {
      print_params(stdout);
      printf("\n");
      print_opfiles(stdout);
      printf("\n");
      }

  /* initialize any stats for run */
   if (run_start() == ERROR)  return ERROR;

  /* evaluate initial population */
   pop_eval();
   if (Print_pop)
      {
      printf(" Current population (gen %d):\n", Gen.index);
      print_pop(Pop, 0, Pop_size-1);
      }  /* if */

  /* calculate gen stats */
   if (gen_stats() == ERROR)  return ERROR;

#ifdef STEP
      printf(" Generation %d %d\n", Gen.index);
      print_population(Pop, 0, Pop_size-1);
      printf(" --- Press any key to continue ---\n");  fgetc(stdin);
#endif

   /* make offspring and put in Kids array */
   if (reproduce() == ERROR)  return ERROR;

  /* end of gen calculations */
   if (gen_end() == ERROR)  return ERROR;

#ifdef DEBUG
   printf(" ---end ga_init()---\n");
#endif
   return OK;
   }  /* ga_init */

/********** ga_loop **********/
/* parameters:
   called by:	ga_start(), ga.c
   actions:	This is the main ga loop.
*/
int ga_loop()
   {
   POPULATION temp;

#ifdef DEBUG
   printf(" ---in ga_loop()---\n");
#endif

   while (ga_continue())
      {
     /* first thing to do is to make the offspring the current generation */
      temp = Pop;
      Pop = Kids;
      Kids = temp;

     /* start of gen stats and calculations & increment generation counter */
      if (gen_start() == ERROR)  return ERROR;

      pop_eval();

     /* calculate gen stats relating to fitness values */
      if (gen_stats() == ERROR)  return ERROR;

     /* FXN_SPECIFIC_CALL */
     /* if we are tracing the construction and disruption of
        building blocks, this is where to do it */
     /* at this point, Pop contains the current evaluated generation
        and Kids actually contains the evaluated parent generation */
     /* 12.30.98.AW check call to trace_bb_data()  This will only
        compile if all functions have this file.  Put a dummy function
        into fxn.c (onemax problem) for now. */
      if (!strcmp(Function_name, "FRR"))
         {
         if (rr.trace_bb_data == 1)  trace_bb_data();
         }
      else if (!strcmp(Function_name, "RR"))
         {
         if (rr.trace_bb_data == 1)  trace_bb_data();
         }

      if (Print_pop)
         {
         printf(" Current population (gen %d):\n", Gen.index);
         print_pop(Pop, 0, Pop_size-1);
         }  /* if */

#ifdef STEP
      printf(" Generation %d \n", Gen.index);
      print_population(Pop, 0, Pop_size-1);
      printf(" --- Press any key to continue ---\n");  fgetc(stdin);
#endif

     /* reproduce and put offspring in Kids array */
      if (reproduce() == ERROR)  return ERROR;

/*
printf(" Current population:\n");
print_population(Pop, 0, Pop_size-1);
printf(" Selected parents population:\n");
print_population(Parents, 0, Pop_size-1);
printf(" Kid population:\n");
print_pop(Kids, 0, Pop_size-1);
*/

     /* end of gen operations and outputs */
      if (gen_end() == ERROR)  return ERROR;
      }  /* while */

#ifdef DEBUG
   printf(" ---end ga_loop()---\n");
#endif
   return OK;
   }  /* ga_loop */

/********** ga_continue **********/
/* parameters:
   called by:	ga_loop(), ga.c
   actions:	Checks stopping condition.
		Returns 1 to continue run, 0 to stop run.
*/
int ga_continue()
   {
#ifdef DEBUG
   printf(" ---in ga_continue()---\n");
#endif

   if (Gen.index >= Max_num_gen)  return 0;
   else  
      {
      if (Loop_until_find_opt)
         if (found_solution())  return 0;
         else  return 1;
      else  return 1;
      }

#ifdef DEBUG
   printf(" ---end ga_continue()---\n");
#endif
   return OK;
   }  /* ga_continue */

/********** ga_end **********/
/* parameters:
   called by:	ga_loop(), ga.c
   actions:	Things that need to be done before the GA ends.
		Like deallocating space allocated in ga_init().
*/
int ga_end()
   {
   int i;

#ifdef DEBUG
   printf(" ---in ga_end()---\n");
#endif

  /* finish up any run stats */
   if (run_end() == ERROR) 
      {
      printf(" Error(ga): run_end() ends on error.\n");
      return ERROR;
      }

  /* finish up function stuff and deallocate */
   end_function();

  /* from ga_init() */
   free(Output_path);
   free(Run_num_file);
   free(Init_pop_file);
   // free(Xover_type);
   // free(Mut_type);
   free(Function_name);
   close_pop();

  /* from read_default_opfiles(), params.c */
   for (i=Max_num_output_files-1; i>=0; i--)
      {
      free(Output_file[i].extension);
      if (Output_file[i].on)  free(Output_file[i].filename);
      }  /* for i */
   free(Output_file);

   printf(" End of run %d\n", Run_num);
#ifdef DEBUG
   printf(" ---end ga_end()---\n");
#endif
   return OK;
   }  /* ga_end */
