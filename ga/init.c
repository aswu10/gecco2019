/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* init.c
   07.07.98.AW	Created.
		Stuff for initializing a run.
		For now, mainly run number and random seed stuff.

   Routines:	init_run()		07.07.98.AW
		get_run_num()		07.07.98.AW
		get_random_seed()	07.07.98.AW
*/

#include <stdio.h>
#include "types.h"
#include "extern.h"
#include "init.h"
#include "random.h"

/********** init_run **********/
/* parameters:
   called by:	ga_init(), ga.c
   actions:	general routine to initialize stuff for a run.
		Gets new run number.
		Creates new directory in output directory.
		Gets random seed.
		Calculate global variables here.
*/
int init_run()
   {
#ifdef DEBUG
   printf(" ---in init_run()---\n");
#endif

  /* set up random number generator */
   if (get_random_seed() == ERROR)  return ERROR;

  /* get new run number -- index to identify run */
   if (get_run_num() == ERROR)  return ERROR;
 
  /* create output directory */
   if (make_directory() == ERROR)  return ERROR;

  /* set global variables */
   Found_optimum = 0;
  /* calculate global variables */
   Num_breeding = Pop_size * Pct_breeding;
   Num_bred = Pop_size * Pct_bred;
   if (!strcmp(Base, "binary"))  Alphabet_size = 2;
      else if (!strcmp(Base, "alphabet"))  Alphabet_size = 26;
      else  Alphabet_size = 0;
   if (!strcmp(Base, "binary"))  Ascii_offset = 48;
      else if (!strcmp(Base, "alphabet"))  Ascii_offset = 97;
      else  Ascii_offset = 0;
  /* check and make sure even */
   if (Num_breeding % 2 != 0)  Num_breeding--;
   if (Num_bred % 2 != 0)  Num_bred--;
  /* 10.09.98.AW check Hx_window and Min_gen_len if Variable_gen_len */
  /* May want to get rid of this restriction and just make it so that
     xover doesn't occur if one individual is too short?  Or maybe not. */
   if (Variable_gen_len)
      if (Min_gen_len < Hx_window)
         {
         printf(" Error(init_run): Min_gen_len(%d) < Hx_window(%d)",
			Min_gen_len, Hx_window);
         printf(" not allowed\n");
         return ERROR;
         }

#ifdef DEBUG
   printf(" ---end init_run()---\n");
#endif
   return OK;
   }  /* init_run */

/********** get_run_num **********/
/* parameters:
   called by:	init_run(), init.c
   actions:	Get run number for this run.
		Reads index of previous run from Run_num_file,
		increments to get current run number, and writes
		current run number back to Run_num_file.

		The program will always assign a new number to each
		new run, even if the user is rerunning an old run.
		Rerunning an old run simply makes the program read
		in the random seed for the run from a user specified
		file.  The run still gets a new number so that it will
		not conflict with existing files and directories in
		the output directory.
*/
int get_run_num()
   {
   FILE *fp;

  /* read previous run number */
   fp = fopen(Run_num_file, "r");
   if (fp == NULL)
      {
      printf(" Error(get_run_num): cannot open run file: %s\n", Run_num_file);
      return ERROR;
      }
   fscanf(fp, "%d", &Run_num);
   fclose(fp);

  /* assign and write back new number */
   Run_num = Run_num + 1;

   fp = fopen(Run_num_file, "w");
   if (fp == NULL)
      {
      printf(" Error(get_run_num): cannot open run file: %s\n", Run_num_file);
      return ERROR;
      }
   fprintf(fp, "%d\n", Run_num);
   fclose(fp);
   
   return OK;
   }  /* get_run_num */

/********** make_directory **********/
/* parameters:
   called by:   init_run(), init.c
   actions:	make directory for output files from this run
		within the Output_path directory.
		Directory will be called run.<Run_num>
*/
int make_directory()
   {
   FILE *fp;
   char tempstring[LONG_LINE_LEN];
   int error;

   sprintf(tempstring, "mkdir %s/run.%d", Output_path, Run_num);
   error = system(tempstring);
   while (error != 0) // if run_num exists, redo to create a need seed/run_num
      {
      printf("Seed exists, getting new seed");
      get_random_seed();
      get_run_num();
      
      sprintf(tempstring, "mkdir %s/run.%d", Output_path, Run_num);
      error = system(tempstring);
      
      printf("Error(make_directory): with system command: %s\n", tempstring);

      sleep(1);
      }
   
   /* print seed to random file for this run if "random" opfiles
     flag is turned on */
   if (Output_file[get_file_pointer("random")].on)
      {
      sprintf(tempstring, "%s/run.%d/run.%d.random",
             Output_path, Run_num, Run_num);
      fp = fopen(tempstring, "w");
      if (fp == NULL)
         {
         printf(" Error(get_random_seed): cannot open file: %s\n", tempstring);
         return ERROR;
         }
      fprintf(fp, "%ld\n", Seed);
      fclose(fp);
      }  /* if print to random file */
      
   return OK;
   }  /* make_directory */

/********** get_random_seed **********/
/* parameters:
   called by:   init_run(), init.c
   actions:     Get random seed for this run.
		If Rerun = -1, the generate a new random seed for this run.
		If Rerun >= 0, then its values specifies the run that we
		want to rerun.  Read the random seed from that run's random
		seed file, i.e. run.<Rerun>.random.

		Need to have already called get_run_num() before calling
		this program because need run number to write random seed
		to file called run.<Run_num>.random.
*/
int get_random_seed()
   {
   FILE *fp;
   char tempstring[LONG_LINE_LEN];
   long temp;

   if (Rerun < 0)
      {
      Seed = seed_random((long)Rerun);
      }  /* if generate new seed */
   else
      {
      sprintf(tempstring, "%s/run.%d/run.%d.random",
		Output_path, Rerun, Rerun);
      fp = fopen(tempstring, "r");
      if (fp == NULL)
         {
         printf(" Error(get_random_seed): cannot open file: %s\n", tempstring);
         return ERROR;
         }
      fscanf(fp, "%ld", &temp);
      Seed = seed_random(temp);
      fclose(fp);
      }  /* else read old seed */

  // /* print seed to random file for this run if "random" opfiles
     // flag is turned on */
   // if (Output_file[get_file_pointer("random")].on)
      // {
      // sprintf(tempstring, "%s/run.%d/run.%d.random",
             // Output_path, Run_num, Run_num);
      // fp = fopen(tempstring, "w");
      // if (fp == NULL)
         // {
         // printf(" Error(get_random_seed): cannot open file: %s\n", tempstring);
         // return ERROR;
         // }
      // fprintf(fp, "%ld\n", Seed);
      // fclose(fp);
      // }  /* if print to random file */

   return OK;
   }  /* get_random_seed */
