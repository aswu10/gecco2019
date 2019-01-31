/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* main.c
   06.18.98.AW	Created.
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "types.h"
#include "ga.h"
/*
#include "main.h"
#include "init.h"
#include "reproduce.h"
#include "crossover.h"
#include "util.h"
#include "output.h"
#include "foutput.h"
#include "feval.h"
#include "eval.h"
#include "trace.h"
*/

#define PRINT_RANDOM    1	/* print random seed for each run */ 
#define GEN_POPULATION 	1	/* prints out pop. each generation */
#undef GEN_POPULATION

/* To shorten the output that is printed each generation,
   undefining PRINT_BEST_EACH_GEN will prevent printing
   the actual individual.  All the other data for each
   generation is still printed though.
*/
#define PRINT_BEST_EACH_GEN 1

main(int argc, char **argv)
   {
   int error;
   char time_string[100];
   char *s;
   time_t *tstart, *tend;

   if (argc < 4)
      {
      printf(" Usage: %s <parameter file> <opfiles file> <function_file>\n",
		argv[0]);
      return ERROR;
      }  /* if */

   printf(" parameter file: %s\n", argv[1]);
   printf(" opfiles file: %s\n", argv[2]);
   printf(" function file: %s\n", argv[3]);
 
// 19.01.30.AW Time collection moved into ga_start() in ga.c so that
// we can print out the time in an output file.
//  /* get start time */
//   tstart = (time_t *)malloc(sizeof(time_t));
//   time(tstart);
//   strftime(time_string, 99, " %m/%d/%y  %H:%M:%S ", localtime(tstart)); 
//   printf(" Start time: %s\n", time_string);


  /* start the ga -- this is a separate function call to make
     it easier to integrate into other programs, should that
     ever become necessary.  main() is simply a shell. */
   error = ga_start(argv[1], argv[2], argv[3], strtod(argv[4], NULL), strtod(argv[5], NULL), strtod(argv[6], NULL), argv[7], argv[8]);
   if (error == -1)  printf(" Error(main): ga_start ends on error.\n");
   else              printf(" --- ga_start ends ok. ---\n");

// 19.01.30.AW Time collection moved into ga_start() in ga.c so that
// we can print out the time in an output file.
//  /* get end time */
//   tend = (time_t *)malloc(sizeof(time_t));
//   time(tend);
//   strftime(time_string, 99, " %m/%d/%y  %H:%M:%S ", localtime(tend));
//   printf(" End time: %s\n", time_string);

  /* Added on 1.19.96 AW */
//   free(tstart);
//   free(tend);

   return OK;
   }  /* main */
