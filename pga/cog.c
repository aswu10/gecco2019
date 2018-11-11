/**** Copyright 2001, Annie S. Wu, All Rights Reserved ****/

/* cog.c
   010412AW	Created.  Routines that calculate center of gravity.
*/

#include <stdio.h>
#include <math.h>
#include "types.h"
#include "extern.h"
#include "cog.h"

#define DEBUG_COG 1
#undef DEBUG_COG

/********** init_cog **********/
/* parameters:
   created:	010412AW
   called by:	ga_init(), ga.c
   actions:	init stuff for recording CoG - make cog directory
		and open file pointers.
*/
int init_cog()
   {
   char tempstring[LONG_LINE_LEN];
   int error;
 
#ifdef DEBUG
   printf(" ---in init_cog---\n");
#endif
 
  /* make the cog directory */
   sprintf(tempstring, "mkdir %s/run.%d/cog", Output_path, Run_num);

#ifdef CMD
   printf(" CMD: %s\n", tempstring);
#endif
   error = system(tempstring);
   if (error != 0)
      printf(" Error(init_cog): with system command: %s\n", tempstring);
 
#ifdef DEBUG
   printf(" ---end init_cog---\n");
#endif
   return OK;
   }  /* init_cog */

/********** cog_indv **********/
/* parameters:  indv
   created:     010412AW
   called by:
   actions:     Calculate center of gravity for each xter on indv.
                If xter not on indv, set to -1.0.
*/
void cog_indv(INDIVIDUAL *indv)
   {
   int i, x;

#ifdef DEBUG_COG
printf(" %s\n", indv->genome);
#endif
 
   for (i=Alphabet_size-1; i>=0; i--)
      {
      indv->cog_count[i] = 0;
      indv->cog_avg[i] = 0.0;
      indv->cog_low[i] = -0.1;
      indv->cog_high[i] = -0.1;
      }  /* for i */

  /* count number of each xter, sum locations, mark end locations */
   for (i=indv->length-1; i>=0; i--)
      {
      x = which_xter(indv->genome, i);
#ifdef DEBUG_COG
printf(" which_xter returned %d\n", x);
#endif
      if (x >= 0)
         {
         indv->cog_count[x]++;
         indv->cog_avg[x] = indv->cog_avg[x] + i;
         indv->cog_low[x] = i;
         if (i > indv->cog_high[x])  indv->cog_high[x] = i;
         }  /* if */
      }  /* for i */
 
  /* normalize with respect to length of individual */
   for (i=Alphabet_size-1; i>=0; i--)
      {
      if (indv->cog_count[i] == 0)
         {
         indv->cog_avg[i] = -0.1;
         }  /* if */
      else
         {
         indv->cog_avg[i] = indv->cog_avg[i]/indv->cog_count[i];
        /* normalize to length of indv */
         indv->cog_avg[i] = indv->cog_avg[i]/indv->length;
         indv->cog_low[i] = indv->cog_low[i]/indv->length;
         indv->cog_high[i] = indv->cog_high[i]/indv->length;
         }  /* else */
#ifdef DEBUG_COG
printf(" %lf", indv->cog_avg[i]);
#endif
      }  /* for i */
#ifdef DEBUG_COG
printf("\n");
#endif
   }  /* cog_indv */

/********** which_xter **********/
/* parameters:
   created:	010412AW
   called by:	cog_indv()
   actions:	returns integer of xter in Xter.
*/
int which_xter(char *genome, int which)
   {
   int i;

   for (i=Alphabet_size-1; i>=0; i--)
      {
#ifdef DEUBG_COG
printf(" compare %s:%d and %s:%d\n", genome, which, Xters, i);
#endif
      if (genome[which] == Xters[i])
         return i;
      }  /* for i */
   return(-1);
   }  /* which_xter */

/********** cog_gen_end **********/
/* parameters:
   created:	010413AW
   called by:	gen_end(), stats.c ???
   actions:	outputs stuff at end of generation
*/
int cog_gen_end()
   {
   FILE *fp;
   char tempstring[INPUT_LINE_LEN];
   int i;

#ifdef DEBUG
   printf(" ---in cog_gen_end---\n");
#endif

  /* print to file trace/cog/cog.genx */
   sprintf(tempstring, "%s/run.%d/cog/cog.gen%d",
                Output_path, Run_num, Gen.index);

   fp = fopen(tempstring, "w");
   fprint_gen_cog(fp);
/*
   fprint_gen_cog_allsame(fp);
*/
   fclose(fp);

#ifdef DEBUG
   printf(" ---end cog_gen_end---\n");
#endif
   return OK;
   }  /* cog_gen_end */

/********** fprint_gen_cog **********/
/* parameters:	fp	file to print to
   created:	010413AW
   called by:	cog_gen_end()
   actions:	prints cog data for current generation.
		must print for each indv in generation.

   NOTE:	Can distinguish the cog's of the various xters.
*/
void fprint_gen_cog(FILE *fp)
   {
   int i, j;

   for (i=0; i<Pop_size; i++)
      {
      fprintf(fp, " %4d ", i);
      for (j=0; j<Alphabet_size; j++)
         {
         fprintf(fp, " %lf %lf %lf ", Pop[i]->cog_avg[j],
                 Pop[i]->cog_low[j], Pop[i]->cog_high[j]);
         }  /* for j */
      fprintf(fp, "\n");
      }  /* for i */

   }  /* fprint_gen_cog */

/********** fprint_gen_cog_allsame **********/
/* parameters:  fp      file to print to
   created:     010413AW
   called by:   cog_gen_end()
   actions:     prints cog data for current generation.
                must print for each indv in generation.

   NOTE:	CANNOT distinguish the cog's of the various xters.
*/
void fprint_gen_cog_allsame(FILE *fp)
   {
   int i, j;

   for (i=0; i<Pop_size; i++)
      {
      for (j=0; j<Alphabet_size; j++)
         {
         fprintf(fp, " %4d  %lf %lf %lf ", i, Pop[i]->cog_avg[j],
                 Pop[i]->cog_low[j], Pop[i]->cog_high[j]);
         }  /* for j */
      }  /* for i */

   }  /* fprint_gen_cog */

