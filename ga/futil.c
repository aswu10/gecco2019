/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* Utility functions 
   07/29/94 AW  Created.
   07/30/94 AW  print_goto prints a binary linked tree depth first.
		print_ordered prints a 1-D linked list.
*/

#include <stdio.h>
#include "types.h"
#include "futil.h"

void print_goto(struct state_type *tree)
   {
   int i;

   if (tree == NULL) return;

   if (tree != tree->g[0])
      print_goto(tree->g[0]);

   if (tree->g[0] != NULL)
      {
      printf(" state %d -> [0] -> state %d ", tree->number, tree->g[0]->number);
      printf(" Outputs: ");
      for (i=0; i<tree->g[0]->num_outputs; i++)
         {
         printf("%d ", tree->g[0]->output[i]);
         }  /* for */
      printf("\n");
      }  /* if */
   if (tree->g[1] != NULL)
      {
      printf(" state %d -> [1] -> state %d ", tree->number, tree->g[1]->number);
      printf(" Outputs: ");
      for (i=0; i<tree->g[1]->num_outputs; i++)
         {
         printf("%d ", tree->g[1]->output[i]);
         }  /* for */
      printf("\n");
      }  /* if */

   if (tree != tree->g[1])
      print_goto(tree->g[1]);
   }  /* print_goto */

void print_ordered(struct state_type *head)
   {
   struct state_type *current;
   int i;

   current = head;
   while (current != NULL)
      {
      printf(" State %d:\n", current->number);

      if (current->f == NULL)  printf("    Fail = NULL");
      else  printf("    Fail = %d", current->f->number);

      printf("    Output: ");
      for (i=0; i<current->num_outputs; i++)
         printf("%d ", current->output[i]);
      printf("\n");

      printf("    Length: ");
      for (i=0; i<current->num_outputs; i++)
         printf("%d ", current->length[i]);
      printf("\n");

      if (current->g[0] != NULL)
         {
         printf("    state %d -> [0] -> state %d\n", current->number,
                current->g[0]->number);
         }  /* if */
      if (current->g[1] != NULL)
         {
         printf("    state %d -> [1] -> state %d\n", current->number,
                current->g[1]->number);
         }  /* if */

      current = current->next;
      }  /* while */
   }  /* print_ordered */

void print_queue(struct state_type *head)
   {
   struct state_type *current;

#ifdef DEBUG
   printf(" ---in print_queue---\n");
#endif
   current = head;
   while (current != NULL)
      {
      printf(" state %d, fail = %d\n", current->number,
             current->f->number);
      current = current->fnext;
      }  /* while */
#ifdef DEBUG
   printf(" ---end print_queue---\n");
#endif
   }  /* print_queue */
