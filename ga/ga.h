/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* ga.h
   06.18.98 AW  Created.
*/

/* prototypes */
int ga_start(char *params_file, char *opfiles_file, char *fxn_file, double xover_rate, double mut_rate, double uniform_xover_rate, char *xover_op, char *mut_op);
int ga_init();
int ga_loop();
int ga_continue();
