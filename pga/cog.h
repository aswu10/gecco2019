/**** Copyright 2001, Annie S. Wu, All Rights Reserved ****/

/* cog.h */

int init_cog();
void cog_indv(INDIVIDUAL *indv);
int which_xter(char *genome, int which);
int cog_gen_end();
void fprint_gen_cog(FILE *fp);
void fprint_gen_cog_allsame(FILE *fp);
