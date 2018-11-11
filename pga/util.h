/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* 09.19.93 AW  Created.
*/

/* prototypes */
/***** 031208AW Comment out and use same fxn in math.h
double fmin(double a, double b);
double fmax(double a, double b);
*/
char *get_date_time();
int int_in_list(int val, int *alist, int len);
int which_parent(INDIVIDUAL *indv, int bit);
int bin_to_int_pos(char *bstring, int left_bit, int right_bit);
int ipow(int x, int y);
int bin_to_int(char *bstring, int left_bit, int right_bit);
