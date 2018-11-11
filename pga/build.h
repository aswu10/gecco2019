/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

#ifndef _BUILD_H
#define _BUILD_H

struct state_type *build_goto(int num_keywords, GA_WORD *list);

void enter_word(int length, char *keyword, int bb_name, int *state_ctr,
                struct state_type *state, int num_keywords,
                struct state_type **ordered);

struct state_type *make_new_state(int *state_ctr, int num_keywords);

void build_fail(struct state_type *state);

#endif /* _BUILD_H */
