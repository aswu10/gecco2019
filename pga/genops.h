/**** Copyright 1998, Annie S. Wu, All Rights Reserved ****/

/* genop.h
   07.24.98.AW	Created.
*/

/* prototypes */
int crossover_pop();
void mutate_pop();
void onept_crossover(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
			INDIVIDUAL *kid1, INDIVIDUAL *kid2);
void twopt_crossover(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
			INDIVIDUAL *kid1, INDIVIDUAL *kid2);
void uniform_crossover(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
			INDIVIDUAL *kid1, INDIVIDUAL *kid2);
void switch_crossover(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
			INDIVIDUAL *kid1, INDIVIDUAL *kid2);
void mutate(INDIVIDUAL *indv);
int homologous_crossover(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
			INDIVIDUAL *kid1, INDIVIDUAL *kid2);
int get_parent2_xpt(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
			int parent1_xpt);
int valid_crossover(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
			int parent1_xpt, int parent2_xpt);
void mutate_pop();
void mutate_alpha(INDIVIDUAL *indv);
void mutate_multichar(INDIVIDUAL *indv);
void mutate_int(INDIVIDUAL *indv);
int poisson(double lambda);
