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
void position_based_crossover(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
			INDIVIDUAL *kid1, INDIVIDUAL *kid2);
void edge_recombination_crossover(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
			INDIVIDUAL *kid1, INDIVIDUAL *kid2);
void partially_mapped_crossover(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
			INDIVIDUAL *kid1, INDIVIDUAL *kid2);
            
void get_edges_from_parent(EDGE_MAP *maps, INDIVIDUAL *parent);
void edge_recombination_to_child(INDIVIDUAL *child, EDGE_MAP *maps, INDIVIDUAL *parent1, INDIVIDUAL *parent2);
void remove_from_maps(int value, EDGE_MAP *maps, int length, int *unvisited, int num_unvisited);


void random_mutate(INDIVIDUAL *indv);
void displacement_mutate(INDIVIDUAL *indv, int start, int length, int invert, int per_gene);
int homologous_crossover(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
			INDIVIDUAL *kid1, INDIVIDUAL *kid2);
int get_parent2_xpt(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
			int parent1_xpt);
int valid_crossover(INDIVIDUAL *parent1, INDIVIDUAL *parent2,
			int parent1_xpt, int parent2_xpt);
