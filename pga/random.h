extern long seed_random();
extern double knuth_random();

/* Use a seed value less than 0 for the seeder to pick its own seed. */

#define NO_SEED                -1L             

/* 
 * These 3 return 1 if it is time to mutate or crossover
 * or to choose the fitter individual in tournament selection.
 * They have to be used in the context of a CONTEXT * called "context".
 */
#define DO_MUTATE              (knuth_random() < context->mutation_prob)
#define DO_CROSSOVER           (knuth_random() < context->crossover_prob)
#define DO_TOURNAMENT          (knuth_random() < context->tournament_prob)
#define DO_UNIFORM             (knuth_random() < context->uniform_prob)

/* This returns [0, n - 1] */
/*
#define uniform(n)             ((int) (knuth_random() * (double)(n)))
*/

/* 
 * Returns true if a random double is less than the given probability. 
 */
#define biased_coin_flip(prob) (knuth_random() < (prob))

/******* the above stuff is Terry's *******/
/******* prototypes *******/

int uniform(int n);
double funiform(double x);
double xfuniform(double x);
