#  define RANDIN  seed_in((long *)NULL)
#  define RANDOUT seed_out((long *)NULL)
#  define UNIF unif_rand()
#  define Salloc(n, t) (t *)S_alloc(n, sizeof(t))
#  define S_EVALUATOR

typedef double singl;

