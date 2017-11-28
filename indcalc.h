#ifndef INDCALC_H
#define INDCALC_H

#include <gmp.h>
#include "core.h"
#include "tdiv.h"

typedef struct {
	ic_opts *opt;
	matrix *m;
	int thread_index;
} thread_data_t;

matrix *find_relations (ic_opts *opt, int target);
bool verify_solution (mpz_t x, ic_opts *opt) ;

#endif
