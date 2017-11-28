#ifndef TDIV_H
#define TDIV_H

#include <gmp.h>
#include "relation.h"

#define USE_EARLY_ABORT 
#define EARLY_ABORT_FRAC 0.08
#define EARLY_ABORT_THRESH 0.18

int get_ea_idx (ic_opts *opt);
int get_ea_thresh (ic_opts *opt);
relation *tdiv (relation *res, mpz_t gk, uint32_t k, ic_opts *opt) ;

#endif
