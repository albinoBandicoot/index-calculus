#ifndef TDIV_H
#define TDIV_H

#include <gmp.h>
#include "relation.h"

relation *tdiv (relation *res, mpz_t gk, uint32_t k, ic_opts *opt) ;

#endif
