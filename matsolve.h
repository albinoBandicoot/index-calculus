#ifndef MATSOLVE_H
#define MATSOLVE_H

#include <stdint.h>
#include <stdio.h>
#include "relation.h"
#include "core.h"

extern int INV_CT;

typedef struct {
	mpz_t *d;
	mpz_t p;	// factor of q-1 that this solution corresponds to
	uint32_t size;
} solution;

solution *new_solution (mpz_t p, int len);
void free_solution (solution *sol) ;
void print_solution (solution *sol);

int generic_inversion_test (mpz_t a, mpz_t n) ;
int safe_prime_inversion_test (mpz_t a, mpz_t n) ;

solution *solve (matrix *m, mpz_t n, int (*has_inverse)(mpz_t, mpz_t));

#endif
