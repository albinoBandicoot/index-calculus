#ifndef MATSOLVE_H
#define MATSOLVE_H

#include <stdint.h>
#include <stdio.h>
#include "relation.h"
#include "core.h"

typedef uint32_t matelem_ui;

typedef struct {
	matelem_ui *d;
	uint64_t p;
	uint64_t pk;
	uint32_t e;
	uint32_t size;
} solution_ui;

typedef struct {
	mpz_t *d;
	mpz_t p;	// factor of q-1 that this solution corresponds to
	mpz_t pk;
	uint32_t e;	// exponent on p. 
	uint32_t size;
} solution;

solution_ui *new_solution_ui (uint64_t p, uint32_t e, int len);
solution *new_solution (mpz_t p, uint32_t e, int len);
solution *convert_up (solution_ui *);
void free_solution_ui (solution_ui *sol) ;
void free_solution (solution *sol) ;
void print_solution_ui (solution_ui *sol);

solution_ui *solve_ui (matrix *m, uint64_t p, uint32_t e);
solution *solve (matrix *m, mpz_t p, uint32_t e);

void solution_ui_merge_ui ();
void solution_merge (solution *dest, solution *src);

solution *solve_all (matrix *m, ic_opts *opt, FILE *factors) ;
solution_ui *solve_all_ui (matrix *m, ic_opts *opt, FILE *factors) ;
#endif
