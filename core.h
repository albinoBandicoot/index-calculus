#ifndef COMMON_H
#define COMMON_H

#include <gmp.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NEW(X) (X*) calloc(sizeof(X), 1)
typedef struct {
	mpz_t g;	// base of logarithm
	mpz_t q;	// prime modulus
	mpz_t h;	// target
} instance;

typedef struct {
	instance in;
	uint32_t *fb;
	uint32_t *fblogs;
	uint32_t nbits;
	uint32_t fblen;
	uint32_t nthreads;
} ic_opts;

ic_opts *init (mpz_t g, mpz_t q, mpz_t h, uint32_t fbb, uint32_t nthreads);
void generate_fb (ic_opts *opt, uint32_t fbb);



#endif
