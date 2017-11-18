#include <gmp.h>
#include <stdio.h>

int main (int argc, const char *argv[]) {
	int nbits = atoi(argv[1]);
	mpz_t p, q;
	mpz_inits (p, q, NULL);
	gmp_randstate_t rand;
	gmp_randinit_default (rand);

	mpz_urandomb (p, rand, nbits-1);
	do {
		mpz_nextprime(p,p);
		mpz_mul_ui (q, p, 2);
		mpz_add_ui (q, q, 1);
	} while (!mpz_probab_prime_p(q, 30));
	gmp_printf ("%Zd\n", q);
	mpz_clears(p,q,NULL);
}
