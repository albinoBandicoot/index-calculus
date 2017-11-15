#include "core.h"

ic_opts *init (mpz_t g, mpz_t q, mpz_t h, uint32_t fbb) {
	ic_opts *opt = NEW(ic_opts);
	mpz_inits (opt->in.g, opt->in.q, opt->in.h, NULL);
	mpz_set (opt->in.g, g);
	mpz_set (opt->in.q, q);
	mpz_set (opt->in.h, h);
	generate_fb (opt, fbb);
	return opt;
}

/* The familiar sieve of Eratosthenes */
void era_sieve (char *vals, uint32_t fbb){
	// assumes vals has been allocated and cleared with sufficient space for fb_bound chars.
	// 0 is prime, 1 is composite.
	int search_max = (int) (sqrt(fbb) + 1);
	for (int skip=2; skip < search_max; skip++){
		if (vals[skip-2] == 1){
			continue;
		}
		for (int pos=2*skip; pos < fbb; pos += skip){
			vals[pos-2] = 1;
		}
	}
}

/* Extract the values of the primes from the sieve. */
void extract (ic_opts *opt, char *vals, uint32_t fbb){
	int count = 0;
	for (int i=0; i < fbb-2; i++){
		if (vals[i] == 0){
			count++;
		}
	}
	// now we know how much space to allocate for our list
	opt->fblen = count;
	opt->fb = (uint32_t *)(malloc(count * sizeof(uint32_t)));
	int w = 0;
	for (int i=0; i < fbb-2; i++){
		if (vals[i] == 0) {
			opt->fb[w] = i+2;
			w++;
		}
	}
}

void generate_fb (ic_opts *opt, uint32_t fbb) {
	char *erasieve_vals = (char *) calloc (fbb, 1);
	era_sieve (erasieve_vals, fbb);
	extract (opt, erasieve_vals, fbb);
	free (erasieve_vals);
}


