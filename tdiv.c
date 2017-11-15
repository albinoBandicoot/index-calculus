#include "tdiv.h"

// clobbers gk
relation *tdiv (relation *res, mpz_t gk, uint32_t k, ic_opts *opt) {
	clear_relation (res);
	res->k = k;
	for (int i=0; i < opt->fblen; i++) {
		while (mpz_divisible_ui_p (gk, opt->fb[i])) {
			mpz_divexact_ui (gk, gk, opt->fb[i]);
			res->exp[i] ++;
		}
		if (mpz_cmp_ui (gk, 1) == 0) return res;
	}
	return NULL;
}

