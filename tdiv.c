#include "tdiv.h"
int get_ea_idx (ic_opts *opt) {
	return (int) (opt->fblen * EARLY_ABORT_FRAC);
}

int get_ea_thresh (ic_opts *opt) {
	return (int) (opt->nbits * EARLY_ABORT_THRESH);
}

// clobbers gk
relation *tdiv (relation *res, mpz_t gk, uint32_t k, ic_opts *opt) {
	clear_relation (res);
	res->k = k;
	int log_approx = 0;
	int early_abort_idx = get_ea_idx(opt);
	for (int i=0; i < opt->fblen; i++) {
		while (mpz_divisible_ui_p (gk, opt->fb[i])) {
			mpz_divexact_ui (gk, gk, opt->fb[i]);
			res->exp[i] ++;
#ifdef USE_EARLY_ABORT
			log_approx += opt->fblogs[i];
#endif
		}
#ifdef USE_EARLY_ABORT
		if (i == early_abort_idx) {
			if (log_approx < get_ea_thresh(opt)) return NULL;
		}
#endif
		if (mpz_cmp_ui (gk, 1) == 0) return res;
	}
	return NULL;
}

