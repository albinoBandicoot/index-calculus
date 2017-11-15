#include "indcalc.h"
#include "matsolve.h"

matrix *find_relations (ic_opts *opt) {
	matrix *m = new_matrix(opt->fblen * 3);
	mpz_t gk, gksave;
	mpz_inits (gk, gksave, NULL);
	mpz_set (gk, opt->in.g);
	mpz_set (gksave, gk);
	uint32_t k = 1;

	relation *rel = new_relation(opt);

	while (1) {
//		gmp_printf ("k = %d; g^k (mod q) = %Zd\n", k, gk);
		relation *r = tdiv(rel, gk, k, opt);
		if (r != NULL) {
			if (add_row(m, r)) {
				mpz_clears (gk, gksave, NULL);
				return m;
			}
			printf ("k = %d; have [%d/%d]\n", k, m->next_slot, m->nrows);
			rel = new_relation(opt);
		}
		mpz_mul (gk, gksave, opt->in.g);
		mpz_mod (gk, gk, opt->in.q);
		mpz_set (gksave, gk);
		k++;
	}
}

void deduce_solution (mpz_t res, relation *r, solution *fbsol, uint64_t s) {
	mpz_t tmp;
	mpz_init (tmp);
	mpz_set_ui (res, 0);
	for (int i=0; i < r->size; i++) {
		mpz_mul_ui (tmp, fbsol->d[i], r->exp[i]);
		mpz_add (res, res, tmp);
	}
	mpz_sub_ui (res, res, s);
	mpz_clear (tmp);
}

void trial_factor_gsh (mpz_t res, ic_opts *opt, solution *fbsol) {
	mpz_t gsh, gshsave;
	mpz_inits (gsh, gshsave, NULL);
	mpz_set (gsh, opt->in.h);
	mpz_set (gshsave, gsh);
	uint64_t s = 0;

	relation *rel = new_relation(opt);
	while (1) {
		printf ("Trying to factor g^%llu * h\n", s);
		relation *r = tdiv(rel, gsh, s, opt);
		if (r != NULL) {
			// g^s * h factored over the factor base.
			printf ("SMOOTH!\n");
			print_relation (r);
			deduce_solution (res, r, fbsol, s);
			mpz_mod (res, res, opt->in.q);
			return;
		}
		mpz_mul (gsh, gshsave, opt->in.g);
		mpz_mod (gsh, gsh, opt->in.q);
		mpz_set (gshsave, gsh);
		s++;
	}
}

bool verify_solution (mpz_t x, ic_opts *opt) {
	mpz_t tmp;
	mpz_init (tmp);
	mpz_powm (tmp, opt->in.g, x, opt->in.q);
	gmp_printf ("g^x (mod q) = %Zd, should be %Zd\n", tmp, opt->in.h);
	bool res = mpz_cmp (tmp, opt->in.h) == 0;
	mpz_clear (tmp);
	return res;
}

int main (int argc, const char *argv[]) {
	mpz_t g,q,h,x;
	mpz_inits (g,q,h,x,NULL);

	mpz_set_str(g, argv[1], 10);
	mpz_set_str(q, argv[2], 10);
	mpz_set_str(h, argv[3], 10);

	int fbb = atoi (argv[4]);

	ic_opts *opt = init(g,q,h,fbb);
	printf ("Sieving done; found %d primes below %d\n", opt->fblen, fbb);

	matrix *m = find_relations (opt);
	print_matrix(m);

	//solution_ui *soln = solve_all_ui (m, opt, NULL);
	solution_ui *soln = solve_ui (m, mpz_get_ui(q) - 1, 0);

	print_solution_ui (soln);
	solution *sol = convert_up (soln);
	printf ("Converted up\n");

	trial_factor_gsh (x, opt, sol);
	gmp_printf ("Found solution: %Zd^%Zd = %Zd (mod %Zd)\n", g, x, h, q);
	if (verify_solution (x, opt)) {
		printf ("Verified.\n");
	} else {
		printf ("Incorrect.\n");
	}
}
