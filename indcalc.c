#include <pthread.h>
#include <time.h>
#include "indcalc.h"
#include "matsolve.h"
#include "tdiv.h"

void *run_worker_thread (void *args) {
	thread_data_t *tdata = (thread_data_t *) args;
	ic_opts *opt = tdata->opt;
	matrix *m = tdata->m;
	int thread_idx = tdata->thread_index;

	mpz_t gk, gnt, gksave;
	mpz_inits (gk, gnt, gksave, NULL);
	printf ("Starting k at %d, will increment by %d\n", thread_idx + 1, opt->nthreads);
	mpz_pow_ui (gk, opt->in.g, thread_idx + 1);
	mpz_pow_ui (gnt, opt->in.g, opt->nthreads);
	mpz_mod (gk, gk, opt->in.q);
	mpz_mod (gnt, gnt, opt->in.q);
	mpz_set (gksave, gk);
	uint32_t k = thread_idx + 1;

	relation *rel = new_relation(opt);

	while (1) {
		if (is_matrix_full(m)) {
			mpz_clears (gk, gnt, NULL);
			printf("\n");
			return NULL;
		}
//		gmp_printf ("k = %d; g^k (mod q) = %Zd\n", k, gk);
		relation *r = tdiv(rel, gk, k, opt);	// trashes gk
		if (r != NULL) {
			if (add_row(m, r)) {
				mpz_clears (gk, gnt, NULL);
				printf("\n");
				return NULL;
			}
			printf ("k = %d; have [%d/%d]\r", k, m->next_slot, m->nrows);
			fflush(stdout);
			rel = new_relation(opt);
		}
		mpz_mul (gksave, gksave, gnt);
		mpz_mod (gksave, gksave, opt->in.q);
		mpz_set (gk, gksave);
		k += opt->nthreads;
	}
}

matrix *find_relations (ic_opts *opt, int target) {
	if (target == 0) target = opt->fblen * 3;
	matrix *m = new_matrix(target);
	pthread_t workers[opt->nthreads];
	thread_data_t tdata[opt->nthreads];
	for (int i=0; i < opt->nthreads; i++) {
		tdata[i].opt = opt;
		tdata[i].m = m;
		tdata[i].thread_index = i;
		pthread_create (&workers[i], NULL, run_worker_thread, (void *) &tdata[i]);
	}
	for (int i=0; i < opt->nthreads; i++) {
		pthread_join (workers[i], NULL);
	}
	return m;
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

void verify_fbsol (solution *fbsol, ic_opts *opt) {
	mpz_t tmp;
	mpz_init (tmp);
	for (int i=0; i < opt->fblen; i++) {
		mpz_powm (tmp, opt->in.g, fbsol->d[i], opt->in.q);
		if (mpz_cmp_ui (tmp, opt->fb[i]) != 0) {
			gmp_printf ("Log of fb[%d] is wrong. Our solution: %Zd; gets value %Zd\n", i, fbsol->d[i], tmp);
		}
	}
}

void trial_factor_gsh (mpz_t res, ic_opts *opt, solution *fbsol) {
	mpz_t gsh, gshsave, qm1;
	mpz_inits (gsh, gshsave, qm1, NULL);
	mpz_set (gsh, opt->in.h);
	mpz_set (gshsave, gsh);
	mpz_sub_ui (qm1, opt->in.q, 1);
	uint64_t s = 0;

	relation *rel = new_relation(opt);
	while (1) {
//		printf ("Trying to factor g^%llu * h\n", s);
		relation *r = tdiv(rel, gsh, s, opt);
		if (r != NULL) {
			// g^s * h factored over the factor base.
			printf ("g^%llu * h was smooth!\n", s);
			print_relation (r);
			deduce_solution (res, r, fbsol, s);
			mpz_mod (res, res, qm1);
	//		if (verify_solution (res, opt)) return;
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

	if (argc < 6) {
		printf ("Usage:  ind g q h fbb nthreads\nWill find x such that g^x = h (mod q)\n");
		exit(0);
	}

	mpz_set_str(g, argv[1], 10);
	mpz_set_str(q, argv[2], 10);
	mpz_set_str(h, argv[3], 10);

	int fbb = atoi (argv[4]);

	int nt = atoi (argv[5]);

	int target = 0;
	if (argc == 7) {
		target = atoi(argv[6]);
	}

	long start = clock();
	long init_time;
	long stage1_time;
	long stage2_time;
	long stage3_time;

	long time = clock();
	ic_opts *opt = init(g,q,h,fbb,nt);
	printf ("Sieving done; found %d primes below %d\n", opt->fblen, fbb);
	printf ("Will use early-abort index of %d (p=%d); candidate must have accumulated %d bits by then\n", get_ea_idx(opt), opt->fb[get_ea_idx(opt)], get_ea_thresh(opt));
	init_time = clock() - time;
	time = clock();

	matrix *m = find_relations (opt, target);

	stage1_time = clock() - time;
	time = clock();

	//print_matrix(m);

	mpz_t n;
	mpz_init_set (n, q);
	mpz_sub_ui (n,n,1);
	solution *sol = solve (m, n, generic_inversion_test);
	//solution *sol = solve_all (m, opt, NULL);
	//solution_ui *soln = solve_all_ui (m, opt, NULL);
	printf ("Did %d modular inverses\n", INV_CT);

//	print_solution (sol);
	//print_solution_ui (soln);
	//solution *sol = convert_up (soln);
	
	stage2_time = clock() - time;
	time = clock();

	verify_fbsol (sol, opt);
	trial_factor_gsh (x, opt, sol);
	gmp_printf ("Found solution: %Zd^%Zd = %Zd (mod %Zd)\n", g, x, h, q);
	if (verify_solution (x, opt)) {
		printf ("Verified.\n");
	} else {
		printf ("Incorrect.\n");
	}
	stage3_time = clock() - time;

	printf ("Time summary: \n");
	printf ("\tInitialization:          %ldms\n", init_time/1000);
	printf ("\tStage 1: (rel finding):  %ldms\n", stage1_time/1000);
	printf ("\tStage 2: (lin. alg.):    %ldms\n", stage2_time/1000);
	printf ("\tStage 3: (deduce log):   %ldms\n", stage3_time/1000);
	printf ("\tTOTAL:                   %ldms\n", (clock()-start)/1000);
}
