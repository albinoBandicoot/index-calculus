#include "matsolve.h"
solution_ui *new_solution_ui (uint64_t p, uint32_t e, int len) {
	solution_ui *sol = NEW(solution_ui);
	sol->d = malloc (8 * len);
	sol->pk = p;
	sol->p = p;
	sol->e = e;
	sol->size = len;
	for (int i=0; i < e; i++) {
		sol->pk *= p;
	}
	return sol;
}
	
solution *new_solution (mpz_t p, uint32_t e, int len) {
	solution *sol = NEW(solution);
	sol->d = malloc (sizeof(mpz_t) * len);
	mpz_init_set (sol->pk, p);
	mpz_init_set (sol->p, p);
	sol->e = e;
	sol->size = len;
	for (int i=1; i < e; i++) {
		mpz_mul (sol->pk, sol->pk, p);
	}
	return sol;
}

solution *convert_up (solution_ui *sol) {
	mpz_t p;
	mpz_init_set_ui (p, sol->p);
	printf ("A\n");
	solution *s = new_solution (p, sol->e, sol->size);
	printf ("B\n");
	for (int i=0; i < sol->size; i++) {
		mpz_init_set_ui (s->d[i], sol->d[i]);
		printf ("%d\n", i);
	}
	mpz_clear(p);
	return s;
}

void free_solution_ui (solution_ui *sol) {
	free (sol->d);
	free (sol);
}

void free_solution (solution *sol) {
	mpz_clears (sol->p, sol->pk, NULL);
	for (int i=0; i < sol->size; i++) {
		mpz_clear (sol->d[i]);
	}
	free (sol->d);
	free (sol);
}

void print_solution_ui (solution_ui *sol) {
	for (int i=0; i < sol->size; i++) {
		printf ("L_%d = %u (mod %llu)\n", i, sol->d[i], sol->p);
	}
}

void mul_row_ui (matelem_ui *row, uint64_t mul, uint64_t modulus, int len) {
	for (int i=0; i < len; i++) {
		row[i] = (matelem_ui) (((uint64_t) row[i] * mul) % modulus);
	}
}

void add_row_multiple_ui (matelem_ui *dest, matelem_ui *src, uint64_t mul, uint64_t modulus, int len) {
	for (int i=0; i < len; i++) {
		dest[i] = (matelem_ui) (((uint64_t) dest[i] + src[i]*mul) % modulus);
	}
}

uint64_t mod_inverse (uint64_t x, uint64_t p) {
	// FIXME: do this properly.
	mpz_t a, b;
	mpz_inits (a,b,NULL);
	mpz_set_ui(a, x);
	mpz_set_ui(b, p);
	int has_inverse = mpz_invert(a, a, b);
	if (has_inverse == 0) {
		mpz_clears(a,b,NULL);
		return 0;
	}
	uint64_t res = mpz_get_ui (a);
	mpz_clears (a,b, NULL);
	return res;
}

solution_ui *solve_ui (matrix *mat, uint64_t p, uint32_t e) {
	printf ("Will solve matrix mod %llu\n", p);
	int nr = mat->nrows;
	int nc = mat->rows[0]->size + 1;
	matelem_ui *m = (matelem_ui *) calloc (nr*nc, 8);
	for (int r = 0; r < nr; r++) {
		for (int c=0; c < nc-1; c++) {
			m[r*nc + c] = mat->rows[r]->exp[c] % p;
		}
		m[r*nc + nc-1] = mat->rows[r]->k % p;
	}
	
	for (int r = 0; r < nr; r++) {
		for (int c=0; c < nc; c++) {
			printf ("%u ", m[r*nc + c]);
		}
		printf ("\n");
	}
		printf ("\n");


	for (int c=0; c < nc-1; c++) {
//		printf ("Column %d\n", c);
		if (m[c*nc + c] == 0) {	// this is bad: a 0 on the diagonal.
//			printf ("  leading zero, must fix.\n");
			// try to find a row below this with a nonzero coefficient in this row, and add it to this row.
			for (int r = c+1; r < nr; r++) {
//				printf ("    trying row %d\n", r);
				if (m[r*nc + c] != 0) {
					uint64_t inv = mod_inverse (m[r*nc + c], p);
					if (inv == 0) continue;
					add_row_multiple_ui (&m[c*nc], &m[r*nc], inv, p, nc);	// may as well make the leading coeff a 1 while we're at it.
					break;
				}
			}
			if (m[c*nc + c] == 0) {	// error
				printf ("  no rows below with nonzero entries!\n");
				printf ("NO INFORMATION ABOUT fb[%d]\n", c);
				continue;	// skip this column
			}
		}
		if (m[c*nc + c] != 1) {	// multiply this row by the multiplicative inverse of the leading entry.
			uint64_t inv = mod_inverse (m[c*nc + c], p);
//			printf ("  correcting leading entry to = 1; multiplying by %llu\n", inv);
			mul_row_ui (&m[c*nc], inv, p, nc);
		}
		// for all other rows, add a multiple to cancel out the term in this column.
		for (int r=0; r < nr; r++) {
			if (r == c) continue;
			if (m[r*nc + c] != 0) {
//				printf ("   adding %llu * row %d to row %d\n", p-m[r*nc+c], c, r);
				add_row_multiple_ui (&m[r*nc], &m[c*nc], p - m[r*nc+c], p, nc);
			}
		}
	}
	for (int r = 0; r < nr; r++) {
		for (int c=0; c < nc; c++) {
			printf ("%u ", m[r*nc + c]);
		}
		printf ("\n");
	}
	solution_ui *sol = new_solution_ui (p, e, nc-1);
	for (int r = 0; r < nc-1; r++) {
		sol->d[r] = m[r*nc + nc-1];
	}
	free(m);
	return sol;
}

solution *solve (matrix *m, mpz_t p, uint32_t e) {
	return NULL;
}

matelem_ui mod (int64_t x, uint64_t n) {
	if (x > 0) return x % n;
	uint64_t z = n + (x % (int64_t) n);
	if (z == n) return 0;
	return z;
}

void solution_ui_merge_ui (solution_ui *dest, solution_ui *src) {
	mpz_t tmp, m1, m2, p, q;
	mpz_inits (tmp,m1,m2,p,q,NULL);
	mpz_set_ui (p, dest->p);
	mpz_set_ui (q, src->p);
	uint64_t n = dest->p * src->p;
	mpz_gcdext (tmp, m1, m2, p, q);
	printf ("p = %llu, q = %llu; m1 = %ld, m2 = %ld\n", dest->p, src->p, mpz_get_si(m1), mpz_get_si(m2));
	for (int i=0; i < dest->size; i++) {
		//gmp_printf ("p = %Zd, q = %Zd, m1 = %Zd, m2 = %Zd\n", p, q, m1, m2);
		int64_t x = ((dest->p * (int64_t) mpz_get_si(m1) * (int64_t) src->d[i]) +
					 src->p  * (int64_t) mpz_get_si(m2) * (int64_t) dest->d[i]) ;
		dest->d[i] = mod(x, n);
	}
	dest->p *= src->p;
	mpz_clears (tmp, m1, m2, p, q, NULL);
	free_solution_ui (src);
}

void solution_merge (solution *dest, solution *src) {
	mpz_t tmp, m1, m2;
	mpz_inits (tmp,m1,m2,NULL);
	for (int i=0; i < dest->size; i++) {
		mpz_gcdext (tmp, m1, m2, dest->p, src->p);
		mpz_mul (tmp, dest->p, m1);
		mpz_mul (tmp, tmp, src->d[i]);
		mpz_mul (src->d[i], src->p, m2);
		mpz_mul (src->d[i], src->d[i], dest->d[i]);
		mpz_add (dest->d[i], tmp, src->d[i]);
	}
	mpz_mul (dest->p, dest->p, src->p);
	mpz_clears(tmp, m1, m2, NULL);
	free_solution (src);
}

solution *solve_all (matrix *m, ic_opts *opt, FILE *factors) {
	if (factors == NULL) {
	} else {
		mpz_t p;
		mpz_init(p);
		solution *prev = NULL;
		while (gmp_fscanf(factors, "%Zd\n", p) != 0) {
			solution *sol = solve(m, p, 0);
			if (prev != NULL) {
				solution_merge (prev, sol);
			} else {
				prev = sol;
			}
		}
		mpz_clear(p);
		return prev;
	}
}	
solution_ui* solve_all_ui (matrix *m, ic_opts *opt, FILE *factors) {
	if (factors == NULL) {
		uint64_t t = 2;
		uint64_t n = mpz_get_ui (opt->in.q) - 1;
		solution_ui *prev = NULL;
		while (n != 1) {
			if (n % t == 0) {
				while (n % t == 0) {
					n /= t;
				}
				solution_ui *sol = solve_ui (m, t, 0);
				print_solution_ui (sol);
				if (prev != NULL) {
					solution_ui_merge_ui (prev, sol);
				} else {
					prev = sol;
				}
			}
			t = (t == 2) ? 3 : t+2;
		}
		return prev;
	} else {
		uint64_t p;
		solution_ui *prev = NULL;
		while (fscanf(factors, "%llu\n", &p) != 0) {
			solution_ui *sol = solve_ui(m, p, 0);
			if (prev != NULL) {
				solution_ui_merge_ui (prev, sol);
			} else {
				prev = sol;
			}
		}
		return prev;
	}
}
