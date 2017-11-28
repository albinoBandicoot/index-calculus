#include "matsolve.h"
int INV_CT = 0;

solution *new_solution (mpz_t p, int len) {
	solution *sol = NEW(solution);
	sol->d = malloc (sizeof(mpz_t) * len);
	mpz_init_set (sol->p, p);
	sol->size = len;
	return sol;
}

void free_solution (solution *sol) {
	mpz_clear (sol->p);
	for (int i=0; i < sol->size; i++) {
		mpz_clear (sol->d[i]);
	}
	free (sol->d);
	free (sol);
}

void print_solution (solution *sol) {
	for (int i=0; i < sol->size; i++) {
		gmp_printf ("L_%d = %Zd (mod %Zd)\n", i, sol->d[i], sol->p);
	}
}

void mul_row (mpz_t *row, mpz_t mul, mpz_t modulus, int len) {
	for (int i=0; i < len; i++) {
		mpz_mul (row[i], row[i], mul);
		mpz_mod (row[i], row[i], modulus);
	}
}

void add_row_multiple (mpz_t *dest, mpz_t *src, mpz_t mul, mpz_t modulus, int len) {
	mpz_t tmp;
	mpz_init (tmp);
	for (int i=0; i < len; i++) {
		mpz_mul(tmp, src[i], mul);
		mpz_add(dest[i], dest[i], tmp);
		mpz_mod(dest[i], dest[i], modulus);
	}
	mpz_clear(tmp);
}

void swap_rows (mpz_t *a, mpz_t *b, int len) {
	for (int i=0; i < len; i++) {
		mpz_swap (a[i], b[i]);
	}
}

#define HAS_INVERSE 1
#define NO_INVERSE 0
int generic_inversion_test (mpz_t a, mpz_t n) {
	if (mpz_cmp_ui (a, 0) == 0) return NO_INVERSE;
	mpz_t temp;
	mpz_init (temp);
	mpz_gcd (temp, a, n);
	if (mpz_cmp_ui (temp, 1) == 0) {
		mpz_clear(temp);
		return HAS_INVERSE;
	}
	mpz_clear(temp);
	return NO_INVERSE;
}

int safe_prime_inversion_test (mpz_t a, mpz_t n) {
	if ((mpz_get_ui(a) & 1) == 0) return NO_INVERSE;
	mpz_mul_ui(a, a, 2);
	if (mpz_cmp (a, n) == 0) {
		mpz_divexact_ui (a, a, 2);
		return NO_INVERSE;
	}
	mpz_divexact_ui (a, a, 2);
	return HAS_INVERSE;
}

solution *solve (matrix *mat, mpz_t n, int (*has_inverse)(mpz_t, mpz_t)) {
	gmp_printf ("Will solve matrix mod %Zd (mpz)\n", n);
	int nr = mat->nrows;
	int nc = mat->rows[0]->size + 1;
	mpz_t *m = (mpz_t *) calloc (nr*nc, sizeof(mpz_t));
	for (int r = 0; r < nr; r++) {
		for (int c=0; c < nc-1; c++) {
			mpz_init_set_ui(m[r*nc+c], mat->rows[r]->exp[c]);	// should technically be mod n, but exponents won't be larger than 2^32
		}
		mpz_init_set_ui(m[r*nc + nc-1], mat->rows[r]->k);
		mpz_mod (m[r*nc + nc-1], m[r*nc + nc-1], n);	// only necessary when n > 2^32 and k > n
	}
	
#ifdef PRINT_LARGEMATS
	for (int r = 0; r < nr; r++) {
		for (int c=0; c < nc; c++) {
			gmp_printf ("%Zd ", m[r*nc + c]);
		}
		printf ("\n");
	}
		printf ("\n");
#endif


	mpz_t inv;
	mpz_init (inv);
	for (int c=0; c < nc-1; c++) {
//		printf ("Column %d\n", c);
		if (!has_inverse (m[c*nc + c], n)) {	// non invertable (perhaps 0) on diagonal
			// find a row with an invertable entry in this column and swap with it. 
			for (int r = c+1; r < nr; r++) {
				if (has_inverse(m[r*nc + c], n)) {
					swap_rows (&m[c*nc], &m[r*nc], nc);
					break;
				}
			}
			if (mpz_cmp_ui(m[c*nc + c],0) == 0) {	// error
				printf ("  no rows below with nonzero entries!\n");
				printf ("NO INFORMATION ABOUT fb[%d]\n", c);
				continue;	// skip this column
			}
		}
		if (mpz_cmp_ui(m[c*nc + c], 1) != 0) {	// multiply this row by the multiplicative inverse of the leading entry.
			mpz_invert (inv, m[c*nc + c], n);
//			printf ("  correcting leading entry to = 1; multiplying by %llu\n", inv);
			mul_row (&m[c*nc], inv, n, nc);
		}
		// for all other rows, add a multiple to cancel out the term in this column.
		for (int r=0; r < nr; r++) {
			if (r == c) continue;
			if (mpz_cmp_ui(m[r*nc + c], 0) != 0) {
//				printf ("   adding %llu * row %d to row %d\n", p-m[r*nc+c], c, r);
				mpz_sub (inv, n, m[r*nc+c]);
				add_row_multiple (&m[r*nc], &m[c*nc], inv, n, nc);
			}
		}
	}
#ifdef PRINT_LARGEMATS
	for (int r = 0; r < nr; r++) {
		for (int c=0; c < nc; c++) {
			gmp_printf ("%Zd ", m[r*nc + c]);
		}
		printf ("\n");
	}
#endif
	solution *sol = new_solution (n, nc-1);
	for (int r = 0; r < nc-1; r++) {
		mpz_init_set (sol->d[r], m[r*nc + nc-1]);
	}
	mpz_clear (inv);
	for (int r=0; r < nr-1; r++) {
		for (int c=0; c < nc-1; c++) {
			mpz_clear (m[r*nc + c]);
		}
	}
	free(m);
	return sol;
}
