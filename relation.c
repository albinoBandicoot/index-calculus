#include <pthread.h>
#include "relation.h"
#include "core.h"

relation *new_relation (ic_opts *opt) {
	relation *res = (relation *) malloc(sizeof(relation));
	res->size = opt->fblen;
	res->k = 0;
	res->exp = (uint8_t *) calloc(1, opt->fblen);
	return res;
}

void clear_relation (relation *r) {
	memset (r->exp, 0, r->size);
}

void print_relation (relation *r) {
	for (int i=0; i < r->size; i++) {
		printf ("%d ", r->exp[i]);
	}
	printf ("| %d\n", r->k);
}

void print_matrix (matrix *m) {
	for (int i=0; i < m->nrows; i++) {
		print_relation (m->rows[i]);
	}
}

matrix *new_matrix (int nrows) {
	matrix *m = NEW(matrix);
	m->nrows = nrows;
	m->next_slot = 0;
	m->rows = (relation **) malloc(nrows * sizeof(relation*));
	return m;
}

bool is_matrix_full (matrix *m) {
	return m->next_slot == m->nrows;
}

pthread_mutex_t MAT_MUTEX = PTHREAD_MUTEX_INITIALIZER;

// returns whether the matrix is full
bool add_row (matrix *m, relation *r) {
	pthread_mutex_lock (&MAT_MUTEX);
	if (m->next_slot == m->nrows) {
		pthread_mutex_unlock(&MAT_MUTEX);
		return true;
	}
	m->rows[m->next_slot] = r;
	m->next_slot++;
	bool matrix_full = is_matrix_full(m);
	pthread_mutex_unlock (&MAT_MUTEX);
	return matrix_full;
}
