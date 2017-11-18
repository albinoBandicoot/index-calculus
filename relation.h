#ifndef RELATION_H
#define RELATION_H

#include <gmp.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include "core.h"

typedef struct {
	uint8_t *exp;
	uint32_t k;
	uint32_t size;
} relation;

typedef struct {
	relation **rows;
	uint32_t nrows;
	uint32_t next_slot;
} matrix;

relation *new_relation (ic_opts*);
void clear_relation (relation *r);
void print_relation (relation *r);
void print_matrix (matrix *m);
bool is_matrix_full (matrix *m);

matrix *new_matrix (int nrows);
bool add_row (matrix *m, relation *r);


#endif
