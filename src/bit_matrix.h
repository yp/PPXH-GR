/**
 *                PPXH-GR
 * A heuristic algorithm for the Pure Parsimony Xor Haplotyping problem
 *
 * Copyright (C) 2008  Yuri Pirola <yuri.pirola(-at-)gmail.com>
 *
 * Distributed under the terms of the GNU General Public License (GPL)
 *
 *
 * This file is part of PPXH-GR.
 *
 * PPXH-GR is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PPXH-GR is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PPXH-GR.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
#ifndef __BIT_MATRIX_H__
#define __BIT_MATRIX_H__

#include <stdbool.h>
#include <stddef.h>

#define _BTYPE unsigned int
#define _LBTYPE (sizeof(_BTYPE)*8)

struct _bit_mat
{
  _BTYPE* arr;
  int nrow, ncol;
  int* perm_row;
  size_t ncelle;
  size_t celle_vettore;
};

typedef struct _bit_mat * pbit_mat;

typedef void (*print_row_names)(int r);


pbit_mat BM_create(int nrow, int ncol);

void BM_destroy(pbit_mat bm);

void BM_set(pbit_mat bm, int r, int c, bool value);

void BM_set_block(pbit_mat bm, int r, int c, _BTYPE block);

void BM_copy_row(const pbit_mat dest, const int r_dest,
					  const pbit_mat src, const int r_src);

bool BM_get(pbit_mat bm, int r, int c);

_BTYPE BM_get_block(pbit_mat bm, int r, int c);

void BM_print_with_row_names(pbit_mat bm, print_row_names frow);

void BM_print(pbit_mat bm);

void BM_clear_perm_row(pbit_mat bm);

void BM_sort_row(pbit_mat bm);

pbit_mat BM_clone(pbit_mat bm);

void BM_copy(pbit_mat ris, pbit_mat bm);


#endif
