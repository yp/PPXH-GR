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
#ifndef __BIT_VECTOR_H__
#define __BIT_VECTOR_H__

#include <stdbool.h>
#include <stddef.h>

#define _BTYPE unsigned int
#define _LBTYPE (sizeof(_BTYPE)*8)

struct _bit_vect
{
  _BTYPE* arr;
  int n;
  size_t ncelle;
};

typedef struct _bit_vect * pbit_vect;


pbit_vect BV_create(int n);

void BV_destroy(pbit_vect bv);

void BV_set(pbit_vect bv, int i, bool value);

void BV_set_block(pbit_vect bv, int i, _BTYPE block);

bool BV_get(pbit_vect bv, int i);

_BTYPE BV_get_block(pbit_vect bv, int i);

_BTYPE BV_get_unaligned_block(pbit_vect v, const int i, const int l);

void BV_print(pbit_vect bv);

pbit_vect BV_clone(pbit_vect bv);

void BV_copy(pbit_vect ris, pbit_vect bv);


#endif
