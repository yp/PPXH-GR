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
#include "bit_vector.h"

#include "util.h"
#include "abort.h"
#include <stdio.h>
#include <limits.h>
#include <string.h>

#define _ASSERT_VALID_BV( bv )						\
  my_assert(bv!=NULL);

#define _ASSERT_VALID_POS( bv, i )				\
  my_assert(0<=(i) && (i)<(bv)->n)

pbit_vect BV_create(int n)
{
  my_assert(n>0);
  pbit_vect bv= PALLOC(struct _bit_vect);
  bv->n= n;
  bv->ncelle= (n/_LBTYPE)+1;
  bv->arr= NPALLOC(_BTYPE, bv->ncelle);
  int i;
  for (i= 0; i<bv->ncelle; ++i) {
	 bv->arr[i]= (_BTYPE)0;
  }
  return bv;
}

void BV_destroy(pbit_vect bv)
{
  _ASSERT_VALID_BV(bv);
  pfree(bv->arr);
  bv->arr= NULL;
  pfree(bv);
}

#define _LINT unsigned long

static void
transform_coord(_LINT i, _LINT* cella, _BTYPE* mask)
{
  *cella= i/_LBTYPE;
  *mask= 1L<<(i%_LBTYPE);
}

void BV_set(pbit_vect bv, int i, bool value)
{
  _ASSERT_VALID_BV(bv);
  _ASSERT_VALID_POS(bv, i);
  _LINT cella;
  _BTYPE mask;
  transform_coord(i, &cella, &mask);
  if (value) {
	 bv->arr[cella]= bv->arr[cella] | mask;
  } else {
	 bv->arr[cella]= bv->arr[cella] & (~mask);
  }
}

void BV_set_block(pbit_vect bv, int i, _BTYPE block)
{
  _ASSERT_VALID_BV(bv);
  _ASSERT_VALID_POS(bv, i);
  my_assert(i%_LBTYPE==0);
  _LINT cella;
  _BTYPE mask;
  transform_coord(i, &cella, &mask);
  bv->arr[cella]= block;
}

bool BV_get(pbit_vect bv, int i)
{
  _ASSERT_VALID_BV(bv);
  _ASSERT_VALID_POS(bv, i);
  _LINT cella;
  _BTYPE mask;
  transform_coord(i, &cella, &mask);
  return (bv->arr[cella] & mask)!=0;
}

_BTYPE BV_get_block(pbit_vect bv, int i)
{
  _ASSERT_VALID_BV(bv);
  _ASSERT_VALID_POS(bv, i);
  my_assert(i%_LBTYPE==0);
  _LINT cella;
  _BTYPE mask;
  transform_coord(i, &cella, &mask);
  return bv->arr[cella];
}

_BTYPE
BV_get_unaligned_block(pbit_vect v, const int i, const int l)
{
  _ASSERT_VALID_BV(v);
  _ASSERT_VALID_POS(v, i);
  my_assert(l<=_LBTYPE);
  _LINT pb1= (i/_LBTYPE)*_LBTYPE;
  _LINT pb2= pb1+_LBTYPE;
  _LINT mod= i%_LBTYPE;
  _BTYPE br, tmp1, tmp2, mask;
  if (mod==0)
	 br= BV_get_block(v, i);
  else {
	 tmp1= BV_get_block(v, pb1);
	 if ((pb1+l<pb1)||(pb2>v->n))
		tmp2= 0L;
	 else
		tmp2= BV_get_block(v, pb2);
//	 print_block("tmp1 ", tmp1, _LBTYPE);
//	 print_block("tmp2 ", tmp2, _LBTYPE);
	 br= (tmp1>>mod) | (tmp2<<(_LBTYPE-mod));
  }
//  print_block("br   ", br, l);
  int j;
  mask= 0L;
  for (j= 0; j<l; ++j) {
	 mask= mask << 1;
	 mask= mask | 1L;
  }
//  printf("L= %d\n", l);
//  print_block("mask ", mask , _LBTYPE);
  br= br & mask;
//  print_block("br   ", br, _LBTYPE);
  return br;
}




void BV_print(pbit_vect bv)
{
  _ASSERT_VALID_BV(bv);
  int i;
  for (i= 0; i<bv->n; ++i) {
	 printf("%c", BV_get(bv, i)?'1':'0');
  }
  printf("\n");
}


pbit_vect BV_clone(pbit_vect bv)
{
  pbit_vect ris= PALLOC(struct _bit_vect);
  ris->n= bv->n;
  ris->ncelle= bv->ncelle;
  ris->arr= NPALLOC(_BTYPE, bv->ncelle);
  memcpy(ris->arr, bv->arr, bv->ncelle*sizeof(_BTYPE));
  return ris;
}

void BV_copy(pbit_vect ris, pbit_vect bv)
{
  memcpy(ris->arr, bv->arr, bv->ncelle*sizeof(_BTYPE));
}

#undef _ASSERT_VALID_BV
#undef _ASSERT_VALID_POS
