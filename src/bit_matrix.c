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
#include "bit_matrix.h"

#include "util.h"
#include "abort.h"
#include <stdio.h>
#include <limits.h>
#include <string.h>

#define _ASSERT_VALID_BM( bm )						\
  my_assert(bm!=NULL);

#define _ASSERT_VALID_POS( bm, nr, nc )				\
  my_assert(0<=(nr) && (nr)<(bm)->nrow &&					\
			0<=(nc) && (nc)<(bm)->ncol)

pbit_mat BM_create(int nrow, int ncol)
{
  my_assert(nrow>0 && ncol>0);
  pbit_mat bm= PALLOC(struct _bit_mat);
  bm->nrow= nrow;
  bm->ncol= ncol;
  bm->celle_vettore= (ncol/_LBTYPE)+1;
  bm->ncelle= (nrow*bm->celle_vettore);
  bm->arr= NPALLOC(_BTYPE, bm->ncelle);
  bm->perm_row= NPALLOC(int, nrow);
  int i;
  for (i= 0; i<bm->ncelle; ++i) {
	 bm->arr[i]= (_BTYPE)0;
  }
  for (i=0; i<nrow; ++i)
	 bm->perm_row[i]= i;
  return bm;
}

void BM_destroy(pbit_mat bm)
{
  _ASSERT_VALID_BM(bm);
  pfree(bm->arr);
  bm->arr= NULL;
  pfree(bm->perm_row);
  bm->perm_row= NULL;
  pfree(bm);
}

#define _LINT unsigned long

static void
transform_coord(pbit_mat bm, _LINT r, _LINT c, _LINT* cella, _BTYPE* mask)
{
  r= bm->perm_row[r];
  *cella= (r*bm->celle_vettore)+(c/_LBTYPE);
  *mask= 1L<<(c%_LBTYPE);
}

void BM_set(pbit_mat bm, int r, int c, bool value)
{
  _ASSERT_VALID_BM(bm);
  _ASSERT_VALID_POS(bm, r, c);
  _LINT cella;
  _BTYPE mask;
  transform_coord(bm, r, c, &cella, &mask);
  if (value) {
	 bm->arr[cella]= bm->arr[cella] | mask;
  } else {
	 bm->arr[cella]= bm->arr[cella] & (~mask);
  }
}

void BM_set_block(pbit_mat bm, int r, int c, _BTYPE block)
{
  _ASSERT_VALID_BM(bm);
  _ASSERT_VALID_POS(bm, r, c);
  my_assert(c%_LBTYPE==0);
  _LINT cella;
  _BTYPE mask;
  transform_coord(bm, r, c, &cella, &mask);
  bm->arr[cella]= block;
}

void BM_copy_row(const pbit_mat dest, const int r_dest,
					  const pbit_mat src, const int r_src)
{
  _ASSERT_VALID_BM(dest);
  _ASSERT_VALID_POS(dest, r_dest, 0);
  _ASSERT_VALID_BM(src);
  _ASSERT_VALID_POS(src, r_src, 0);
  my_assert(src->ncol==dest->ncol);
  int c;
  for (c= 0; c<src->ncol; c+= _LBTYPE)
	 BM_set_block(dest, r_dest, c, BM_get_block(src, r_src, c));
}

bool BM_get(pbit_mat bm, int r, int c)
{
  _ASSERT_VALID_BM(bm);
  _ASSERT_VALID_POS(bm, r, c);
  _LINT cella;
  _BTYPE mask;
  transform_coord(bm, r, c, &cella, &mask);
  return (bm->arr[cella] & mask)!=0;
}

_BTYPE BM_get_block(pbit_mat bm, int r, int c)
{
  _ASSERT_VALID_BM(bm);
  _ASSERT_VALID_POS(bm, r, c);
  my_assert(c%_LBTYPE==0);
  _LINT cella;
  _BTYPE mask;
  transform_coord(bm, r, c, &cella, &mask);
  return bm->arr[cella];
}

void BM_print(pbit_mat bm)
{
  _ASSERT_VALID_BM(bm);
  int r, c;
  for (r= 0; r<bm->nrow; ++r) {
	 for (c= 0; c<bm->ncol; ++c) {
		printf("%c", BM_get(bm, r, c)?'1':'0');
	 }
	 printf("\n");
  }
}


void BM_print_with_row_names(pbit_mat bm, print_row_names frow)
{
  _ASSERT_VALID_BM(bm);
  int r, c;
  for (r= 0; r<bm->nrow; ++r) {
	 frow(r);
	 for (c= 0; c<bm->ncol; ++c) {
		printf("%c", BM_get(bm, r, c)?'1':' ');
	 }
	 printf("\n");
  }
}


void BM_clear_perm_row(pbit_mat bm)
{
  _ASSERT_VALID_BM(bm);
  int i;
  for (i= 0; i<bm->nrow; ++i)
	 bm->perm_row[i]= i;
}

void BM_sort_row(pbit_mat bm)
{
  _ASSERT_VALID_BM(bm);
  BM_clear_perm_row(bm);
  int* newperm= NPALLOC(int, bm->nrow);
  int col;
  int* vn0= NPALLOC(int, bm->nrow);
  int* vn1= NPALLOC(int, bm->nrow);
  int i, j;
  memset(vn0, 0, bm->ncol*sizeof(int));
  memset(vn1, 0, bm->ncol*sizeof(int));
  for (i= 0; i<bm->nrow; ++i) {
	 for (col= 0; col<bm->ncol; col+= _LBTYPE) {
		_BTYPE block= BM_get_block(bm, i, col);
		for (j= col; j<bm->ncol; ++j, block= block>>1) {
		  if ((block & 1L)!=0)
			 ++vn1[j];
		  else
			 ++vn0[j];
		}
	 }
  }
  int n0, n1;
  for (col= bm->ncol-1; col>=0; --col) {
	 n0= vn0[col];
	 n1= n0+vn1[col];
	 for (i= bm->nrow-1; i>=0; --i) {
		if (BM_get(bm, i, col)) {
		  --n1;
		  newperm[n1]= bm->perm_row[i];
		} else {
		  --n0;
		  newperm[n0]= bm->perm_row[i];
		}
	 }
	 int * tmp_perm= bm->perm_row;
	 bm->perm_row= newperm;
	 newperm= tmp_perm;
  }
  pfree(vn0);
  pfree(vn1);
  pfree(newperm);
}



pbit_mat BM_clone(pbit_mat bm)
{
  pbit_mat ris= PALLOC(struct _bit_mat);
  ris->nrow= bm->nrow;
  ris->ncol= bm->ncol;
  ris->ncelle= bm->ncelle;
  ris->celle_vettore= bm->celle_vettore;
  ris->arr= NPALLOC(_BTYPE, bm->ncelle);
  ris->perm_row= NPALLOC(int, bm->nrow);
  memcpy(ris->arr, bm->arr, bm->ncelle*sizeof(_BTYPE));
  memcpy(ris->perm_row, bm->perm_row, bm->nrow*sizeof(int));
  return ris;
}

void BM_copy(pbit_mat ris, pbit_mat bm)
{
  memcpy(ris->arr, bm->arr, bm->ncelle*sizeof(_BTYPE));
  memcpy(ris->perm_row, bm->perm_row, bm->nrow*sizeof(int));
}



#undef _ASSERT_VALID_BM
#undef _ASSERT_VALID_POS
