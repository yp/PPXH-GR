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
#include "util.h"

#include <stdio.h>
#include <stdlib.h>
#include "abort.h"

void* palloc(size_t size) {
  void* p= malloc(size);
  if (p==NULL) {
	 fprintf(stderr, "Allocation memory error. "
				"Trying to allocate %zu bytes.", size);
	 my_abort();
  }
  return p;
}

char* c_palloc(size_t size)
{
  return (char*)palloc(size*sizeof(char));
}

void pfree(void* p)
{
  if (p==NULL) {
	 fprintf(stderr, "Freeing a NULL pointer is not permitted.");
	 my_abort();
  }
  free(p);
}

void noop_free(void* el) {
  if (el!=NULL) el= NULL; // Istruzione inutile -> to make compiler happy
}

