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
#ifndef __RND_H__
#define __RND_H__

#include <stdbool.h>

struct _rnd_gen;

typedef struct _rnd_gen * prnd_gen;

prnd_gen
RG_create(void);

prnd_gen
RG_create_seed(unsigned long int);

void
RG_destroy(prnd_gen rg);

unsigned long
RG_next_int(prnd_gen rg);

unsigned long
RG_next_int_less_than(prnd_gen rg, unsigned long max);

unsigned long
RG_next_int_between(prnd_gen rg, unsigned long min, unsigned long max);

bool
RG_next_bool(prnd_gen rg);

double
RG_next_probability(prnd_gen rg);

#endif //__RND_H__
