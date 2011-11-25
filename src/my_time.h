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
#ifndef __MYTIME_H__
#define __MYTIME_H__

#include <stdbool.h>
#include <time.h>
#include <sys/time.h>


struct _mytime
{
  struct timeval start;
  struct timeval stop;
  bool active;
};

typedef struct _mytime* pmytime;

pmytime
MYTIME_create(void);

void
MYTIME_destroy(pmytime pt);

void
MYTIME_start(pmytime pt);

void
MYTIME_print_reset(pmytime pt);

void
MYTIME_print_stop(pmytime pt);

void
MYTIME_print_continue(pmytime pt);

#endif // __MYTIME_H__
