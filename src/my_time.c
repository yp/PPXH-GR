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
#include "my_time.h"
#include "util.h"
#include <stdio.h>
#include <assert.h>

#define TOUT stdout
#define DTYPE unsigned long long

static DTYPE
diff_usec(struct timeval start, struct timeval stop)
{
  DTYPE start_usec= ((DTYPE)start.tv_sec*(DTYPE)1000000)+(DTYPE)start.tv_usec;
  DTYPE stop_usec= ((DTYPE)stop.tv_sec*(DTYPE)1000000)+(DTYPE)stop.tv_usec;
  return stop_usec-start_usec;
}

static void
print_difference(pmytime pt)
{
  DTYPE diff= diff_usec(pt->start, pt->stop);
  long double ld= ((long double)(diff)/(long double)(1000.0));
  fprintf(TOUT, "@time: %Lfmsec\n", ld);
}


pmytime
MYTIME_create(void)
{
  pmytime pt= PALLOC(struct _mytime);
  pt->active= false;
  return pt;
}

void
MYTIME_destroy(pmytime pt)
{
  assert(pt!=NULL);
  pfree(pt);
}

void
MYTIME_start(pmytime pt)
{
  assert(pt!=NULL);
  pt->active= true;
  gettimeofday(&pt->start, NULL);
}

void
MYTIME_print_reset(pmytime pt)
{
  assert(pt!=NULL);
  gettimeofday(&(pt->stop), NULL);
  pt->active= false;
  print_difference(pt);
  MYTIME_start(pt);
}

void
MYTIME_print_stop(pmytime pt)
{
  assert(pt!=NULL);
  gettimeofday(&(pt->stop), NULL);
  pt->active= false;
  print_difference(pt);
}

void
MYTIME_print_continue(pmytime pt)
{
  assert(pt!=NULL);
  gettimeofday(&(pt->stop), NULL);
  print_difference(pt);
}


#undef TOUT
#undef DTYPE
