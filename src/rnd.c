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
#include "rnd.h"

#include "util.h"
#include "abort.h"

#include <gsl/gsl_rng.h>
#include <time.h>

struct _rnd_gen
{
  gsl_rng* rng;
};

prnd_gen
RG_create(void)
{
  prnd_gen pgen= PALLOC(struct _rnd_gen);
  pgen->rng= gsl_rng_alloc(gsl_rng_mt19937);
  abort_if(pgen->rng==NULL);
  gsl_rng_set(pgen->rng, (unsigned long int)time(NULL));
  gsl_rng_set(pgen->rng, (unsigned long int)1981);
  return pgen;
}

prnd_gen
RG_create_seed(unsigned long int seed)
{
  prnd_gen pgen= PALLOC(struct _rnd_gen);
  pgen->rng= gsl_rng_alloc(gsl_rng_mt19937);
  abort_if(pgen->rng==NULL);
  gsl_rng_set(pgen->rng, seed);
  return pgen;
}

void
RG_destroy(prnd_gen rg)
{
  abort_if(rg==NULL);
  gsl_rng_free(rg->rng);
  pfree(rg);
}

unsigned long
RG_next_int(prnd_gen rg)
{
  abort_if(rg==NULL);
  return gsl_rng_get(rg->rng);
}


unsigned long
RG_next_int_less_than(prnd_gen rg, unsigned long max)
{
  abort_if(rg==NULL);
  my_assert(max>0);
  return gsl_rng_uniform_int(rg->rng, max);
}

unsigned long
RG_next_int_between(prnd_gen rg, unsigned long min, unsigned long max)
{
  abort_if(rg==NULL);
  my_assert(min<max);
  return gsl_rng_uniform_int(rg->rng, max-min)+min;
}

bool
RG_next_bool(prnd_gen rg)
{
  abort_if(rg==NULL);
  return gsl_rng_uniform_int(rg->rng, 2)==1;
}

double
RG_next_probability(prnd_gen rg)
{
  abort_if(rg==NULL);
  return gsl_rng_uniform(rg->rng);
}

