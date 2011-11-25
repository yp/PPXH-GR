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
#ifndef __ABORT_H__
#define __ABORT_H__

#include <stdbool.h>

#ifndef NDEBUG

void my_abort(void);

void abort_if(bool cond);

void my_assert(bool cond);

#else // NDEBUG

#define my_abort( ) ;

#define abort_if( c ) ;

#define my_assert( c ) ;


#endif // NDEBUG

#endif // __ABORT_H__
