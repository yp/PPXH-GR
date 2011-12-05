/**
* PPXH-GR
* A heuristic algorithm for the Pure Parsimony Xor Haplotyping problem
*
* Copyright (C) 2008 Yuri Pirola <yuri.pirola(-at-)gmail.com>
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
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with PPXH-GR. If not, see <http://www.gnu.org/licenses/>.
*
**/

#define _GNU_SOURCE

#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

// getopt
#include <getopt.h>

#include "abort.h"
#include "bit_matrix.h"
#include "bit_vector.h"
#include "my_time.h"
#include "rnd.h"
#include "util.h"

typedef pbit_mat xgen_mat;

typedef pbit_mat hap_mat;

struct _bmatrix {
  bool** mat;
  int* perm_row;
  int* perm_col;
  int* name_row;
  int* name_col;
  int ncol;
  int nrow;
};

typedef struct _bmatrix* pbmatrix;


pbmatrix
bmatrix_create(const int nrow, const int ncol) {
  pbmatrix pb= PALLOC(struct _bmatrix);
  pb->ncol= ncol;
  pb->nrow= nrow;
  pb->mat= NPALLOC(bool*, nrow);
  for (int i= 0; i<nrow; ++i) {
	 pb->mat[i]= NPALLOC(bool, ncol);
	 for (int j= 0; j<ncol; ++j) {
		pb->mat[i][j]= false;
	 }
  }
  pb->perm_row= NPALLOC(int, nrow);
  pb->name_row= NPALLOC(int, nrow);
  for (int i= 0; i<nrow; ++i) {
	 pb->perm_row[i]= i;
	 pb->name_row[i]= i;
  }
  pb->perm_col= NPALLOC(int, ncol);
  pb->name_col= NPALLOC(int, ncol);
  for (int i= 0; i<ncol; ++i) {
	 pb->perm_col[i]= i;
	 pb->name_col[i]= i;
  }
  return pb;
}

void
bmatrix_destroy(pbmatrix pb) {
  for (int i= 0; i<pb->nrow; ++i) {
	 pfree(pb->mat[i]);
  }
  pfree(pb->mat);
  pfree(pb->perm_col);
  pfree(pb->name_col);
  pfree(pb->perm_row);
  pfree(pb->name_row);
  pfree(pb);

}

bool
bget(pbmatrix pb, int r, int c) {
  assert(pb!=NULL);
  assert(0<=r&&r<pb->nrow);
  assert(0<=c&&c<pb->ncol);

  return pb->mat[pb->perm_row[r]][pb->perm_col[c]];
}

void
bset(pbmatrix pb, int r, int c, bool v) {
  assert(pb!=NULL);
  assert(0<=r&&r<pb->nrow);
  assert(0<=c&&c<pb->ncol);

  pb->mat[pb->perm_row[r]][pb->perm_col[c]]= v;
}


static void
bmatrix_basic_print_row(FILE* f, char* D, pbmatrix pb, int r) {
  for (int i= 0; i< pb->ncol; ++i) {
	 fprintf(f, D, bget(pb, r, i)?1:0);
  }
  fprintf(f, "\n");
}


void
bmatrix_print_row(FILE* f, pbmatrix pb, int r) {
  bmatrix_basic_print_row(f, "%3d", pb, r);
}

void
bmatrix_basic_print(FILE* f, pbmatrix pb) {
  for (int i= 0; i< pb->nrow; ++i) {
	 bmatrix_basic_print_row(f, "%d", pb, i);
  }
}

void
bmatrix_print(FILE* f, pbmatrix pb) {
  fprintf(f, "****");
  for (int i= 0; i< pb->ncol; ++i) {
	 fprintf(f, "%3d", pb->name_col[i]);
  }
  fprintf(f, "\n");
  fprintf(f, "****");
  for (int i= 0; i< pb->ncol; ++i) {
	 fprintf(f, "---");
  }
  fprintf(f, "--\n");
  for (int i= 0; i< pb->nrow; ++i) {
	 fprintf(f, "%3d|",pb->name_row[i]);
	 bmatrix_print_row(f, pb, i);
  }
}

pbmatrix
bmatrix_copy(pbmatrix src) {
  pbmatrix dest= bmatrix_create(src->nrow, src->ncol);

  for (int r= 0; r<src->nrow; ++r) {
//	 dest->perm_row[r]= src->perm_row[r];
	 dest->name_row[r]= src->name_row[r];
  }
  for (int c= 0; c<src->ncol; ++c) {
//	 dest->perm_col[c]= src->perm_col[c];
	 dest->name_col[c]= src->name_col[c];
  }

  for (int r= 0; r<src->nrow; ++r) {
	 for (int c= 0; c<src->ncol; ++c) {
		bset(dest, r, c, bget(src, r, c));
	 }
  }

  return dest;
}

void
bmatrix_sum_row(pbmatrix pb, int rris, int r1, int r2) {
  for (int i= 0; i<pb->ncol; ++i) {
	 bset(pb, rris, i, bget(pb, r1, i)!=bget(pb, r2, i));
  }
}


#define SWAP( a, b ) \
  {						\
  int swap_temp= a;	\
  a= b;					\
  b= swap_temp;		\
  }

#define min( a, b ) (((a)<(b))?(a):(b))
#define max( a, b ) (((a)>=(b))?(a):(b))

// gauss modifica la matrice!!
int
bmatrix_gauss(pbmatrix pb) {
  int rango= -1;
  int limite= min(pb->ncol, pb->nrow);
// Rendo triangolare superiore
  for (int r= 0; r<limite; ++r) {
// assicuro un 1 in [r,r]
	 if (!bget(pb, r, r)) {
// cerco la prima cella con un 1 nella restante parte della matrice
		int r1= r, c= r;
		bool one_found= false;
		for (r1= r; r1<pb->nrow && !one_found; ++r1) {
		  for (c= r; c<pb->ncol && !one_found; ++c) {
// do nothing
			 one_found= one_found || bget(pb, r1, c);
		  }
		}
		if (one_found) {
// scambio le colonne e le righe
		  SWAP(pb->perm_col[r], pb->perm_col[c-1]);
		  SWAP(pb->name_col[r], pb->name_col[c-1]);
		  SWAP(pb->perm_row[r], pb->perm_row[r1-1]);
		  SWAP(pb->name_row[r], pb->name_row[r1-1]);
		  assert(bget(pb, r, r));
		} else {
		  rango= r;
		  break;
		}
	 }
// azzero la parte restante della colonna
	 for (int r2= r+1; r2<pb->nrow; ++r2) {
		if (bget(pb, r2, r)) {
		  bmatrix_sum_row(pb, r2, r, r2);
		}
	 }
  }
  if (rango==-1)
	 rango= limite;
// Rendo "matrice identica"
  for (int r= rango-1; r>0; --r) {
	 for (int r2= 0; r2<r; ++r2) {
		if (bget(pb, r2, r)) {
		  bmatrix_sum_row(pb, r2, r, r2);
		}
	 }
  }
  return rango;
}

pbmatrix
bmatrix_transpose_first_m_column(pbmatrix src, int m) {
  assert(src!=NULL);
  assert(0<m&&m<=src->ncol);
  pbmatrix dest= bmatrix_create(m, src->nrow);

  for (int r= 0; r<src->nrow; ++r) {
//	 dest->perm_col[r]= src->perm_row[r];
	 dest->name_col[r]= src->name_row[r];
  }
  for (int c= 0; c<m; ++c) {
//	 dest->perm_row[c]= src->perm_col[c];
	 dest->name_row[c]= src->name_col[c];
  }

  for (int r= 0; r<src->nrow; ++r) {
	 for (int c= 0; c<m; ++c) {
		bset(dest, c, r, bget(src, r, c));
	 }
  }

  return dest;
}


pbmatrix
XGM_to_bmatrix(xgen_mat X) {
  pbmatrix ris= bmatrix_create(X->nrow, X->ncol);
  for (int r= 0; r<X->nrow; ++r) {
	 for (int c= 0; c<X->ncol; ++c) {
		bset(ris, r, c, BM_get(X, r, c));
	 }
  }
  return ris;
}

typedef struct _row {
  int v1;
  int v2;
  int uniq_label;
  bool used;
  char* label;
} srow;

typedef struct _gr* pgr;


static void
execute_GR(pbmatrix Xt, bool * gen_gr) {
  bool gr_found= false;
// Intervallo di colonne su cui cercare la GR [base_c, end_c)
  int base_c= Xt->nrow;
  int next_g= base_c+1;
  for (int g= 0; g< Xt->ncol; ++g) {
	 gen_gr[g]= false;
  }
  int last_g= base_c;
  gen_gr[last_g]= true;
  do {
	 FILE* F= fopen("hypergraph.txt", "w");
	 for (int c= base_c; c<Xt->ncol; ++c) {
		if (gen_gr[c]) {
		  for (int r= 0; r < Xt->nrow; ++r) {
			 if (r>0)  fprintf(F, " ");
			 fprintf(F, bget(Xt, r, c)?"1":"0");
		  }
		  fprintf(F, "\n");
		}
	 }
	 fclose(F);
	 int exit_result= system("./exec-gr.sh");
	 gr_found= exit_result==0;
	 if (!gr_found) {
		gen_gr[last_g]= false;
		fprintf(stderr, "GR not found %d becomes FALSE.\n", last_g);
	 }
	 if (next_g<Xt->ncol) {
//		fprintf(stderr, "# Set %d to TRUE.\n", next_g);
		gen_gr[next_g]= true;
		last_g= next_g;
	 }
	 ++next_g;
  } while(next_g<=Xt->ncol);
  FILE* F= fopen("hypergraph.txt", "w");
  for (int c= base_c; c<Xt->ncol; ++c) {
	 if (gen_gr[c]) {
		for (int r= 0; r < Xt->nrow; ++r) {
		  if (r>0)  fprintf(F, " ");
		  fprintf(F, bget(Xt, r, c)?"1":"0");
		}
		fprintf(F, "\n");
	 }
  }
  fclose(F);
  int exit_result= system("./exec-gr.sh");
  assert(exit_result==0);
}


int
cerca_riga_genotipo(pbmatrix X, int g) {
  int real_r= -1;
  for (int r1= 0; r1<X->nrow && real_r==-1; ++r1) {
	 if (X->name_row[r1]==g)
		real_r= r1;
  }
  return real_r;
}

static void
read_GR(pbmatrix X,
		  pbmatrix Xt,
		  pbmatrix H,
		  int * nhap,
		  bool * solved_gen) {
  srow* realiz= NPALLOC(srow, Xt->nrow+1);
  for (int temp= 0; temp<Xt->nrow+1; ++temp)
	 realiz[temp].used= false;
  FILE* F= fopen("realization.txt", "r");
  char* buffer;
  size_t lbuff= 0;
  int riga= 0, maxv= -1;
  fprintf(stderr, "Realization (as read)\n");
  while (!feof(F)) {
	 int bread=getline(&buffer, &lbuff, F);
	 if (bread<=0)
		continue;
	 fprintf(stderr, "%s", buffer);
// legge il primo numero
	 int i= 0;
	 while(isdigit(buffer[i]))
		++i;
	 buffer[i]= '\0';
	 realiz[riga].v1= atoi(buffer);
	 int pi= i+1;
	 i= pi;
	 while(isdigit(buffer[i]))
		++i;
	 buffer[i]= '\0';
	 realiz[riga].v2= atoi(buffer+pi);
	 pi= i+1;
	 i= pi;
	 realiz[riga].label= NPALLOC(char, strlen(buffer+pi));
	 while(isdigit(buffer[i])||buffer[i]==',') {
		realiz[riga].label[i-pi]= buffer[i];
		++i;
	 }
	 realiz[riga].label[i-pi]= '\0';
	 maxv= max(maxv, realiz[riga].v1);
	 maxv= max(maxv, realiz[riga].v2);
	 ++riga;
  }
  free(buffer);
// suddivide le label
  int righe_lette= riga;
  for (int r= 0; r<righe_lette; ++r) {
	 realiz[r].used= false;
	 int v1= realiz[r].v1;
	 int i= 0;
	 int basep= 0;
	 do {
		while(isdigit(realiz[r].label[i]))
		  ++i;
		if (realiz[r].label[i]==',') {
		  realiz[r].label[i]= '\0';
		  ++maxv;
		  realiz[riga].v1= v1;
		  realiz[riga].v2= maxv;
		  realiz[riga].used= false;
		  v1= maxv;
		  realiz[riga].uniq_label= atoi(realiz[r].label+basep);
		  ++riga;
		  ++i;
		  basep= i;
		} else {
		  assert(realiz[r].label[i]=='\0');
		  realiz[r].v1= v1;
		  realiz[r].uniq_label= atoi(realiz[r].label+basep);
		}
	 } while (realiz[r].label[i]!='\0');
  }
  fprintf(stderr, "Realization (as interpreted)\n");
  for (int r= 0; r<riga; ++r) {
	 fprintf(stderr, "n%3d--n%3d   (g%3d)\n",
				realiz[r].v1, realiz[r].v2,
				Xt->name_col[realiz[r].uniq_label]);
  }
// estraggo gli aplotipi
  int * key= NPALLOC(int, riga+1);
  int * value= NPALLOC(int, riga+1);
  bool * visited= NPALLOC(bool, riga+1);
  for (int r= 0; r<=riga; ++r) {
	 key[r]= 0;
	 value[r]= 0;
	 visited[r]= false;
  }
  int scorr= 0;
  while (scorr<riga) {
	 if (!realiz[scorr].used) {
		int npila= 1;
		key[0]= realiz[scorr].v1;
		value[0]= 0;
		visited[realiz[scorr].v1]= true;
		realiz[scorr].used= true;
		fprintf(stderr, "Haplotype matrix construction "
				  "starting from %d\n", realiz[scorr].v1);
		while (npila>0) {
		  --npila;
		  int ckey= key[npila];
		  int cvalue= value[npila];
		  for (int r= 0; r<riga; ++r) {
			 int h1= -1, h2= -1, g= -1;
			 if (realiz[r].v1==ckey || realiz[r].v2==ckey) {
				realiz[r].used= true;
			 }
			 if (realiz[r].v1==ckey && !visited[realiz[r].v2]) {
				h1= cvalue;
				h2= realiz[r].v2;
				g= realiz[r].uniq_label;
			 } else if (realiz[r].v2==ckey && !visited[realiz[r].v1]) {
				h1= cvalue;
				h2= realiz[r].v1;
				g= realiz[r].uniq_label;
			 }
			 if (g!=-1 && solved_gen[Xt->name_col[g]]) {
				fprintf(stderr, "h%d = g%d xor h%d\n",
						  *nhap, Xt->name_col[g], h1);
				visited[h2]= true;
// cerco la riga del genotipo etichettato con Xt->name_col[g]
				int real_r= cerca_riga_genotipo(X, Xt->name_col[g]);
				assert(real_r>=0);
				for (int c= 0; c<X->ncol; ++c) {
				  bset(H, (*nhap), c, bget(X, real_r, c)!=bget(H, h1, c));
				}
				key[npila]= h2;
				value[npila]= (*nhap);
				++npila;
				++(*nhap);
			 } else if (g!=-1) {
				fprintf(stderr, "genotype %d is not solved. ignoring h%d e h%d\n",
						  Xt->name_col[g], h1, h2);
			 }
		  }
		}
	 }
	 ++scorr;
  }
  pfree(key);
  pfree(value);
  pfree(visited);
}

bool
is_solved_by(pbmatrix X, int g,
				 pbmatrix H, int h1, int h2) {
  bool solved= true;
  for (int c= 0; c<X->ncol && solved; ++c) {
	 solved= solved && bget(X, g, c) == (bget(H, h1, c)!=bget(H, h2, c));
  }
  return solved;
}


void
get_biggest_GR(pbmatrix X,
					pbmatrix Xt,
					pbmatrix H,
					int * nhap,
					bool* solved_gen,
					int * nsolved_gen) {
#ifndef NDEBUG
  {
// le prime nrow colonne di Xt sono la matrice identica
	 assert(Xt->nrow < Xt->ncol);
// i genotipi che intestano le colonne non sono risolti
	 for (int c= 0; c<Xt->ncol; ++c) {
		assert(!solved_gen[Xt->name_col[c]]);
	 }
// il numero di genotipi risolti coincide con il numero
// di true nell'array
	 int real_nsolved_gen= 0;
	 for (int g= 0; g<X->nrow; ++g) {
		if (solved_gen[g])
		  ++real_nsolved_gen;
	 }
	 assert(real_nsolved_gen==*nsolved_gen);
  }
#endif
  bool* gen_gr= NPALLOC(bool, Xt->ncol);
  execute_GR(Xt, gen_gr);
  fprintf(stderr, "Realization of columns:\n");
  for (int c= 0; c<Xt->ncol; ++c) {
	 fprintf(stderr, "%d (g%d) -> %d\n", c, Xt->name_col[c], gen_gr[c]?1:0);
  }
// determino i genotipi risolti
  for (int r= 0; r<Xt->nrow; ++r) {
	 bool used_r= false;
	 for (int c= Xt->nrow; c<Xt->ncol && !used_r; ++c)
		used_r= used_r || (gen_gr[c] && bget(Xt, r, c));
	 if (used_r) {
		assert(!solved_gen[Xt->name_col[r]]);
		solved_gen[Xt->name_col[r]]= true;
		++(*nsolved_gen);
	 }
  }
// aggiorno gli altri genotipi risolti
  for (int c= Xt->nrow; c<Xt->ncol; ++c) {
	 if (gen_gr[c]) {
		assert(!solved_gen[Xt->name_col[c]]);
		solved_gen[Xt->name_col[c]]= true;
		++(*nsolved_gen);
	 }
  }
#ifndef NDEBUG
// stampa genotipi risolti
  fprintf(stderr, "Genotipi risolti\n");
  for (int g= 0; g<X->nrow; ++g) {
	 fprintf(stderr, "g%3d) %d\n", g, solved_gen[g]?1:0);
  }
  fprintf(stderr, "Totale genotipi risolti: %d\n", *nsolved_gen);
#endif
// leggo la realizzazione
  read_GR(X, Xt, H, nhap, solved_gen);
// aggiungo genotipi risolti indirettamente
  for (int g= 0; g<X->nrow; ++g) {
	 if (!solved_gen[g]) {
		int real_r= cerca_riga_genotipo(X, g);
		bool s= false;
		for (int h1= 0; h1<*nhap && !s; ++h1) {
		  for (int h2= h1+1; h2<*nhap && !s; ++h2) {
			 s= s||is_solved_by(X, real_r, H, h1, h2);
			 if (s) {
				fprintf(stderr, "Genotype %d indirectly solved by h%d and h%d\n",
						  g, h1, h2);
			 }
		  }
		}
		if (s) {
		  fprintf(stderr, "Genotype %d indirectly solved\n", g);
		  solved_gen[g]= true;
		  ++(*nsolved_gen);
		}
	 }
  }
}


void
solve_genotypes(pbmatrix X,
					 pbmatrix Xp,
					 pbmatrix H,
					 int * nhap,
					 bool* solved_gen,
					 int * nsolved_gen) {
  pbmatrix Xpc= bmatrix_copy(Xp);
  fprintf(stderr, "Solving genotypes...\n");
//  bmatrix_print(stderr, Xpc);
  int rango= bmatrix_gauss(Xpc);
  fprintf(stderr, "Matrix rank %d\n", rango);
//  bmatrix_print(stderr, Xpc);
// Copio la permutazione delle colonne su Xp
  for (int i= 0; i<Xp->ncol; ++i) {
	 Xp->perm_col[i]= Xpc->perm_col[i];
	 Xp->name_col[i]= Xpc->name_col[i];
  }
  if (rango<Xp->nrow) {
	 pbmatrix Xt= bmatrix_transpose_first_m_column(Xp, rango);
	 fprintf(stderr, "Applying Gauss algorithm on X^t...\n");
	 bmatrix_gauss(Xt);
//	 bmatrix_print(stderr, Xt);
	 get_biggest_GR(X, Xt, H, nhap, solved_gen, nsolved_gen);
//	 fprintf(stderr, "Partial H\n");
//	 bmatrix_print(stderr, H);
	 if (X->nrow>*nsolved_gen) {
		fprintf(stderr, "Genotype to solve %d\n", X->nrow-*nsolved_gen);
		pbmatrix Xpp= bmatrix_create(X->nrow-*nsolved_gen, X->ncol);
		for (int g= 0, new_g= 0; g<Xp->nrow; ++g) {
		  if (!solved_gen[Xp->name_row[g]]) {
// Copio il genotipo non risolto
			 int real_r= cerca_riga_genotipo(X, Xp->name_row[g]);
			 fprintf(stderr, "Genotype %d in row %d not yet solved\n",
						Xp->name_row[g], real_r);
			 for (int c= 0; c<X->ncol; ++c) {
				bset(Xpp, new_g, c, bget(X, real_r, c));
			 }
			 Xpp->name_row[new_g]= Xp->name_row[g];
			 ++new_g;
		  }
		}
		solve_genotypes(X, Xpp, H, nhap, solved_gen, nsolved_gen);
	 }
  } else {
// Aggiungo i genotipi rimanenti come aplotipi
	 for (int i= 0; i<Xp->nrow; ++i) {
		int real_r= cerca_riga_genotipo(X, Xp->name_row[i]);
		fprintf(stderr, "Add genotype %d from row %d\n",
				  Xp->name_row[i], real_r);
		for (int c= 0; c<X->ncol; ++c) {
		  bset(H, (*nhap), c, bget(X, real_r, c));
		}
		++(*nhap);
	 }
  }
}

static void
shuffle_rows(pbmatrix X, pbmatrix Xp) {
//  unsigned long int seed=  1228234799;//  1228229164;
  unsigned long int seed= (unsigned long int)time(NULL);
  fprintf(stderr, "shuffle seed %lu\n", seed);
  prnd_gen rg= RG_create_seed(seed);
  for (int r= 0; r<X->nrow; ++r) {
	 int new_r= RG_next_int_less_than(rg, X->nrow);
	 if (r!=new_r) {
		SWAP(X->name_row[r], X->name_row[new_r]);
		SWAP(X->perm_row[r], X->perm_row[new_r]);
		SWAP(Xp->name_row[r], Xp->name_row[new_r]);
		SWAP(Xp->perm_row[r], Xp->perm_row[new_r]);
	 }
  }
  RG_destroy(rg);
}


static void
run(xgen_mat Xgm) {
  pbmatrix X= XGM_to_bmatrix(Xgm);
  pbmatrix Xp= XGM_to_bmatrix(Xgm);
  shuffle_rows(X, Xp);
  pbmatrix H= bmatrix_create(Xp->nrow+1, Xp->ncol);
  int nhap= 1;
  bool* solved_gen= NPALLOC(bool, X->nrow);
  for (int g= 0; g<X->nrow; ++g) {
	 solved_gen[g]= false;
  }
  int nsolved_gen= 0;
  solve_genotypes(X, Xp, H, &nhap, solved_gen, &nsolved_gen);
//  fprintf(stderr, "Genotype Matrix\n");
//  bmatrix_print(stderr, X);
//  fprintf(stderr, "Haplotype Matrix\n");
//  bmatrix_print(stderr, H);
  printf("==Genotype Matrix\n");
  bmatrix_basic_print(stdout, X);
  printf("==END Genotype Matrix\n");
  printf("==Haplotype Matrix\n");
  bmatrix_basic_print(stdout, H);
  printf("==END Haplotype Matrix\n");
  printf("n. of genotypes: %d\n", X->nrow);
  printf("n. of snps: %d\n", X->ncol);
  printf("n. of distinct haplotypes: %d\n", nhap);
}

struct nodo {
  char* line;
  struct nodo* next;
};

static
ssize_t my_getline(char **lineptr, size_t *n, FILE *stream) {
  ssize_t ris= getline(lineptr, n, stream);
  while (ris>0 && !isprint((*lineptr)[ris-1])) {
	 (*lineptr)[ris-1]= '\0';
	 --ris;
  }
  return ris;
}

static
xgen_mat read_genmat_from_file(char * fc) {
  FILE *f;
  if (strcmp(fc, "-")==0) {
	 f= stdin;
  } else {
	 f= fopen(fc, "r");
  }
  int nrow= 0;
  int ncol= 0;
  struct nodo* first= PALLOC(struct nodo);
  first->line= NULL;
  first->next= NULL;
  struct nodo* cur= first;
  while (!feof(f)) {
	 char* buffer= NULL;
	 size_t size= 0;
	 int br= my_getline(&buffer, &size, f);
//	 printf("%s\n", buffer);
	 if (br==0 || strlen(buffer)==0) {
		free(buffer);
		continue;
	 }
	 if (strlen(buffer)!=ncol) {
		if (ncol==0) {
		  ncol= strlen(buffer);
		} else {
		  fprintf(stderr,
					 "Error in reading file!\n"
					 "The previous rows have %d columns, while the current row has %zu columns.\n"
					 "Aborting...", ncol, strlen(buffer));
		  my_abort();
		  exit(1);
		}
	 }
	 cur->next= PALLOC(struct nodo);
	 cur->next->line= buffer;
	 cur->next->next= NULL;
	 cur= cur->next;
	 nrow++;
  }
  printf("Read %d genotypes on %d snps\n", nrow, ncol);
  xgen_mat X= BM_create(nrow, ncol);
  cur= first;
  int ng= 0;
  while (cur->next!=NULL) {
	 cur= cur->next;
	 int ns= 0;
	 while (ns<ncol) {
		BM_set(X, ng, ns, cur->line[ns]=='2');
		++ns;
	 }
	 ++ng;
  }
  return X;
}


int main(int argc, char ** argv) {
  int nrun= 10;
  char* file= NULL;
  int c;
// leggo i parametri
  while ((c = getopt (argc, argv, "r:f:")) != -1) {
	 switch (c) {
		case 'r':
		  nrun = atoi(optarg);
		  break;
		case 'f':
		  file = optarg;
		  break;
		case '?':
		  if (isprint (optopt))
			 fprintf(stderr, "Unknown option `-%c'.\n", optopt);
		  else
			 fprintf(stderr, "Unknown character `\\x%x'.\n", optopt);
		  return 1;
		default:
		  my_abort();
		  exit(1);
	 }
  }

  xgen_mat Xgm;
  if (file!=NULL) {
	 Xgm= read_genmat_from_file(file);
  } else {
	 fprintf(stderr, "File not specified!\n");
	 fprintf(stderr, "usage: %s -f GENOTYPE_MATRIX_FILENAME -r NO_OF_RUNS\n", argv[0]);
	 fprintf(stderr, "Aborting...\n");
	 my_abort();
	 exit(1);
  }

  pmytime time= MYTIME_create();
  MYTIME_start(time);
  for (int i= 0; i<nrun; ++i) {
	 run(Xgm);
  }
  MYTIME_print_stop(time);
  MYTIME_destroy(time);
  BM_destroy(Xgm);
}
