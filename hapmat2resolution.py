#!/usr/bin/env python

########################################################################
#                 PPXH-GR
# A heuristic algorithm for the Pure Parsimony Xor Haplotyping problem
#
# Copyright (C) 2008 Yuri Pirola <yuri.pirola(-at-)gmail.com>
#
# Distributed under the terms of the GNU General Public License (GPL)
#
#
# This file is part of PPXH-GR.
#
# PPXH-GR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PPXH-GR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PPXH-GR. If not, see <http://www.gnu.org/licenses/>.
#
########################################################################

########################################################################
#
#  hapmat2resolution.py
#
#  A program that reads the haplotype matrixes computed by PPXH-GR
#  (possibly in several runs) and outputs, for each xor-genotype, the
#  pair(s) of haplotypes that solves it.
#
########################################################################


import itertools
import json
import logging
import optparse
import os
import os.path
import sys


parser= optparse.OptionParser(usage="usage: "
                              "%prog -i <XOR-GENOTYPES FILE> -o <HAPLOTYPE MATRIXES FILE> [-v]")
parser.add_option("-i", "--input",
                  action="store", dest="xorgen",
                  type="string", default=None,
                  help="the file containing the input xor-genotype matrix.",
                  metavar="FILE")
parser.add_option("-o", "--output",
                  action="store", dest="hapmat",
                  type="string", default=None,
                  help="the file containing the results computed by PPXH-GR "
                  "(i.e. the haplotype matrixes)",
                  metavar="FILE")
parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose",
                  default=False,
                  help="print additional log messages")

(options, args)= parser.parse_args()



log_level= logging.DEBUG if options.verbose else logging.INFO

logging.basicConfig(level=log_level,
                    format='%(levelname)-6s [%(asctime)s]  %(message)s')


logging.info("PPXH-GR -- hapmat2resolution")

if (  options.xorgen is None
      or not os.path.isfile(options.xorgen)
      or not os.access(options.xorgen, os.R_OK)):
    logging.fatal("Input file for xor-genotypes not specified or invalid. Given: '%s'.",
                  options.xorgen)
    sys.exit("File '{}' not found.".format(options.xorgen))

if (  options.hapmat is None
      or not os.path.isfile(options.hapmat)
      or not os.access(options.hapmat, os.R_OK)):
    logging.fatal("Input file for haplotype matrices not specified or invalid. Given: '%s'.",
                  options.hapmat)
    sys.exit("File '{}' not found.".format(options.hapmat))


####
#
# Read xor-genotypes
#
####

# 0/1 -> homozygous
# 2   -> heterozygous
gen_enc= { "0": 0,
           "1": 0,
           "2": 1 }

logging.info("Reading xor-genotypes file '%s'...", options.xorgen)
xgens= []
with open(options.xorgen, "r") as fxg:
    for line in fxg:
        line= line.strip()
        xgens.append([ gen_enc[g] for g in line if g in ('0', '1', '2') ])

n_gen= len(xgens)
assert n_gen>0, "The xor-genotypes file does not contain valid xor-genotypes."
l_gen= len(xgens[0])
assert all(( len(xg) == l_gen for xg in xgens )), \
    "The xor-genotypes file contains xor-genotypes with different lengths."

logging.info("Read %d xor-genotypes of length %d.", n_gen, l_gen)

logging.info("Reading haplotype matrixes file '%s'...", options.hapmat)
hapmats= []
with open(options.hapmat, "r") as fhm:
    in_mat= False
    for line in fhm:
        line= line.strip()
        if line == "==Haplotype Matrix":
            assert not in_mat, "Malformed haplotype matrixes file."
            in_mat= True
            hapmats.append(set())
            logging.debug("Reading haplotype matrix #%d...", len(hapmats))
        elif line == "==END Haplotype Matrix":
            assert in_mat, "Malformed haplotype matrixes file."
            in_mat= False
        elif in_mat:
            assert len(line)==l_gen, \
                "The haplotype length is not equal to the xor-genotype length."
            assert all(( h in ('0', '1') for h in line)), \
                "The haplotype has invalid characters."
            hapmats[-1].add(tuple( int(h) for h in line ))
        else:
            # Not-interesting line
            pass
logging.info("Read %d haplotype matrixes.", len(hapmats))

min_hap= min((len(hm) for hm in hapmats))

logging.info("A most parsimonious matrix has %d haplotypes.", min_hap)

# Get the first most parsimonious matrix
hapmat= next(( hm for hm in hapmats if len(hm) == min_hap))
hapmat= [ h for h in hapmat ]

results= {}
results['xor-genotypes']= [ "".join((str(g) for g in xg)) for xg in xgens ]
results['haplotypes']= [ "".join((str(h) for h in hapl)) for hapl in hapmat ]

def xg_solve(xg, h1, h2):
    return all([ abs(h1i-h2i)==xgi
                 for xgi, h1i, h2i in zip(xg, h1, h2) ])

logging.info("Computing resolutions...")
resolutions=[]
for i,xg in itertools.izip(range(n_gen), xgens):
    logging.debug("Considering xor-genotype %d (of %d)...", i, n_gen)
    resxg={ 'xor-genotype-id': i,
            'xor-genotype-str': results['xor-genotypes'][i] }
    hapl_pairs= [ { 'haplotype1-id': h1id,
                    'haplotype1-str': results['haplotypes'][h1id],
                    'haplotype2-id': h2id,
                    'haplotype2-str': results['haplotypes'][h2id] }
                  for h1id, h1 in itertools.izip(range(min_hap), hapmat)
                  for h2id, h2 in itertools.izip(range(min_hap), hapmat)
                  if h1id <= h2id and xg_solve(xg, h1, h2) ]
    assert len(hapl_pairs)>0, \
        "Xor-genotype {} has no resolutions.".format(i)

    logging.debug("    ...xor-genotype %d has %d resolutions.", i, len(hapl_pairs))

    resxg['haplotype-pairs']= hapl_pairs


    resolutions.append(resxg)

results['resolutions']= resolutions


json.dump(results, sys.stdout, indent=4)
logging.info("PPXH-GR -- hapmat2resolution -- Completed")
