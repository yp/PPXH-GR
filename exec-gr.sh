#!/bin/bash


rm -f RESULT*

echo "" | LD_LIBRARY_PATH=. ./grealLinux -i hypergraph.txt > /dev/null

if fgrep "NO REALIZATION" RESULT_realization.txt > /dev/null; then
exit 1;
else
sed 's/;/\n/g' < RESULT_neato_format.txt | grep '\-\-' | sed 's/n\([0-9]\+\)--n\([0-9]\+\) .label=.\([0-9,]\+\).\+/\1 \2 \3/' > realization.txt
exit 0;
fi;
