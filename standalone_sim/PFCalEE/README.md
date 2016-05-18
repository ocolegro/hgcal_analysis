Check out Anne-Marie's code on github.
Follow her well thought out instructions.
To generate the same events that I used in my studies, please run:

for i in `seq 0 1000`; do python submitProd.py -s 1nd -q 2nd -t hexaV02-01-01 -g -r ${i} -v 30 -m 0 -e /store/cmst3/group/hgcal/Geant4 -o /afs/cern.ch/work/o/ocolegro/electron_single -d e- -n 2500 -a 10000 ; done
