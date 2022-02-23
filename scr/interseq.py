# Petit script permettant faire les intersections entre les
# lignes de deux fichiers

import numpy as np 
import sys


Arg = sys.argv[:]

if len(Arg) not in [4]:
    print("Use : "+Arg[0]+ " seq1.txt seq2.txt outputPrefix ")
else:
    inter = set()
    seq2Not = set()
    with open(Arg[1],'r') as seq1:
        Sq = set(seq1.readlines())
    with open(Arg[2],'r') as seq2:
        for line in seq2:
            if line in Sq:
                inter.add(line)
                Sq.remove(line)
            else:
                seq2Not.add(line)

    with open(Arg[3]+"intersection.txt",'w') as out:
        for line in inter:
            out.write(line)

    with open(Arg[3]+"seq1Missing.txt",'w') as out:
        for line in Sq:
            out.write(line)

    with open(Arg[3]+"seq2Missing.txt",'w') as out:
        for line in seq2Not:
            out.write(line)