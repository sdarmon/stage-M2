import numpy as np 
import sys


Arg = sys.argv[:]

if len(Arg) not in [5]:
    print("Use : "+Arg[0]+ " edges.txt inputRead.txt threshold output.txt")
    exit()
if len(Arg) == 5:
    node = set()
    with open(Arg[1],'r') as f:
        for line in f:
            if len(line) <2:
                break
            L = line.split("\t")
            node.add(int(L[0]))
            node.add(int(L[1]))
    with open(Arg[2],'r') as f:
        with open(Arg[4],'w') as o:
            th = int(Arg[3])
            for line in f:
                if len(line) <2:
                    break
                pos = int(line.split('\t')[0])
                weight = int(line.split('\t')[2][:-1])
                if weight >= th:
                    if pos in node:
                        o.write(line[:-1]+'\t1\n')
                else:
                    if pos in node:
                        o.write(line[:-1]+'\t0\n')