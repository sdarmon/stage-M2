import sys

Arg = sys.argv[:]
seq = []
with open(Arg[1], 'r') as f:
    with open(Arg[2], 'w') as o:
        for line in f:
            L=line.split("\t")
            o.write("\t".join(L[:8]+["".join(L[8:])]))