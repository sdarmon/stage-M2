import numpy as np 
import sys


Arg = sys.argv[:]

if len(Arg) != 4:
    print("Use : "+Arg[0]+ " input.txt output.fq threshold")
else:
    seqs = []
    t = Arg[3]
    with open(Arg[1],'r') as f:
        for line in f:
            if (len(line)<2):
                break
            L = line.split('\t')
            if (int(L[2][:-1])> t):
                seqs.append(int(L[1]))
    with open(Arg[2],'w') as f:
        compt = 0
        for seq in seqs:
            f.write("@SEQ_"+str(compt))
            compt+=1
            f.write(seq)
            f.write("+")
            f.write("J"*len(seq))
