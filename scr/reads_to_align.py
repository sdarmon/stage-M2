import numpy as np 
import sys


Arg = sys.argv[:]

if len(Arg) != 4:
    print("Use : "+Arg[0]+ " input.txt output.fq threshold")
else:
    seqs = []
    t = int(Arg[3])
    with open(Arg[1],'r') as f:
        for line in f:
            if (len(line)<2):
                break
            L = line.split('\t')
            if (int(L[2][:-1])> t):
                seqs.append(L[1])
    with open(Arg[2],'w') as f:
        compt = 0
        for seq in seqs:
            if len(seq)> 600:
                s1 = seq[:len(seq)//2]
                s2 = seq[len(seq)//2:]
                f.write("@SEQ_"+str(compt)+"\n")
                compt+=1
                f.write(s1+"\n")
                f.write("+"+"\n")
                f.write("J"*len(s1)+"\n")
                f.write("@SEQ_"+str(compt)+"\n")
                compt+=1
                f.write(s2+"\n")
                f.write("+"+"\n")
                f.write("J"*len(s2)+"\n")
            else:
                f.write("@SEQ_"+str(compt)+"\n")
                compt+=1
                f.write(seq+"\n")
                f.write("+"+"\n")
                f.write("J"*len(seq)+"\n")
