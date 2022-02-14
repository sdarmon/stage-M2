import numpy as np 
import sys


Arg = sys.argv[:]

if len(Arg) not in [4,6]:
    print("Use : "+Arg[0]+ " input.txt output.fq threshold -reverse ref.txt")
    exit()
if len(Arg) == 4:
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
else:
    seqs = []
    ref = []
    t = int(Arg[3])
    with open(Arg[1],'r') as f:
        for line in f:
            if (len(line)<2):
                break
            L = line.split('\t')
            if (int(L[2][:-1])> t):
                seqs.append(line)
    with open(Arg[5],'r') as f2:
        for line in f2:
            if (len(line)<2):
                break
            L = line.split('\t')
            ref.append(int(L[12].split('_')[1]))
    with open(Arg[2],'w') as f:
        compt = 0
        for seq in seqs:
            if len(seq.split('\t')[1]) > 600:
                if compt in ref or (compt+1) in ref:
                   f.write(seq)
                compt+=2
            elif compt in ref:
                f.write(seq)
                compt+=1