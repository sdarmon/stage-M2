# Cette fonction permet de convertir les unitigs ayant un poids 
# supérieur au threshold au format fastq. Note à moi-même, cela 
# sera utilisé par STAR qui prend également en entrée les fichiers
# au format fasta qui sont donc bien moins gros et compliqué à faire!!!
# L'option "-reserse" permet de récupérer les unitigs qui sont issus du
# fichier d'intersectionKissNoDouble.txt



import numpy as np 
import sys


Arg = sys.argv[:]


def isPoly(seq):
    """
    Entrée : une séquence de nucléotides
    Sortie : s'il existe une queue poly
    """
    n = len(seq)
    m = 0
    compt = 1
    last = seq[0]
    for i in range(1,n):
        if seq[i] == last:
            compt+=1
        else:
            m = max(m,compt)
            compt = 1
            last = seq[i]
    m = max(m,compt)
    if (m > 14):
        return(True)
    return(False)

if len(Arg) not in [4,5,6]:
    print("Use : "+Arg[0]+ " input.txt output.fq threshold [-reverse ref.txt] [-clean] \n Ne pas utiliser -reverse et -clean en même temps!")
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
            if isPoly(seq):
                continue
            f.write(">SEQ_"+str(compt)+"\n")
            compt+=1
            f.write(seq+"\n")
elif len(Arg) == 5: #i.e. argument -clean
    seqs = []
    t = int(Arg[3])
    with open(Arg[1],'r') as f:
        for line in f:
            if (len(line)<2):
                break
            L = line.split('\t')
            if (int(L[2][:-1])> t):
                seqs.append(line)
    with open(Arg[2],'w') as f:
        compt = 0
        for seq in seqs:
            if isPoly(seq):
                continue
            f.write(seq)
else:
    seqs = []
    ref = set()
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
            ref.add(int(L[12].split('_')[1]))
    with open(Arg[2],'w') as f:
        compt = 0
        for seq in seqs:
            if isPoly(seq.split('\t')[1]):
                continue
            if compt in ref:
                f.write(seq)
            compt+=1