# Cette fonction permet de convertir les unitigs ayant un poids 
# supérieur au threshold au format fastq. Note à moi-même, cela 
# sera utilisé par STAR qui prend également en entrée les fichiers
# au format fasta qui sont donc bien moins gros et compliqué à faire!!!
# L'option "-reserse" permet de récupérer les unitigs qui sont issus du
# fichier d'intersectionKissNoDouble.txt



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
            if len(seq)> 1200:
                s1 = seq[:len(seq)//4]
                s2 = seq[len(seq)//4:2*len(seq)//4]
                s3 = seq[2*len(seq)//4:3*len(seq)//4]
                s4 = seq[3*len(seq)//4:]
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
                f.write("@SEQ_"+str(compt)+"\n")
                compt+=1
                f.write(s3+"\n")
                f.write("+"+"\n")
                f.write("J"*len(s3)+"\n")
                f.write("@SEQ_"+str(compt)+"\n")
                compt+=1
                f.write(s4+"\n")
                f.write("+"+"\n")
                f.write("J"*len(s4)+"\n")
            
            if len(seq)> 900:
                s1 = seq[:len(seq)//3]
                s2 = seq[len(seq)//3:2*len(seq)//3]
                s3 = seq[2*len(seq)//3:]
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
                f.write("@SEQ_"+str(compt)+"\n")
                compt+=1
                f.write(s3+"\n")
                f.write("+"+"\n")
                f.write("J"*len(s3)+"\n")
            elif len(seq)> 600:
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
            if len(seq.split('\t')[1]) > 1200:
                if compt in ref or (compt+1) in ref or (compt+2) in ref or (compt+3) in ref:
                   f.write(seq)
                compt+=4
            elif len(seq.split('\t')[1]) > 900:
                if compt in ref or (compt+1) in ref or (compt+2) in ref:
                   f.write(seq)
                compt+=3
            elif len(seq.split('\t')[1]) > 600:
                if compt in ref or (compt+1) in ref:
                   f.write(seq)
                compt+=2
            else:
                if compt in ref:
                    f.write(seq)
                compt+=1