import numpy as np 
import sys


Arg = sys.argv[:]

if len(Arg) not in [4]:
    print("Use : "+Arg[0]+ " inputRead.txt selectedReads.txt output.txt")
    exit()
if len(Arg) == 4:
    with open(Arg[1],'r') as f:
        with open(Arg[2],'r') as read:
            with open(Arg[3],'w') as o:
                reads = read.readlines()
                index = 0
                for line in f:
                    if len(line) <2 or len(reads) == index:
                        break
                    pos = int(line.split('\t')[0])
                    if pos == int(reads[index].split('\t')[0]):
                        o.write(line[:-1]+'\t1\n')
                        index+=1
                    else:
                        o.write(line[:-1]+'\t0\n')