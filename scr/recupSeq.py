#


import sys

Arg = sys.argv[:]

if len(Arg) not in [5]:
    print("Use : " + Arg[0] + " prefixCompTe prefixComp nbComponent pos") #Pos = 12 pour le chien
    exit()
if len(Arg) == 5:
    for i in range(int(Arg[3])):
        seq = set()
        with open(Arg[1]+str(i)+".txt",'r') as f:
            for line in f:
                l = int(line.split("\t")[int(Arg[4])][4:])
                seq.add(l)
        with open(Arg[2]+str(i)+".txt",'r') as f:
            index = 0
            for line in f:
                if line[0]==">":
                    continue
                if index in seq:
                    print(line[:-1])
                index+=1
