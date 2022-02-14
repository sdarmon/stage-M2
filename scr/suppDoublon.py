import numpy as np 
import sys


Arg = sys.argv[:]

if len(Arg) not in [3,5]:
    print("Use : "+Arg[0]+ " input.txt output.txt -s pos")
elif len(Arg) == 5 and Arg[3] == "-s":
    with open(Arg[1],'r') as f:
        with open(Arg[2],'w') as o:
            Vu = []
            for line in f:
                if line[:2] == "##":
                    continue
                target = int((line.split("\t")[int(Arg[4])]).split("_")[1])
                if target >= len(Vu):
                    for i in range(len(Vu),target+10):
                        Vu.append(False)
                if Vu[target]:
                    continue
                else:
                    o.write(line)
                    Vu[target] = True
elif len(Arg) == 5 and Arg[3] == "-t":
    with open(Arg[1],'r') as f:
        with open(Arg[2],'w') as o:
            Vu = set()
            for line in f:
                if line[:2] == "##":
                    continue
                target = line.split("\t")[int(Arg[4])]
                if target in Vu:
                    continue
                else:
                    o.write(line)
                    Vu.add(target[:])
else:
    with open(Arg[1],'r') as f:
        with open(Arg[2],'w') as o:
            oldline = ""
            for line in f:
                if oldline == "":
                    oldline = line
                    o.write(oldline)
                elif line == oldline:
                    continue
                else:
                    o.write(line)
                    oldline=line