import numpy as np 
import sys


Arg = sys.argv[:]

if len(Arg) != 3:
    print("Use : "+Arg[0]+ " input.txt output.txt")
elif Arg[1] == "intersectionKiss.txt":
    with open(Arg[1],'r') as f:
        with open(Arg[2],'w') as o:
            oldTarget = ""
            for line in f:
                target = line.split("\t")[8]
                print(target)
                if oldline == "":
                    oldTarget = target
                    o.write(oldline)
                elif target == oldTarget:
                    continue
                else:
                    o.write(line)
                    oldTarget = target
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