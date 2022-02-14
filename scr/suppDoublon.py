import numpy as np 
import sys


Arg = sys.argv[:]

if len(Arg) not in [3,5]:
    print("Use : "+Arg[0]+ " input.txt output.txt -t pos")
elif len(Arg) == 5:
    with open(Arg[1],'r') as f:
        with open(Arg[2],'w') as o:
            oldTarget = ""
            for line in f:
                target = line.split("\t")[int(Arg[4])]
                if oldTarget == "":
                    oldTarget = target
                    o.write(line)
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