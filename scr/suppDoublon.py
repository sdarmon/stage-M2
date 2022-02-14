import numpy as np 
import sys


Arg = sys.argv[:]

if len(Arg) != 3:
    print("Use : "+Arg[0]+ " input.txt output.txt")
else:
    with open(Arg[1],'r') as f:
        with open(Arg[2],'w') as o:
            oldline = ""
            for line in f:
                if oldline == "":
                    oldline = line[:-1]
                    o.write(oldline)
                elif line[:-1] == oldline:
                    continue
                else:
                    o.write(line)
                    oldline=line[:-1]