import sys

Arg = sys.argv[:]
step = -1
with open(Arg[1], 'r') as f:
    for line in f:
        if line[0] in ['-','_','#']:
            continue
        step = (step+1)%4
        if step == 0:
            L = line.split("|")
            id = L[0][1:]+"|"+L[1]
            print(id)
