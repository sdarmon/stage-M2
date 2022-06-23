import sys

Arg = sys.argv[:]
bubble = set()
with open(Arg[1], 'r') as f:
    for line in f:
        L=line.split('\t')
        bubble.add(L[:-4])
with open(Arg[2], 'r') as f:
    for line in f:
        if line[:-1] in bubble:
            bubble.remove(line[:-1])
with open(Arg[3], 'r') as f:
    for line in f:
        if line[:-1] in bubble:
            bubble.remove(line[:-1])
with open(Arg[4], 'r') as f:
    for line in f:
        if line[:-1] in bubble:
            bubble.remove(line[:-1])
for el in bubble:
    print(el)