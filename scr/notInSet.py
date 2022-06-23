import sys

Arg = sys.argv[:]
bubble = set()
correspondance = dict()
with open(Arg[1], 'r') as f:
    for line in f:
        L=line.split('\t')
        bubble.add(L[-4])
        correspondance[bubble]=line[:-1]
with open(Arg[2], 'r') as f:
    for line in f:
        b=line.split('\t')[-1][:-1]
        if b in bubble:
            bubble.remove(b)
with open(Arg[3], 'r') as f:
    for line in f:
        b=line.split('\t')[-1][:-1]
        if b in bubble:
            bubble.remove(b)
with open(Arg[4], 'r') as f:
    for line in f:
        b=line.split('\t')[-1][:-1]
        if b in bubble:
            bubble.remove(b)
for el in bubble:
    print(correspondance[el])