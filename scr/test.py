import sys

Arg = sys.argv[:]
def reverseC(seq):
    s = ""
    for el in seq:
        if el == 'A':
            s = 'T' + s
        elif el == 'T':
            s = 'A' + s
        elif el == 'C':
            s = 'G' + s
        elif el == 'G':
            s = 'C' + s
    return s

with open(Arg[1],'r') as f:
    for line in f:
        seq = line[:-1]
        print(reverseC(seq))