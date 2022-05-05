import sys

Arg = sys.argv[:]

with open(Arg[1],'r') as f:
    for line in f:
        print(line)
        print(line[:-1])