import numpy as np 
import sys


Arg = sys.argv[1:]

if len(Arg) == 0:
    print("No argument given....")
else:
    y = []
    m = 0
    with open(Arg[0]) as f:
        for line in f:
            if (len(line)<2):
                break
            L = line.split('\t')
            y.append(int(L[1][:-1]))
            m = max(m,y[-1])
        x = np.arange(1,m+1,1)
        Y = [0 for i in range(m)]
        for el in y:
            Y[el-1] +=100/len(y)
        for el in Y:
            print(el)
