

import numpy as np 
import sys
import matplotlib.pyplot as plt


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
            y.append(int(L[2][:-1]))
            m = max(m,y[-1])
        x = np.arange(1,m+1,1)
        Y = [0 for i in range(m)]
        for el in y:
            Y[el-1] +=100/len(y)
    xb = []
    yb = []

    xg = []
    yg = []

    for i in range(len(x)):
        if Y[i] == float(0):
            yb.append(Y[i])
            xb.append(x[i])
        else:
            yg.append(Y[i])
            xg.append(x[i])

    fig, axs = plt.subplots(2)

    axs[0].plot(xg,yg, 'bo', label = 'Positive Frequency')
    axs[0].legend()
    axs[1].plot(xb,yb, 'ro', label = 'Null Frequency')
    axs[1].legend()
    axs[0].plot(x,Y, 'green')
    axs[1].plot(x,Y, 'green')
    axs[0].set(ylabel= "Frequencies (in %)")
    axs[1].set(ylabel= "Frequencies (in %)")
    plt.xlabel("Sizes of weight")
    fig.suptitle("Frequencies of the different sizes of weight (Case radius = 200)")
    plt.show()
