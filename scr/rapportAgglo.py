# Cette fonction permet d'obtenir le rapport post agglomeration'
# à partir des intersections avec les TE de chacune des composantes
# En particulier, cela affiche un histogramme récapitulatif

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys

Arg = sys.argv[:]

if len(Arg) not in [4]:
    print("Use : " + Arg[0] + " TE.gff prefixCompTe nbComponent")
    exit()
if len(Arg) == 4:

    dicTE = {}
    with open(Arg[1], 'r') as f:
        for line in f:
            dicTE[line[:-1]] = []

    for i in range(int(Arg[3])):
        with open(Arg[2] + str(i) + ".txt", 'r') as f:
            for line in f:
                dicTE[line[:-1]].append(i)

    X = [i for i in range(int(Arg[3]))]
    Y = [[0 for i in range(int(Arg[3]))] ]
    for lst in dicTE.values():
        if not lst:
            continue
        freq = len(lst)
        for comp in lst:
            if len(Y) <= freq:
                for i in range(freq - len(Y) + 1):
                    Y.append([0 for i in range(int(Arg[3]))])
            Y[freq][comp] += 1

    for i in range(1,len(Y)):
        plt.bar(X,Y[i],label = str(i), color = cm.hsv(i/len(Y)))
    plt.title("Nombre de TE dans chaque composante")
    plt.ylabel("Nombre de TE distincts")
    plt.xlabel("Numéro de composante")
    plt.legend(title = "Nombre de matchs dans des composantes distinctes")
    plt.show()