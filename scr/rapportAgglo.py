# Cette fonction permet d'obtenir le rapport post agglomeration'
# à partir des intersections avec les TE de chacune des composantes
# En particulier, cela affiche un histogramme récapitulatif

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import sys

Arg = sys.argv[:]

if len(Arg) not in [5]:
    print("Use : " + Arg[0] + " TE.gff prefixCompTe nbComponent [-chien/-target]")
    exit()
if len(Arg) == 5:

    dicTE = {}
    with open(Arg[1], 'r') as f:
        for line in f:
            if len(line) > 2 and line[0] != '#':
                if Arg[4] == "-target":
                    target = line.split("\t")[8].split(" ")[0]
                else:
                    target = line.split("\t")[8].split(";")[0]

                if '(' in target and ')n' in target:
                    continue
                dicTE[target] = []

    for i in range(int(Arg[3])):
        if int(Arg[3]) == 1:  # Cas où l'on ne veut analyser qu'un seul fichier
            file = Arg[2]
        else:
            file = Arg[2] + str(i) + ".txt"  # Cas on veut bien analyser plusieurs fichiers
        with open(file, 'r') as f:
            for line in f:
                if len(line) > 2:
                    if Arg[4] == "-target":
                        target = line.split("\t")[8].split(" ")[0]
                    else:
                        target = line.split("\t")[8].split(";")[0]
                    if '(' in target and ')n' in target:
                        continue
                    if dicTE[target] == [] or dicTE[target][-1] != i:
                        dicTE[target].append(i)

    X = [i for i in range(int(Arg[3]))]
    Y = [[0 for i in range(int(Arg[3]))]]
    for target, lst in dicTE.items():
        if not lst:
            continue
        freq = len(lst)
        for comp in lst:
            if len(Y) <= freq:
                for i in range(freq - len(Y) + 1):
                    Y.append([0 for j in range(int(Arg[3]))])
            Y[freq][comp] += 1
    Z = [[] for el in Y]
    for target, lst in dicTE.items():
        Z[len(lst)].append(target)
    for i in range(1, len(Z)):
        print("Found in " + str(i) + " componants :", Z[i])
    for i in range(1, len(Y)):
        plt.bar(X, Y[i], label=str(i), color=cm.hsv(i / len(Y)), bottom=np.sum(Y[0:i], axis=0))
    plt.title("Nombre de TE dans chaque composante")
    plt.ylabel("Nombre de TE distincts")
    plt.xlabel("Numéro de composante")
    plt.legend(title="TE contenus dans X composante(s)")
    plt.show()
