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

    dicTE = {} #On crée un dictionnaire qui va associer à des TE les numéros des composantes dont au moins l'un de ses
    #sommets s'intersecte avec le-dit TE.
    with open(Arg[1], 'r') as f: #On ouvre le fichier de référence et on initialise le dicTE avec les TE
        for line in f:
            if len(line) > 2 and line[0] != '#':
                if Arg[4] == "-target":
                    target = line.split("\t")[8].split(" ")[0]
                else:
                    target = line.split("\t")[8].split(";")[0]

                if '(' in target and ')n' in target:
                    continue
                dicTE[target] = []

    for i in range(int(Arg[3])): #On boucle sur les composantes
        if int(Arg[3]) == 1:  # Cas où l'on ne veut analyser qu'un seul fichier
            file = Arg[2]
        else:
            file = Arg[2] + str(i) + ".txt"  # Cas on veut bien analyser plusieurs fichiers
        with open(file, 'r') as f: #On ouvre l'intersection de la i-ème composante
            for line in f:
                if len(line) > 2: #on récupère le nom du TE
                    if Arg[4] == "-target":
                        target = line.split("\t")[8].split(" ")[0]
                    else:
                        target = line.split("\t")[8].split(";")[0]
                    if '(' in target and ')n' in target:
                        continue
                    if dicTE[target] == [] or dicTE[target][-1] != i: #Si on n'a pas déjà trouvé d'intersection avec
                        #cette composante, on l'ajoute au dicTE. Remarque: ici on pourrait break et enlever la seconde
                        #condition.
                        dicTE[target].append(i)

    #Variable pour plt.plot le nombre d'éléments afin obtenu
    X = [i for i in range(int(Arg[3]))]
    Y = [[0 for i in range(int(Arg[3]))]]
    for target, lst in dicTE.items(): #On boucle sur le dictionnaire
        if not lst: #Cas de TE n'ayant pas été retrouvé dans les composantes
            continue
        freq = len(lst) #On définit la fréquence comme étant le nombre de composantes dans lesquels ce TE a été intersecté
        for comp in lst:
            if len(Y) <= freq: #Si jamais Y n'est pas assez grand pour rajouter que cette fréquence comme vu dans les
                #composantes de lst
                for i in range(freq - len(Y) + 1):
                    Y.append([0 for j in range(int(Arg[3]))])
            Y[freq][comp] += 1

    #On inverse les dimensions de Y dans Z
    Z = [[] for el in Y]
    for target, lst in dicTE.items():
        Z[len(lst)].append(target)
    for i in range(1, len(Z)): #Afin d'afficher les stastistiques des composantes
        print("Found in " + str(i) + " componants :", Z[i])

    #On plt.plot les histogrammes des fréquences ainsi trouvée (avec des jolies couleurs et les bars se superposant
    for i in range(1, len(Y)):
        plt.bar(X, Y[i], label=str(i), color=cm.hsv(i / len(Y)), bottom=np.sum(Y[0:i], axis=0))
    plt.title("Nombre de TE dans chaque composante")
    plt.ylabel("Nombre de TE distincts")
    plt.xlabel("Numéro de composante")
    plt.legend(title="TE contenus dans X composante(s)")
    plt.show()
