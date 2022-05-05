#Ce programme permet de savoir quelles bulles (d'un fichier bulle.fa) traversent les composantes (du dossier dirComp)
#On indique si les chemins traversent les composantes par le chemin du haut, du bas ou des deux.

import sys

Arg = sys.argv[:]

#Fonction renvoyant le reverse complement d'une séquence de nucléotides
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

#Renvoie une séquence de nucléotides en masjucule
def maj(sequ):
    s = ""
    for el in sequ:
        if el == 'a':
            s = s + 'A'
        elif el == 'c':
            s = s + 'C'
        elif el == 'g':
            s = s + 'G'
        elif el == 't':
            s = s + 'T'
        else:
            s = s + el
    return s


if len(Arg) not in [5,6,7,8]:
    print("Use : " + Arg[0] + " dirComp bulle.fa nbComp threshold -k kmer -rapport\n\t k = 41 by default.")
    exit()
if len(Arg) in [5,6,7,8]:
    titre = "" #Variable contenant le titre de la séquence dans le fichier .fa
    seq = "" #Sequence de nucléotides
    if (len(Arg)>= 7 and Arg[5][1] == 'k') :
        k = int(Arg[6])
    else:
        k = 41
    threshold = int(Arg[4])
    strat = [i*(threshold+1)//6 for i in range(1,6)]
    kmerFrom = [dict() for i in range(1,6)] #Dictionnaire associant à un kmer donné la composante et strat dont il est issu
    #Récupération des séquences des composantes
    for i in range(int(Arg[3])):
        with open(Arg[1] + "comp" + str(i) + ".txt", 'r') as f:
            for line in f:
                if len(line) < 2:
                    continue
                L = line.split("\t")
                #On l'ajoute au dictionnaire ainsi que son complémentaire
                compl= reverseC(L[1])
                for pos in range(len(L[1])-k+1):
                    st = min((int(L[2][:-1]) // ((threshold+1)//6)) -1,4)
                    kmerFrom[st][L[1][pos:k+pos]] = i
                    kmerFrom[st][compl[pos:k+pos]] = i

    #On parcourt maintenant les bulles
    with open(Arg[2], 'r') as f:
        upper = False #Variable booléenne indiquant si on est un upper path ou un under path
        for line in f:
            if len(line) < 2:
                continue
            if titre == "": #Cas où la ligne est un titre
                titre = line[:-1]
                upper = not upper #On change la valeur de upper
                if upper:
                    titreUpper = titre[:-1] #On stocke la titre
                else:
                    titreUnder = titre[:-1]
                continue
            seq = line[:-1]
            seq = maj(seq) #Les séquences ayant des nucléotides écrits en miniscule, on les passe en majuscule

            #On récupère les numéros des composantes contenant des kmers de la séquence
            comp_possible = [set() for i in range(1,6)]
            trouve = 0
            for pos in range(len(seq) - k + 1):
                mer = seq[pos: pos + k]
                for st in range(5):
                    ind = kmerFrom[st].get(mer, -1)
                    if ind != -1:
                        comp_possible[st].add(ind)
                        trouve = 1
            if upper: #Cas chemin du haut
                seqUpper = line[:-1]
                upperComp = [el for el in comp_possible]
                trouveUpper = trouve
            else: #Cas chemin du bas
                seqUnder = line[:-1]
                underComp = [el for el in comp_possible]
                trouveUnder = trouve
                #On peut écrit la bulle et son rapport si c'est intéressant
                if trouveUnder or trouveUpper:
                    if Arg[-1] == "-rapport":
                        print(titreUpper)
                        print(seqUpper)
                        print(titreUnder)
                        print(seqUnder)
                    text = ""
                    printing = False
                    for st in range(1,5):
                        if len(upperComp[st]) != 0 or len(underComp[st]) != 0:
                            text += "In strat ["+str((st+1)*(threshold+1)//6) +"," + str(-1+(st+2)*(threshold+1)//6) + "] :\t"
                            A = upperComp[st] & underComp[st]
                            B = upperComp[st] - underComp[st]
                            C = underComp[st] - upperComp[st]
                            if len(A) != 0:
                                t = ""
                                printing = True
                                for el in A:
                                    t+= " " + str(el)
                                text+= "in both path :" + t + "\t"
                                if len(B) != 0 or len(C) != 0:
                                    printing = False
                            if False :#len(B) != 0:
                                t = ""
                                printing = True
                                for el in B:
                                    t+= " " + str(el)
                                text+= "only in upper :" + t + "\t"
                            if False : #len(C) != 0:
                                t = ""
                                printing = True
                                for el in C:
                                    t+= " " + str(el)
                                text+= "only in under :" + t + "\t"
                            if printing:
                                print(text)
            titre = ""#On part pour la ligne suivante qui sera un titre
