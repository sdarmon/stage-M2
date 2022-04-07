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


if len(Arg) not in [4,5]:
    print("Use : " + Arg[0] + " dirComp bulle.fa nbComp -rapport")
    exit()
if len(Arg) in [4,5]:
    titre = "" #Variable contenant le titre de la séquence dans le fichier .fa
    seq = "" #Sequence de nucléotides
    kmerFrom = dict() #Dictionnaire associant à un kmer donné la composante dont il est issu
    k = 41 #k du kmer, il serait judicieux de le mettre en argument de la fonction

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
                    kmerFrom[L[1][pos:k+pos]] = i
                    kmerFrom[compl[pos:k+pos]] = i

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
            comp_possible = set()
            for pos in range(len(seq) - k + 1):
                mer = seq[pos: pos + k]
                ind = kmerFrom.get(mer, -1)
                if ind != -1:
                    comp_possible.add(ind)
            if upper: #Cas chemin du haut
                seqUpper = line[:-1]
                upperComp = comp_possible
            else: #Cas chemin du bas
                seqUnder = line[:-1]
                underComp = comp_possible
                #On peut écrit la bulle et son rapport si c'est intéressant
                if len(underComp) != 0 or len(upperComp) != 0:
                    if len(Arg)==4:
                        print(titreUpper)
                        print(seqUpper)
                        print(titreUnder)
                        print(seqUnder)
                    print("Both found in components:", upperComp & underComp, "\t Only in upper:",
                          upperComp - underComp, "\t Only in under:", underComp - upperComp)
            titre = ""#On part pour la ligne suivante qui sera un titre
