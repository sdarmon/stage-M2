#

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


if len(Arg) not in [4]:
    print("Use : " + Arg[0] + " dirComp bulle.fa nbComp")
    exit()
if len(Arg) == 4:
    seq = ""
    titre = ""
    comp = [[] for i in range(int(Arg[3]))]
    kmer = dict()
    unitig = dict()
    k = 41
    for i in range(int(Arg[3])):
        with open(Arg[1] + "comp" + str(i) + ".txt", 'r') as f:
            if len(line) < 2:
                continue
            for line in f:
                L = line.split("\t")
                kmer[L[1][:k]] = int(L[0])
                kmer[reverseC([1])[:k]] = int(L[0])
                unitig[int(L[0])] = L[1]
                comp[i].append(int(L[0]))

    with open(Arg[2], 'r') as f:
        upper = 0
        for line in f:
            if len(line) < 2:
                continue
            if titre == "":
                titre = line[:-1]
                upper = (upper + 1) % 2
                if upper:
                    titreUpper = titre[:-1]
                else:
                    titreUnder = titre[:-1]
                continue
            seq = line[:-1]
            seq = maj(seq)
            comp_possible = [0 for x in range(int(Arg[3]))]
            for pos in range(len(seq) - k + 1):
                mer = seq[pos, pos + k]
                L = kmer.get(mer, [])
                for comp in L:
                    comp_possible[comp] += 1
            if upper:
                seqUpper = line[:-1]
                upperComp = set()
                for i in range(len(comp_possible)):
                    if comp_possible[i] != 0:
                        upperComp.add(i)
            else:
                seqUnder = line[:-1]
                underComp = set()
                for i in range(len(comp_possible)):
                    if comp_possible[i] != 0:
                        underComp.add(i)
                if len(underComp) != 0 or len(upperComp) != 0:
                    print(titreUpper)
                    print(seqUpper)
                    print(titreUnder)
                    print(seqUnder)
                    print("Both found in components:", upperComp & underComp, "\t Only in upper:",
                          upperComp - underComp, "\t Only in under:", underComp - upperComp)

            titre = ""
