#

import sys

Arg = sys.argv[:]

if len(Arg) not in [3]:
    print("Use : " + Arg[0] + "dirComp bulle.fa nbComp")
    exit()
if len(Arg) == 3:
    seq = ""
    titre = ""
    comp = [[] for i in range(Arg[3])]

    for i in range(Arg[3]):
        with open(Arg[1] + "comp" + str(i) + ".txt", 'r') as f:
            if len(line) < 2:
                continue
            for line in f:
                comp[i].append(line)

    with open(Arg[2], 'r') as f:
        for line in f:
            if len(line) < 2:
                continue
            if titre == "":
                titre = line
                continue
            seq = line
            for i in range(Arg[3]):
                if not comp[i]:
                    continue
                with open(Arg[2] + "interbulle" + str(i) + ".txt", 'a') as o:
                    for el in comp[i]:
                        if el.split("\t")[1] in seq:
                            o.write(titre)
                            o.write(el)
            titre = ""
