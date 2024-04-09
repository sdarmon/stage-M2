#This code read a a comp%i.txt file and a seq_intersectionTE%i.txt file
#and edit the comp%i.txt file to add the TE reference in the seq_intersectionTE%i.txt file for every line
# that is present in the seq_intersectionTE%i.txt file


import sys

Arg = sys.argv[:]

if len(Arg) not in [3]:
    print("Use : " + Arg[0] + " comp%i.txt seq_intersectionTE%i.txt")
    exit()

#Reading seq_intersectionTE%i.txt into a list
dic_index2TE = {}

with open(Arg[2], 'r') as f:
    for line in f:
        L = line.split('\t')
        TE = L[8].split('"')[3]
        index = int(L[12].split('_')[1])
        if index not in dic_index2TE:
            dic_index2TE[index] = [TE]
        else:
            dic_index2TE[index].append(TE)

#Reading every line of comp%i.txt and adding the TE reference if the index of the line is in dic_index2TE
#Writing the new line in a new file comp%i.txt
i = 0
with open(Arg[1], 'r') as f:
    with open(Arg[1][:-4] + "_TE.txt", 'w') as f_out:
        for line in f:
            if i in dic_index2TE:
                f_out.write(line[:-1] + "\t" + "; ".join(dic_index2TE[i]) + "\n")
            else:
                f_out.write(line[:-1] + "\t" + ";" + "\n")
            i += 1