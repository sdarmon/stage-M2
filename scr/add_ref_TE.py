#This code read a a comp%i.txt file and a seq_intersectionTE%i.txt file
#and edit the comp%i.txt file to add the TE reference in the seq_intersectionTE%i.txt file for every line
# that is present in the seq_intersectionTE%i.txt file


import sys

Arg = sys.argv[:]

if len(Arg) not in [6]:
    print("Use : " + Arg[0] + " comp%i.txt seq_intersectionTE%i.txt seq_intersectionRef%i.txt seq_intersectionConsensus%i.txt file.abundance")
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
        elif TE not in dic_index2TE[index]:
            dic_index2TE[index].append(TE)

#Reading intersectionRef%i.txt into a list to get the gene intersection with the comps
dic_index2gene = {}

with open(Arg[3], 'r') as f:
    for line in f:
        L = line.split('\t')
        gene = L[8].split('"')[1]
        index = int(L[12].split('_')[1])
        if index not in dic_index2gene:
            dic_index2gene[index] = [gene]
        elif gene not in dic_index2gene[index]:
            dic_index2gene[index].append(gene)


#Reading intersectionConsensus%i.txt into a list to get the TE consensus intersection with the comps
dic_index2cons = {}

with open(Arg[4], 'r') as f:
    for line in f:
        L = line.split('\t')
        gene = L[2]
        index = int(L[0].split('_')[1])
        if index not in dic_index2cons:
            dic_index2cons[index] = [gene]
        elif gene not in dic_index2cons[index]:
            dic_index2cons[index].append(gene)


#Reading the abundance file
abundance = []
with open(Arg[5], 'r') as f:
    for line in f:
        abundance.append(line.split('.')[0])


#Reading every line of comp%i.txt and adding the TE reference if the index of the line is in dic_index2TE
#Writing the new line in a new file comp%i.txt
i = 0
all_representants = []
with open(Arg[1], 'r') as f:
    with open(Arg[1].split('.')[0] + "_annoted.nodes", 'w') as f_out:
        for line in f:
            L = line.split('\t')
            node = int(L[0])
            representant = "*"
            representant_index = -1
            if i in dic_index2gene:
                str_gene = "; ".join(dic_index2gene[i])
                representant = dic_index2gene[i][0]
            else:
                str_gene = "*"
            if i in dic_index2TE:
                str_TE = "; ".join(dic_index2TE[i])
                representant = dic_index2TE[i][0].split("$")[0] 
            else :
                str_TE = "*"
            if representant not in all_representants:
                all_representants.append(representant)
            elif representant != "*":
                representant_index = all_representants.index(representant)

            if i  in dic_index2cons:
                str_cons = "; ".join(dic_index2cons[i])
            else:
                str_cons = "*"
            L=line[:-1].split("\t")
            ind=L[0]
            seq=L[1]
            if len(L) == 3:
                weight= L[2]
                d = "0"
            else:
                d = L[2]
                weight = L[3]
            f_out.write(ind +"\t" + seq + "\t" + d + "\t" + weight + "\t" + str_TE + "\t" + str_gene +  "\t" + str(representant_index) + "\t" + str_cons + "\t" +  abundance[int(L[0])] + "\n")

            i += 1
