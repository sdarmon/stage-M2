
import sys

Arg = sys.argv[:]

if len(Arg) not in [4]:
    print("Use : " + Arg[0] + " allNodes.txt Te.fa TeNodes_corrected.txt")
    exit()


#Function that returns the reverse complement of a sequence
def rev_comp(seq):
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    return ''.join([comp[el] for el in seq[::-1]])


#Reading the TE.fa file and storing the 41-mers of the TEs
Te = {}
name = ''
with open(Arg[2], 'r') as f:
    for line in f:
        if line[0] == '>':
            name = line[1:-1]
            Te[name] = ''
        else:
            Te[name] += line.strip()

#Defining TEnodes as the set of 41-mers of the TEs
TeNodes = set()
for TE in Te.values():
    for i in range(len(TE) - 41):
        TeNodes.add(TE[i:i+41])
        TeNodes.add(rev_comp(TE[i:i+41]))

#Reading the allNodes.txt file and writing the corrected file
with open(Arg[1], 'r') as f:
    with open(Arg[3], 'w') as f_out:
        for line in f:
            L = line.split('\t')
            kmers = [L[1][i:i+41] for i in range(0, len(L[1])-41)]
            #Checking if the kmer is in the TEnodes
            for kmer in kmers:
                if kmer in TeNodes:
                    f_out.write(line)
                    break