
import sys

Arg = sys.argv[:]

if len(Arg) not in [4]:
    print("Use : " + Arg[0] + " allNodes.txt TeNodes.txt TeNodes_corrected.txt")
    exit()


#Function that returns the reverse complement of a sequence
def rev_comp(seq):
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    return ''.join([comp[el] for el in seq[::-1]])


#Reading the TEnodes
TeNodes = set()
with open(Arg[2], 'r') as f:
    for line in f:
        L = line.split('\t')
        #dividing L[1] into its every 41-mers
        for i in range(0, len(L[1])-41):
            TeNodes.add(L[1][i:i+41])
            TeNodes.add(rev_comp(L[1][i:i+41]))


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