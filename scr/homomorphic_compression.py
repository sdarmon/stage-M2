import sys

Arg = sys.argv[:]

def homomorphic_compression(read):
    #Removing consecutive duplicates characters
    S = read[0]
    for i in range(1, len(read)):
        if read[i] != read[i-1]:
            S += read[i]
    return S

if len(Arg) not in [3]:
    print("Use : " + Arg[0] + " reads.fq reads_out.fq")
    exit()

cycle = 0
if len(Arg) == 3:
    with open(Arg[1], 'r') as f:
        with open(Arg[2], 'w') as f_out:
            for line in f:
                cycle = (cycle + 1) % 4
                if len(line) > 2:
                    if cycle != 2:
                        f_out.write(line)
                        continue
                    f_out.write(homomorphic_compression(line))