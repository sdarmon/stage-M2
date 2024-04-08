import sys

Arg = sys.argv[:]

def homomorphic_compression(read,w):
    #Removing consecutive duplicates window by window
    S = read[0:w]
    current_w = read[0:w]
    for i in range(1, len(read)-w):
        next_w = read[i:i+w]
        if current_w != next_w:
            S += read[i+w-1]
            current_w = next_w
    return S

if len(Arg) not in [4]:
    print("Use : " + Arg[0] + " reads.fq reads_out.fq window")
    exit()
