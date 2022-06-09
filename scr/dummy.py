import sys

Arg = sys.argv[:]
seq = []
vu = set()
compt = 0
step = 0
with open(Arg[1], 'r') as f:
    for line in f:
        step = (step+1)%4
        if step == 0:
            vu.add(line)
            compt+=1
print("Il y a "+str(len(vu))+" under paths unqiues sur les "+str(compt)+" bulles.")
with open(Arg[2], 'w') as o:
    for el in vu:
        o.write(el)