
import random
N = int(2)
K = int(4)
F = K//4

base=""
for i in range(0,K):
    a = random.randint(0,3)
    if a == 0:
        base+="A"
    elif a == 1:
        base+="T"
    elif a == 2:
        base+="G"
    elif a == 3:
        base+="C"


genome=[]
genome.append(base)

for j in range(0,N):
    b1= genome[j]
    # print(b1)
    for st in ["A", "T","G","C"]:
        genome.append(st+b1[:K-1])
        genome.append(b1[1:]+st)

# print(genome)
print(len(genome))
for st in sorted(genome):
    print(st)