import time
start = time.time()


f = open("ms3-dna-mammuthus.txt", "r")
file = f.read()
file_s = file.split("\n")
dna_list = []
for x in file_s:
    dna_list.append(x)


import itertools
def overlap(str_1,str_2, min_lengh):
    start = 0
    while True:
        start = str_1.find(str_2[:min_lengh], start)
        if start == -1:
            return 0
        if str_2.startswith(str_1[start:]):
            return len(str_1)-start
        start += 1

def pick_maximal_overlap(reads,k):
    reada, readb = None, None
    best_olen = 0
    for a,b in itertools.permutations(reads,2):
        olen = overlap(a,b,min_lengh=k)
        if olen > best_olen:
            reada,readb = a,b
            best_olen = olen
    return reada,readb, best_olen

def greedy_scs(reads,k):
    read_a, read_b, olen = pick_maximal_overlap(reads,k)
    while olen > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a+read_b[olen:])
        read_a , read_b, olen = pick_maximal_overlap(reads,k)
    return ''.join(reads)


read = ["ATTAGACCTG","CCTGCCGGAA","AGACCTGCCG","GCCGGAATAC"]
print(greedy_scs(read,1))
end = time.time()
print("The time of execution of above program is :",(end-start), "s")
