from itertools import permutations


f = open("ms3-dna-mammuthus.txt", "r")
file = f.read()
file_s = file.split("\n")
dna_list = []
for x in file_s:
    dna_list.append(x)


def assemble_genome(dna_list):
    def fcombine(one , two):
        combinations = []
        for j in range(0, len(one)):
            if one[-j:] == two[:j]:
                combinations.append(one + two[j:])
        if combinations == []:
            return one + two
        return min(combinations, key=len)

    def ecombine(one, two):
        combinations = []
        for j in range(0, len(one)):
            if one[:j] == two[-j:]:
                combinations.append(two + one[j:])
        if combinations == []:
            return two + one
        return min(combinations, key=len)

    def fccombine(dnaList):
        combinations = []
        dnaListN = dnaList
        x = 0
        while x < len(dnaListN) - 1:
            combinations.append(fcombine(dnaListN[x], dnaListN[x + 1]))
            dnaListN[0] = combinations[-1]
            dnaListN.pop(x + 1)
        return combinations[-1]

    def eccombine(dnaList):
        combinations = []
        dnaListN2 = dnaList
        x = 0
        while x < len(dnaListN2) - 1:
            combinations.append(ecombine(dnaListN2[x], dnaListN2[x + 1]))
            dnaListN2[0] = combinations[-1]
            dnaListN2.pop(x + 1)
        return combinations[-1]

    pairs = []
    pairsCopy = []
    fcombines = []
    #ecombines = []
    for pair in permutations(dna_list, len(dna_list)):
        pairs.append(list(pair))
    for x in pairs:
        pairsCopy.append(list(x))
    for x in range(0, len(pairs)):
        fcombines.append(fccombine(pairs[x]))
        fmin = min(fcombines, key=len)
    #for x in pairsCopy:
        #ecombines.append(eccombine(x))
        #fmin = min(ecombines, key = len)

    return fmin

# if S1 == S2:
# do nothing, continue to next string

# S1 + S2:
# i: S1: 0; j: S2:len(S1)-1 (if len(S2) > len(S1))
# i: S1: len(S1) - len(S2); j: S2:len(S2)-1 (if len(S2) <= len(S1))
# tempString 1

# S2 + S1:
# i: S2: 0; j: S1: len(S2)-1 (if len(S2) <= len(S1))
# i: S1:len(S1)-1; j: S2:len(S2) - len(S1) (if len(S2) > len(S1))
# tempString 2

# stringToKeep: tempString1 if len(tempString1) < len(tempString2) else tempString2
