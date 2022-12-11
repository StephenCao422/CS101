def find_splice(dna_motif, dna):
    newDna = dna
    foundList = []
    indexAdj = 0
    for x in dna_motif:
        try:
            found = newDna.index(x)  # Check if x of dna_motif is found
        except:
            return []  # Return empty list if full motif not found
        newDna = newDna[found + 1:]  # Shorten the list
        # So only first occurence after each letter in motif exists
        foundList.append(found + indexAdj)
        indexAdj += found + 1  # Adjust the index for the shortening
    return foundList


def shared_motif(dna_list):
    maxList = max(dna_list, key=len)  # Identify the list with max length
    stepTot = 0
    index = 0
    maybeLong = []
    for x in range(1, len(maxList)):  # Iterate max len list, change step length x
        while stepTot < len(maxList):  # Iterate through the max len list with stepping
            yes = 0
            for i in range(0, len(dna_list)):  # Iterate through each given string
                if maxList[index: index + x] in dna_list[i]:  # Check substring in strings
                    yes += 1
            if yes == len(dna_list):  # Check if substring in every string
                maybeLong.append(maxList[index: index + x])
            if index + len(maxList[index: index + x]) >= len(maxList):  # break when string longer than max string
                break
            index += 1  # Change starting index
            stepTot += 1  # Increase the amout of steps taken
        stepTot = 0
        index = 0
    try:
        maxSub = max(maybeLong, key=len)
    except ValueError:
        return ''
    return maxSub


def get_edges(dna_dict):
    match = []
    for x in dna_dict:
        for i in dna_dict:  # Check each item with each item
            if [x, i] in match or [i, x] in match or x == i:  # Check repeats
                continue
            if dna_dict[x][:3] == dna_dict[i][-3:] or dna_dict[x][-3:] == dna_dict[i][:3]:  # Checking if beginning in end or end in beginning
                match.append([x, i])
    return match


from itertools import permutations


def assemble_genome(dna_list):
    def fcombine(one, two):
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
    ecombines = []
    for pair in permutations(dna_list, len(dna_list)):
        pairs.append(list(pair))
    for x in pairs:
        pairsCopy.append(list(x))
    for x in range(0, len(pairs)):
        fcombines.append(fccombine(pairs[x]))
        fmin = min(fcombines, key=len)
    for x in pairsCopy:
        ecombines.append(eccombine(x))
        emin = min(ecombines, key=len)

    return min(fmin, emin)


import math


def perfect_match(rna):
    d = {'A': 0, 'C': 0, 'G': 0, 'U': 0}
    for i in rna:
        d[i] += 1
    if d['A'] == d['U'] and d['C'] == d['G']:  # Find a 'A' for each 'U' and a 'G' for each 'C'
        return math.factorial(d['A']) * math.factorial(
            d['C'])  # Total number = Factorial(number of A) * Factorial(number of C)
    else:
        return 0


def random_genome(dna, gc_content):
    dna = dna.upper()
    CG = 0
    AT = 0
    for i in dna:
        if i == 'C':
            CG += 1
        elif i == 'G':
            CG += 1
        elif i == 'A':
            AT += 1
        elif i == 'T':
            AT += 1
    result = []
    for j in range(0, len(gc_content)):
        probabilities = CG * math.log10(float(gc_content[j]) / 2) + AT * math.log10((1 - float(gc_content[j])) / 2)
        result.append(probabilities)
    return result


def reverse_complement(dna):
    complement = ''
    for x in range(len(dna) - 1, -1, -1):  # Iterate from last index and add new values to beginning index
        if dna[x] == 'A':
            complement += 'T'
        elif dna[x] == 'T':
            complement += 'A'
        elif dna[x] == 'C':
            complement += 'G'
        else:
            complement += 'C'
    return complement


def rev_palindrome(dna):
    result = []
    for i in range(len(dna) - 4):  # Go through the entire string length
        for j in range(i + 3, min(len(dna), i + 12)):  # Go through from string from length 4 to 12
            string = dna[i:j + 1]
            if reverse_complement(dna[i:j + 1]) == string:  # check if the string compliment is a palindrome
                result.append((i, j - i + 1))
    return result
