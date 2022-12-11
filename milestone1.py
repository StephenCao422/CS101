def s(dna):
    dnaDict = {'A': 0, 'C': 0, 'G': 0, 'T': 0}  # New dictionary
    for x in dna:  # Iterate through the letters and update the dictionary
        dnaDict[x] += 1
    return dnaDict


def dna2rna(dna):
    newRNA = ''  # Empty string to fill a new string of dna
    for x in range(0, len(dna)):  # Iterate through the given dna
        if dna[x] != 'T':  # If not a T, add the letter to the new string unchanged
            newRNA += dna[x]
        else:  # Otherwise change the Ts to Us
            newRNA += 'U'
    return newRNA


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


def mendels_law(hom, het, rec):  # Learn some probability before reading
    x, y, z = hom, het, rec
    term1 = y * z / ((x + y + z) * (x + y + z - 1))  # Probability to choose Aa aa for both orders
    term2 = 0.25 * y * (y - 1) / ((x + y + z) * (x + y + z - 1))  # Probability to choose Aa Aa
    term3 = z * (z - 1) / ((x + y + z) * (x + y + z - 1))  # Probability to choose aa aa
    if y != 0 and z != 0:
        return 1 - term1 - term2 - term3  # Instead of add probabilites for dominants, subtract 100% by the probability of recessives
    elif y == 0:
        return 1 - term1 - term3  # Checks if any input is 0
    elif z == 0:
        return 1 - term1 - term2


def fibonacci_rabbits(n, k):
    if n <= 2:  # Stop reccursing when n <= 1
        return 1  # Value of F1 and F2
    else:
        return k * fibonacci_rabbits(n - 2, k) + fibonacci_rabbits(n - 1,
                                                                   k)  # Sum the previous terms and multiply one of the terms by k


def gc_content(dna_list):
    highest, index = 0, 0

    for i in range(len(dna_list)):  # Iterate through the multiple dna strings
        count = 0
        for j in range(len(dna_list[i])):  # Iterate through each letter of the dna strings
            letter = dna_list[i][j]
            if letter == "C" or letter == "G":  # Add the G and Cs
                count += 1

        percentage = count / len(dna_list[i])  # % of G and Cs in the dna
        if percentage > highest:
            highest = percentage
            index = i

    return index, highest * 100


def rna2codon(rna):
    dic = {
        "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
        "UCU": "S", "UCC": "s", "UCA": "S", "UCG": "S",
        "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
        "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
        "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
        "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }

    protein = ""

    for i in range(0, int(len(rna) / 3)):
        a = rna[3 * i:3 * i + 3]
        if dic[a] != "*":
            protein += dic[a]
        else:
            break

    return protein


def locate_substring(dna_snippet, dna):
    # indexes = [i for i in range(len(dna_snippet))
    # if dna_snippet.startswith(dna, i)]
    indexes = []
    if dna_snippet in dna:
        i = 0
        while i < len(dna):  # Iterate through dna
            newDna = dna[i:]  # Remove previous occurences of snippet from dna
            indexAdj = i  # Adjust index to append correct indexes
            try:
                indexes.append(newDna.index(dna_snippet) + indexAdj)
            except:
                break
            i += newDna.index(dna_snippet) + 1  # Resume after occurence of snippet
    return indexes


def hamming_dist(dna1, dna2):
    i = 0
    flag = 0

    for i in range(len(dna1)):
        if i < len(dna1):
            if dna1[i] != dna2[i]:
                flag = flag + 1
            i = i + 1

    return flag


def count_dom_phenotype(genotypes):
    pairDomCount = {'1': 1, '2': 1, '3': 1, '4': 0.75, '5': 0.5,
                    '6': 0}  # Dictionary for the probabilities of each pair
    offSpring = 0
    for x in range(1, len(genotypes)):  # Loop starts at index one to match keys of dom count
        offSpring += 2 * pairDomCount[str(x)] * genotypes[
            x - 1]  # Dominant offspring product = two babies times the probability from dom count and the count in the list
    return offSpring


def source_rna(protein):
    codonCounts = {'A': 4, 'C': 2, 'D': 2, 'E': 2, 'F': 2, 'G': 4, 'H': 2, 'I': 3, 'K': 2, 'L': 6, 'M': 1, 'N': 2,
                   'P': 4, 'Q': 2, 'R': 6, 'S': 6, 'T': 4, 'V': 4, 'W': 1, 'Y': 2,
                   'Stop': 3}  # Amount of possible amino acids for each codon
    proteinList = list(protein)
    prob = 3  # Multiple for the Stop
    for x in proteinList:  # Multiply the probability by the dictionary multiples for each codon
        prob *= codonCounts[x]
    return prob


def splice_rna(dna, intron_list):
    rna1 = dna2rna(dna)
    for i in intron_list:
        rna1 = rna1.replace(dna2rna(i), "")

    return rna2codon(rna1)