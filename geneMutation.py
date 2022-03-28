# input DNA sequence
# randomly generate a single mutation
# detect the mutation
# output the original gene and mutated gene
# output the result of detecting mutation


import random

# basic given data
DNA_codon = {'TTT': 'Phe', 'TTC': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu',
             'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser',
             'TAT': 'Tyr', 'TAC': 'Tyr', 'TAA': 'stop', 'TAG': 'stop',
             'TGT': 'Cys', 'TGC': 'Cys', 'TGA': 'stop', 'TGG': 'Trp',
             'CTT': 'Leu', 'CTC': 'Leu', 'CTA': 'Leu', 'CTG': 'Leu',
             'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
             'CAT': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
             'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
             'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile', 'ATG': 'Met',
             'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
             'AAT': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
             'AGT': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
             'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val',
             'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
             'GAT': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
             'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'}
base = ['A', 'G', 'C', 'T']


# Task 1
# define a function to make single point mutation
def single_point_mutation(sequence):
    while True:
        mutate_point = random.randint(1, len(sequence))  # randomly generate a mutate point
        origin_base = sequence[mutate_point - 1]
        mutated_base = base[random.randint(0, 3)]
        if origin_base != mutated_base:
            break
            # This make sure that the base isn't changed to itself, which means no mutation.

    mutated_sequence = sequence[:(mutate_point - 1)] + mutated_base + sequence[mutate_point:]

    # The following part generate the codon.
    # The remainder of mutate_point/3 was classified in order to find out the codon which include 3 base pairs.
    origin_amino_acid = None
    mutated_amino_acid = None
    if mutate_point % 3 == 1:
        origin_codon = sequence[(mutate_point - 1): (mutate_point + 2)]
        mutated_codon = mutated_base + sequence[mutate_point: (mutate_point + 2)]
        origin_amino_acid = DNA_codon[origin_codon]
        mutated_amino_acid = DNA_codon[mutated_codon]

    elif mutate_point % 3 == 2:
        origin_codon = sequence[(mutate_point - 2): (mutate_point + 1)]
        mutated_codon = sequence[mutate_point - 2] + mutated_base + sequence[mutate_point]
        origin_amino_acid = DNA_codon[origin_codon]
        mutated_amino_acid = DNA_codon[mutated_codon]

    elif mutate_point % 3 == 0:
        origin_codon = sequence[(mutate_point - 3): mutate_point]
        mutated_codon = sequence[(mutate_point - 3): (mutate_point - 1)] + mutated_base
        origin_amino_acid = DNA_codon[origin_codon]
        mutated_amino_acid = DNA_codon[mutated_codon]

    return mutate_point, origin_base, mutated_base, mutated_sequence, origin_amino_acid, mutated_amino_acid
    #      result[0]     result[1]    result[2]     result[3]         result[4]          result[5]
    # The value is returned as tuple


# Task 2 and 3
def single_mutation_count(sequence):
    total_count = 0
    synonymous_count = 0
    nonsynonymous_count = 0
    truncating_count = 0
    temp = 0
    mutated_sequence = sequence

    for i in range(len(mutated_sequence)):
        temp += 1

        for j in range(4):
            if sequence[i] != base[j]:
                mutated_sequence = origin_sequence[:i] + base[j] + origin_sequence[(i + 1):]
                origin_codon = sequence[(i + 1 - temp): (i + 4 - temp)]
                mutated_codon = mutated_sequence[(i + 1 - temp): (i + 4 - temp)]
                origin_amino_acid = DNA_codon[origin_codon]
                mutated_amino_acid = DNA_codon[mutated_codon]

                if origin_amino_acid == mutated_amino_acid:
                    synonymous_count += 1
                if origin_amino_acid != mutated_amino_acid:
                    nonsynonymous_count += 1
                if (mutated_amino_acid == "stop") and ((i - temp + 4) != len(sequence)):
                    truncating_count += 1
                total_count += 1

        if temp == 3:
            temp = 0

    return total_count, synonymous_count, nonsynonymous_count, truncating_count


# Task 4
def single_transition_mutation(sequence):
    temp = 0
    amino_acid_changes = 0
    total_transition_mutation = 0
    mutated_sequence = sequence

    for i in range(len(mutated_sequence)):
        temp += 1
        if sequence[i] == "A":
            mutated_sequence = origin_sequence[:i] + "G" + origin_sequence[(i + 1):]
        if sequence[i] == "G":
            mutated_sequence = origin_sequence[:i] + "A" + origin_sequence[(i + 1):]
        if sequence[i] == "C":
            mutated_sequence = origin_sequence[:i] + "T" + origin_sequence[(i + 1):]
        if sequence[i] == "T":
            mutated_sequence = origin_sequence[:i] + "C" + origin_sequence[(i + 1):]

        origin_codon = sequence[(i + 1 - temp): (i + 4 - temp)]
        mutated_codon = mutated_sequence[(i + 1 - temp): (i + 4 - temp)]
        origin_amino_acid = DNA_codon[origin_codon]
        mutated_amino_acid = DNA_codon[mutated_codon]

        total_transition_mutation += 1
        if origin_amino_acid != mutated_amino_acid:
            amino_acid_changes += 1

        if temp == 3:
            temp = 0

    return total_transition_mutation, amino_acid_changes


def single_transversion_mutation(sequence):
    temp = 0
    amino_acid_changes = 0
    total_transversion_mutation = 0
    origin_amino_acid = None
    mutated_amino_acid = None
    mutated_sequence = sequence

    for i in range(len(mutated_sequence)):
        temp += 1

        if (origin_sequence[i] == base[0]) or (origin_sequence[i] == base[1]):
            for j in range(2, 4):
                mutated_sequence = origin_sequence[:i] + base[j] + origin_sequence[(i + 1):]
                origin_codon = sequence[(i + 1 - temp): (i + 4 - temp)]
                mutated_codon = mutated_sequence[(i + 1 - temp): (i + 4 - temp)]
                origin_amino_acid = DNA_codon[origin_codon]
                mutated_amino_acid = DNA_codon[mutated_codon]

                total_transversion_mutation += 1

        if (origin_sequence[i] == base[2]) or (origin_sequence[i] == base[3]):
            for j in range(2):
                mutated_sequence = origin_sequence[:i] + base[j] + origin_sequence[(i + 1):]
                origin_codon = sequence[(i + 1 - temp): (i + 4 - temp)]
                mutated_codon = mutated_sequence[(i + 1 - temp): (i + 4 - temp)]
                origin_amino_acid = DNA_codon[origin_codon]
                mutated_amino_acid = DNA_codon[mutated_codon]

                total_transversion_mutation += 1

        if origin_amino_acid != mutated_amino_acid:
            amino_acid_changes += 1

        if temp == 3:
            temp = 0

    return total_transversion_mutation, amino_acid_changes


# Task 5
# Insertion mutation
def amino_acid_sequence(sequence):
    scan = 0
    amino_acid_sequence = []
    while True:
        codon = sequence[scan: (scan + 3)]
        amino_acid = DNA_codon[codon]
        amino_acid_sequence.append(amino_acid)
        scan += 3
        if (scan + 3) > len(sequence):
            break

    return amino_acid_sequence


def insertion_mutation(sequence):
    insertion_point = random.randint(0, len(sequence))
    insert_base = base[random.randint(0, 3)]
    mutated_sequence = sequence[0: insertion_point] + insert_base + sequence[insertion_point:]
    mutated_amino_acid_sequence = amino_acid_sequence(mutated_sequence)
    return mutated_sequence, mutated_amino_acid_sequence


origin_sequence = input("please write the origin DNA sequence")
print("The origin sequence is ", origin_sequence)  # test

# Task 1 result
result = single_point_mutation(origin_sequence)
print("The mutated sequence is", result[3])
print("Position", result[0], "is mutated from", result[1], "to", result[2])
print("The origin amino acid is ", result[4])
print("The mutated amino acid is", result[5])

# Task 2 and 3 result
result = single_mutation_count(origin_sequence)
print(result[0], "possible single point mutation(s)")
print(result[1], "synonymous mutation(s)")
print(result[2], "nonsynonymous mutation(s)")
print(result[3], "truncating mutation(s)")
truncating_fraction = result[3] / result[0]
print("The fraction of truncating mutation is", truncating_fraction)

# Task 4
transition_result = single_transition_mutation(origin_sequence)
transversion_result = single_transversion_mutation(origin_sequence)
print("There are(is)", transition_result[0], "possible transition mutation(s).")
print(transition_result[1], "of them change the amino acid.")
print("There are(is)", transversion_result[0], "possible transversion mutation(s).")
print(transversion_result[1], "of them change the amino acid.")

# Task 5
insertion_result = insertion_mutation(origin_sequence)
print("The origin sequence is ", origin_sequence)
print("The mutated sequence is", insertion_result[0])
print("The origin amino acid sequence is ", amino_acid_sequence(origin_sequence))
print("The mutated amino acid sequence is", insertion_result[1])
