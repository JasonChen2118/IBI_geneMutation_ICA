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


origin_sequence = input("please write the origin DNA sequence")
print("The origin sequence is ", origin_sequence)  # test

result = single_point_mutation(origin_sequence)
print("The mutated sequence is", result[3])
print("Position", result[0], "is mutated from", result[1], "to", result[2])
print("The origin amino acid is ", result[4])
print("The mutated amino acid is", result[5])
