# input DNA sequence
# randomly generate a single mutation
# detect the mutation
# output the original gene and mutated gene
# output the result of detecting mutation


import random

DNACodon = {'TTT': 'Phe', 'TTC': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu',
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


def single_point_mutation(sequence):
    mutate_point = random.randint(1, len(sequence))
    origin_base = sequence[mutate_point - 1]
    mutated_base = base[random.randint(0, 3)]

    mutated_sequence = sequence[:(mutate_point - 1)] + mutated_base + sequence[mutate_point:]
    return mutate_point, origin_base, mutated_base, mutated_sequence


origin_sequence = input("please write the origin DNA sequence")
print("The origin sequence is ", origin_sequence)  # test

mutated_sequence = single_point_mutation(origin_sequence)
print("The mutated sequence is", mutated_sequence[3])
print("Position", mutated_sequence[0], "is mutated from", mutated_sequence[1], "to", mutated_sequence[2])
