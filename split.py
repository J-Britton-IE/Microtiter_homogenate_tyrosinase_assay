# Central dogma functions library
import numpy as np
import pandas as pd
import matplotlib as plt
import scipy
from scipy import stats

# 2) check if sequence has correct characters
def valid(sequence, valid_chars):

    for item in sequence:
        if item in valid_chars:
            pass
        elif item not in valid_chars:
            return True
        else:
            return False


# 3) Identify if sequence is DNA, RNA or protein
def identify(user_sequence):

    DNA_chars = ["A", "T", "C", "G"]
    RNA_chars = ["U"]
    Amino_acid_chars = ["A", "R", "N", "D", "B", "C", "Q", "E", "Z", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W",
                        "Y", "V"]
    not_compatible = ["T", "U"]

    sequence_identity = []

    if (set(not_compatible).issubset(set(user_sequence))):
        sequence_identity += "N"

    for item in user_sequence:
        if item in Amino_acid_chars:
            if item not in DNA_chars:
                sequence_identity += "P"
            elif item in DNA_chars:
                sequence_identity += "D"
        elif item in RNA_chars:
            sequence_identity += "R"
        else:
            sequence_identity += "N"

    sequence_identity = "".join(sequence_identity)

    identity = ""

    Protein = "P"
    DNA = "D"
    RNA = "R"
    Non = "N"

    if sequence_identity.find(Non) != -1:
        identity = "Not a valid sequence"

    elif sequence_identity.find(RNA) != -1:

        if sequence_identity.find(Protein) != -1:
            identity = "Not a valid sequence"

        else:
            identity = "RNA"

    elif sequence_identity.find(Protein) != -1:
        identity = "Protein"

    elif sequence_identity.find(DNA) != -1:
        identity = "DNA"

    return identity


# 4) Transcribe DNA into RNA
def transcription(user_sequence):

    user_sequence = list(user_sequence)

    rna_Sequence = []

    for item in user_sequence:
        if item == 'T':
            rna_Sequence += 'U'
        else:
            rna_Sequence += item

    rna_Sequence = ''.join(rna_Sequence)

    return rna_Sequence


# 5) Translate RNA into Amino acid sequence
def translate(user_sequence):
    triplets = 3
# Separate sequence to triplets
    rna_triplets = [user_sequence[i:i+triplets] for i in range(0, len(user_sequence), triplets)]

    amino_acid_rna_dict = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "stop",
    "UAG": "stop",
    "UGU": "C",
    "UGC": "C",
    "UGA": "stop",
    "UGG": "W",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
    }

    amino_acid_list = []

    for item in rna_triplets:
        for key, value in amino_acid_rna_dict.items():
            if item == key:
                amino_acid_list.append(value)

    amino_acid_list = "".join(amino_acid_list)

    return amino_acid_list


# 6) Amino acid molecular weight (kDa)
def mol_weight(user_sequence):

    AA_mol_weight_dict = {
    "A": 0.08909,
    "R": 0.1742,
    "N": 0.1321,
    "D": 0.1331,
    "B": 0.133,
    "C": 0.1212,
    "Q": 0.1461,
    "E": 0.1471,
    "Z": 0.147,
    "G": 0.07507,
    "H": 0.1552,
    "I": 0.1312,
    "L": 0.1312,
    "K": 0.1462,
    "M": 0.1492,
    "F": 0.1652,
    "P": 0.1151,
    "S": 0.10509,
    "T": 0.1191,
    "W": 0.2042,
    "Y": 0.1812,
    "V": 0.1171,
}

    user_sequence = list(user_sequence)

    mol_weight_list = []

    for item in user_sequence:
        for key, value in AA_mol_weight_dict.items():
            if item == key:
                mol_weight_list.append(value)

    mol_weight_list = sum(mol_weight_list)

    mol_weight_list = round(mol_weight_list, 2)

    return str(mol_weight_list)


# 7) Reverse transcribe RNA to DNA
def r_transcribe(user_sequence):

    user_sequence = list(user_sequence)

    dna_Sequence = []

    for item in user_sequence:
        if item == 'U':
            dna_Sequence += 'T'
        else:
            dna_Sequence += item

    dna_Sequence = ''.join(dna_Sequence)

    return dna_Sequence


# 8) Reverse translate Protein to RNA
def r_translate(user_sequence):

    user_sequence = list(user_sequence)

    amino_acid_rna_dict_e_coli = {
        "UUU": "F",
        "UUA": "L",
        "UCU": "S",
        "UAU": "Y",
        "UAA": "stop",
        "UGC": "C",
        "UGA": "stop",
        "UGG": "W",
        "CUG": "L",
        "CCG": "P",
        "CAU": "H",
        "CAG": "Q",
        "CGU": "R",
        "AUU": "I",
        "AUG": "M",
        "ACC": "T",
        "AAC": "N",
        "AAA": "K",
        "AGC": "S",
        "GUU": "V",
        "GCG": "A",
        "GAU": "D",
        "GAA": "E",
        "GGC": "G",
    }

    rna_sequence = []

    for item in user_sequence:
        for key, value in amino_acid_rna_dict_e_coli.items():
            if item == value:
                rna_sequence.append(key)

    rna_sequence = "".join(rna_sequence)

    return rna_sequence


# 9) GC content calculator
def gc(user_sequence):

    sequence_length = len(user_sequence)

    user_sequence = list(user_sequence)

    gsandcs = ["G", "C"]

    gccount = []

    for item in user_sequence:
        if item in gsandcs:
            gccount += item

    gccount = "".join(gccount)

    gclen = len(gccount)

    gc_percent = (gclen / sequence_length)*100

    gc_percent = round(gc_percent, 2)

    return gc_percent


# 10) Flip string
def flip(user_sequence):

    return user_sequence[::-1]


# 11) DNA Complementary strand
def complement(user_sequence):

    user_sequence = list(user_sequence)

    complement_sequence = []

    for item in user_sequence:
        if item == "T":
            complement_sequence += "A"
        elif item == "A":
            complement_sequence += "T"
        elif item == "G":
            complement_sequence += "C"
        elif item == "C":
            complement_sequence += "G"

    complement_sequence = "".join(complement_sequence)

    return complement_sequence


# 12) Reverse complement
def reverse_complement(user_sequence):

    user_sequence = flip(user_sequence)

    user_sequence = list(user_sequence)

    reverse_complement_sequence = []

    for item in user_sequence:
        if item == "T":
            reverse_complement_sequence += "A"
        elif item == "A":
            reverse_complement_sequence += "T"
        elif item == "G":
            reverse_complement_sequence += "C"
        elif item == "C":
            reverse_complement_sequence += "G"

    reverse_complement_sequence = "".join(reverse_complement_sequence)

    return reverse_complement_sequence


# 13) Count number of each base in sequence
def base_counts(user_sequence):

    user_sequence = list(user_sequence)

    tees = []
    gees = []
    cees = []
    aaas = []

    for item in user_sequence:
        if item == "T":
            tees += "T"
        elif item == "G":
            gees += "G"
        elif item == "C":
            cees += "C"
        elif item == "A":
            aaas += "A"

    tees = "".join(tees)
    gees = "".join(gees)
    cees = "".join(cees)
    aaas = "".join(aaas)

    tees = len(tees)
    gees = len(gees)
    cees = len(cees)
    aaas = len(aaas)

    tees = [tees]
    gees = [gees]
    cees = [cees]
    aaas = [aaas]

    joins = aaas + tees + cees + gees

    return joins

# 14) Codon optimisation e.coli - replaces low frequency with high frequency codons.
def e_coli_rna_codon_optimisation(sequence):

    e_coli_codon_optimisation_dict = {
        "A": {"B": '1'},
        "UUU": {"F": '1.9'},
        "UUC": {"F": '1.8'},
        "UUA": {"L": '1'},
        "UUG": {"L": '1.1'},
        "UCU": {"S": '1.1'},
        "UCC": {"S": '1'},
        "UCA": {"S": '0.7'},
        "UCG": {"S": '0.8'},
        "UAU": {"Y": '1.6'},
        "UAC": {"Y": '1.4'},
        "UAA": {"stop": '0.2'},
        "UAG": {"stop": '0.03'},
        "UGU": {"C": '0.4'},
        "UGC": {"C": '0.6'},
        "UGA": {"stop": '0.1'},
        "UGG": {"W": '1.4'},
        "CUU": {"L": '1'},
        "CUC": {"L": '0.9'},
        "CUA": {"L": '0.3'},
        "CUG": {"L": '5.2'},
        "CCU": {"P": '0.7'},
        "CCC": {"P": '0.4'},
        "CCA": {"P": '0.8'},
        "CCG": {"P": '2.4'},
        "CAU": {"H": '1.2'},
        "CAC": {"H": '1.1'},
        "CAA": {"Q": '1.3'},
        "CAG": {"Q": '2.9'},
        "CGU": {"R": '2.4'},
        "CGC": {"R": '2.2'},
        "CGA": {"R": '0.3'},
        "CGG": {"R": '0.5'},
        "AUU": {"I": '2.71'},
        "AUC": {"I": '2.7'},
        "AUA": {"I": '0.4'},
        "AUG": {"M": '2.6'},
        "ACU": {"T": '1.2'},
        "ACC": {"T": '2.4'},
        "ACA": {"T": '0.1'},
        "ACG": {"T": '1.3'},
        "AAU": {"N": '1.6'},
        "AAC": {"N": '2.6'},
        "AAA": {"K": '3.8'},
        "AAG": {"K": '1.2'},
        "AGU": {"S": '0.7'},
        "AGC": {"S": '1.5'},
        "AGA": {"R": '0.2'},
        "AGG": {"R": '0.2'},
        "GUU": {"V": '2'},
        "GUC": {"V": "1.4"},
        "GUA": {"V": '1.2'},
        "GUG": {"V": '2.4'},
        "GCU": {"A": '1.8'},
        "GCC": {"A": '2.3'},
        "GCA": {"A": '2.1'},
        "GCG": {"A": '3.2'},
        "GAU": {"D": '3.3'},
        "GAC": {"D": '2.3'},
        "GAA": {"E": '4.4'},
        "GAG": {"E": '1.9'},
        "GGU": {"G": '2.8'},
        "GGC": {"G": '3'},
        "GGA": {"G": '0.7'},
        "GGG": {"G": '0.9'}
    }

    amino_acid_rna_dict = {
        "UUU": "F",
        "UUC": "F",
        "UUA": "L",
        "UUG": "L",
        "UCU": "S",
        "UCC": "S",
        "UCA": "S",
        "UCG": "S",
        "UAU": "Y",
        "UAC": "Y",
        "UAA": "-",
        "UAG": "-",
        "UGU": "C",
        "UGC": "C",
        "UGA": "-",
        "UGG": "W",
        "CUU": "L",
        "CUC": "L",
        "CUA": "L",
        "CUG": "L",
        "CCU": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "CAU": "H",
        "CAC": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGU": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "AUU": "I",
        "AUC": "I",
        "AUA": "I",
        "AUG": "M",
        "ACU": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "AAU": "N",
        "AAC": "N",
        "AAA": "K",
        "AAG": "K",
        "AGU": "S",
        "AGC": "S",
        "AGA": "R",
        "AGG": "R",
        "GUU": "V",
        "GUC": "V",
        "GUA": "V",
        "GUG": "V",
        "GCU": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "GAU": "D",
        "GAC": "D",
        "GAA": "E",
        "GAG": "E",
        "GGU": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
    }

    triplets = 3

    sequence = [sequence[i:i + triplets] for i in range(0, len(sequence), triplets)]  # split sequence to triplets.

    aas = []

    # t item for item in user list
    for item in sequence:
        for key, value in amino_acid_rna_dict.items():
            if item == key:
                aas += value

    amino_acid__optimal_rna_dict_e_coli = {
        "UUU": "F",
        "UUA": "L",
        "UCU": "S",
        "UAU": "Y",
        "UAA": "-",
        "UGC": "C",
        "UGG": "W",
        "CUG": "L",
        "CCG": "P",
        "CAU": "H",
        "CAG": "Q",
        "CGU": "R",
        "AUU": "I",
        "AUG": "M",
        "ACC": "T",
        "AAC": "N",
        "AAA": "K",
        "AGC": "S",
        "GUU": "V",
        "GCG": "A",
        "GAU": "D",
        "GAA": "E",
        "GGC": "G",
    }

    optimised = []

    for item in aas:
        for key, value in amino_acid__optimal_rna_dict_e_coli.items():
            if item == value:
                optimised.append(key)

    optimised = ''.join(optimised)

    return optimised

##Issue with stop codons -

# 15) Rate calculation sliding window function.

def rate_func(x, y, width):

    n = len(x);

    m = len(y.columns);
    
    steps = n - width + 1;

    rates_coefficients_m    =   np.zeros((steps, m),dtype = int); # Make vector of size  (steps x  number of columns)
    
    rates_coefficients_b    =   np.zeros((steps, m), dtype = int); #Make vector of size  (steps x  number of columns)

    rates_standardError     =   np.zeros((steps, m), dtype = int);
    
    rates_steps             =   steps;

    slopes = [];
    
    c_vals = [];

    for j in range(0, steps + 1): # Iterate through every row
    
        poly_fit = np.polyfit(x[j:(j + width -1)], y[j: (j + width -1)],1); # Calculates slopes in small sets
        
        slopes.append(poly_fit[0]);
        c_vals.append(poly_fit[1]);

    x_plot_vals = x;

    x_plot_vals = x_plot_vals[width -2:]

    x_plot_vals = pd.DataFrame(x_plot_vals);
    
    slopes = pd.DataFrame(slopes);

    rate_data = pd.concat([x_plot_vals, slopes], axis = 1);
    
    return rate_data

# r2 value

def rsquared(x,y):
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x,y)
    return r_value**2
