import sys
import json
myfile = open('/tmp/data.json','w+')
location = {}
mutation_sites = []
mutation_sites1 = []
mutation_sites2 = []
def replace_six_characters(original_string, index, replacement_string):
    if len(original_string) >= index + 6 and len(replacement_string) == 6:
        new_string = original_string[:index] + replacement_string + original_string[index + 6:]
        return new_string

    
res_library = {
        'GGGCCC':  'ApaI',
        'GTGCAC':  'ApaLI',
        'ATTAAT':  'AseI',
        'CCTAGG':  'AvrII',
        'GGATCC':  'BamHI', 
        'TGATCA':  'BclI',
        'AGATCT':  'BglII',
        'GCGCGC':  'BssHII',
        'CGGCCG':  'EagI',
        'GAATTC':  'EcorI',
        'GATATC':  'EcoRV',
        'AAGCTT':  'HindIII',
        'GTTAAC':  'HpaI',
        'GGTACC':  'KpnI',
        'ACGCGT':  'M1uI',
        'GCCGGC':  'NaeI',
        'GGCGCC':  'NarI',
        'CATATG':  'NdeI',
        'GCTAGC':  'NheI',
        'TCGCGA':  'NruI',
        'ATGCAT':  'NsiI',
        'CTGCAG':  'PstI',
        'CGATCG':  'PvuI',
        'CAGCTG':  'PvuII',
        'GAGCTC':  'SacI',
        'CCGCGG':  'SacII',
        'GTCGAC':  'SalI',
        'AGTACT':  'ScaI',
        'CCCGGG':  'SmaI',
        'GTATAC':  'SnaI',
        'TACGTA':  'SnaBI',
        'ACTAGT':  'SpeI',
        'GCATGC':  'SphI',
        'AATATT':  'SspI',
        'AGGCCT':  'StuI',
        'TCTAGA':  'XbaI',
        'CTCGAG':  'XhoI'
    }

def add_spaces(input_string):
    result = ""
    for char in input_string:
        result += char + "   "
    return result.strip()  # Remove the trailing space at the end


def format_output(input_string):
    # Convert all characters to uppercase
    input_string = input_string.upper()

    # Add spaces after every three characters
    formatted_output = ' '.join([input_string[i:i+3] for i in range(0, len(input_string), 3)])

    return formatted_output
  
def count_different_characters(str1, str2):
    # Function to count the number of characters that differ between two strings
    if len(str1) != len(str2):
        raise ValueError("The two strings must have the same length.")
    diff_count = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            diff_count += 1
    return diff_count



def replace_substring(original_string, replacement_substring, index):
    if index < 0 or index + len(replacement_substring) > len(original_string):
        print("Error: Index out of bounds.")
        return original_string

    return original_string[:index] + replacement_substring + original_string[index + len(replacement_substring):]




# Step 1: Translate DNA into amino acids
def translate_dna_to_amino_acids(dna_sequence):
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    amino_acids = []
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3].upper()
        amino_acid = codon_table.get(codon, '')
        # if amino_acid == '':
        #     raise ValueError(f"Invalid codon: {codon}")
        amino_acids.append(amino_acid)
    return ''.join(amino_acids)

# Step 2: Cut the input into strings of 6 characters
def cut_into_six_character_strings(dna_sequence):
    return [dna_sequence[i:i+6] for i in range(len(dna_sequence) - 5)]

# Step 3: Run against the res library for a perfect, one mismatch or two mismatches
def find_mutation_site(dna_sequence):
    num_diff = 0
    increment = 1
    six_character_strings = cut_into_six_character_strings(dna_sequence)
  
    for string in six_character_strings:
        mut_string = ""
        mutation_change = 0
        for string2 in res_library:
            num_diff = count_different_characters(string, string2)
            if (num_diff <= 2):
              mut_string = replace_six_characters(dna_sequence, increment - 1, string2)
              mutation_change = count_different_characters(translate_dna_to_amino_acids(mut_string), translate_dna_to_amino_acids(dna_sequence))
              if (mutation_change > 1):
                mutation_sites2.append([string, string2, increment, num_diff, mutation_change])
              elif (mutation_change == 1):
                mutation_sites1.append([string, string2, increment, num_diff, mutation_change])
              else:
                mutation_sites.append([string, string2, increment, num_diff, mutation_change])
        increment = increment + 1



# Main function to process the DNA sequence
def process_dna_sequence(input_dna_sequence):

    translated_amino_acids = translate_dna_to_amino_acids(input_dna_sequence)
    result = f"Input Protein Sequence:\n{translated_amino_acids}\n\n"
    find_mutation_site(input_dna_sequence.upper())

    result += "Restriction Enzyme Sites with Mutations that Preserve Protein Sequence\n"
    result += "Name       Sequence          Site\n"
    result += "========  =========  =====================\n"
    for sequence, enzyme in res_library.items():
        for seq in mutation_sites:
            if seq[1] == sequence:
                result += f"{enzyme:<9} ({sequence}):"
                result += f"  {seq[2]} ({seq[0]}->{seq[1]})\n"
                result += "\n"

    return result

# Input and output file paths
input_file_path = "input_sequence.txt"  
output_file_path = "output_results.txt" 
if len(sys.argv) != 2:
  print("Usage: python findMut.py <input_file>")
  
input_file_path = sys.argv[1]

with open(input_file_path, "r") as input_file:
    p_data = input_file.read()
  
# Read input sequence from file
input_dna_sequence = p_data

# Process the DNA sequence and get the results
processed_results = process_dna_sequence(input_dna_sequence)

print(processed_results)
