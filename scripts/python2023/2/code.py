# Define a function to convert a codon to binary bits
def codon_to_binary(codon):
    # Determine the binary bits for each nucleotide in the codon
    bit_2 = (codon[1] == 'G' or codon[1] == 'C')  # 1 if G or C, 0 if A or U
    bit_1 = (codon[0] == 'G' or codon[0] == 'C')  # 1 if G or C, 0 if A or U
    bit_3 = (codon[2] == 'G' or codon[2] == 'C')  # 1 if G or C, 0 if A or U
    # Combine the binary bits into a 6-bit binary index
    binary_index = (bit_2 << 2) | (bit_1 << 1) | bit_3  # Left-shift bit_2 by 2, left-shift bit_1 by 1, and OR them with bit_3
    # Return the binary index as a 6-bit Boolean list
    return [bool(int(bit)) for bit in format(binary_index, '06b')]

# Define a function to convert a nucleotide sequence to binary code
def nucleotide_to_binary(nucleotide_sequence):
    # Divide the nucleotide sequence into groups of three
    codons = [nucleotide_sequence[i:i+3] for i in range(0, len(nucleotide_sequence) - len(nucleotide_sequence) % 3, 3)]
    # Convert each codon to a binary index and concatenate the results
    binary_indices = [codon_to_binary(codon) for codon in codons]
    binary_code = [bit for index in binary_indices for bit in index]
    # Return the final binary code as a Boolean list
    return binary_code

# Define the nucleotide sequence
nucleotides = "ATGCGTCGAGCTCGATCGTTAGCTCGATCGACGATCGTCG"

# Ensure that the sequence length is a multiple of three
if len(nucleotides) % 3 != 0:
    nucleotides = nucleotides[:-(len(nucleotides) % 3)]

# Convert the nucleotide sequence to binary code
binary_code = nucleotide_to_binary(nucleotides)

# Print the binary code as a string
binary_string = ''.join(str(int(bit)) for bit in binary_code)
print(binary_string)
