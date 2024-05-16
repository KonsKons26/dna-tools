
def reverse_string(seq):
    """This function returns the reverse of a given DNA sequence."""
    return seq[::-1]


def complement(seq):
    """This function returns the complement of a given DNA sequence.
    Base pairs:
    A: adenosine
    C: cytosine
    G: guanine
    T: thymine
    N: A or C or G or T
    M: A or C
    R: A or G
    W: A or T
    S: C or G
    Y: C or T
    K: G or T
    V: A or C or G; not T
    H: A or C or T; not G
    D: A or G or T; not C
    B: C or G or T; not A
    """
    complement_bases = {
        'A': 'T', 'a': 't', 'T': 'A', 't': 'a', 'C': 'G', 'c': 'g', 'G': 'C', 'g': 'c',
        'N': 'N', 'n': 'n', 'M': 'K', 'm': 'k', 'K': 'M', 'k': 'm', 'R': 'Y', 'r': 'y',
        'Y': 'R', 'y': 'r', 'W': 'W', 'w': 'w', 'S': 'S', 's': 's', 'V': 'B', 'v': 'b',
        'B': 'V', 'b': 'v', 'H': 'D', 'h': 'd', 'D': 'H', 'd': 'h'
    }
    letters = list(seq)
    letters = [complement_bases[base] for base in letters]
    return ''.join(letters)


def reverse_complement(seq):
    """This function returns the reverse complement of a given DNA sequence."""
    seq = reverse_string(seq)
    seq = complement(seq)
    return seq


def count_base(dna, base):
    """This function counts the given base in a dna string."""
    return dna.count(base)


def count_all_letters(dna):
    """This function counts all letter occurrences."""
    letters = set(dna)
    return {letter: count_base(dna, letter) for letter in letters}


def create_dna(n, alphabet='acgt'):
    """This function generates a random DNA sequence."""
    import random
    return ''.join([random.choice(alphabet) for i in range(n)])


def rmv_new_lines_spaces(seq):
    """This function removes all whitecpases of a given DNA sequence."""
    return ''.join(seq.split())


def gc(dna):
    """This function computes the GC percentage of a DNA sequence."""
    nbases = dna.count('N') + dna.count('n')
    gcpercent = float(dna.count('G') + dna.count('g') + dna.count('C') +
                      dna.count('c')) * 100 / (len(dna)-nbases)
    return gcpercent


def get_all_frames(dna):
    """This function returns a dict with all reading frames (froward and reverse)."""
    rev_comp = reverse_complement(dna)
    fframes = [dna[i:] for i in range(0, 3)]
    rframes = [rev_comp[i:] for i in range(0, 3)]
    frames = fframes + rframes
    frame_names = ["f1", "f2", "f3", "r1", "r2", "r3"]
    return {frame_name: frame for frame_name,
            frame in zip(frame_names, frames)}


def has_stop_codon(dna, frame=0):
    """This function checks if a given DNA sequence has an in frame stop codon."""
    stop_codon_found = False
    stop_codons = ['tga', 'tag', 'taa']
    for i in range(frame, len(dna), 3):
        codon = dna[i: i + 3].lower()
        if codon in stop_codons:
            stop_codon_found = True
            break
    return stop_codon_found


def find_stop_codons(dna, frame=0):
    """This function returns the indeces of all in frame stop codons in a DNA sequence."""
    stop_codon_indexes = []
    stop_codons = ['tga', 'tag', 'taa']
    for i in range(frame, len(dna), 3):
        codon = dna[i: i + 3].lower()
        if codon in stop_codons:
            stop_codon_indexes.append(i)
    return stop_codon_indexes


def fasta_parser(file):
    """This function parses a .fasta file and returns a dictionary of
    headers:sequences pairs."""
    if file.endswith(".fasta"):
        file = file
    else:
        file = f"{file}.fasta"
    try:
        with open(file, 'r') as f:
            seqs = {}
            for line in f:
                line = line.rstrip()  # this deletes all newline (\n) characters
                if line.startswith('>'):  # check if this is a header character
                    words = line.split()  # splits the header line on whitespaces
                    # [0] to use the first word only as a key in the dictionary
                    name = words[0][1:]
                    # [1:] to remove the '>' sign used in fasta format from the
                    # first position of the header
                    # generates a new empty entry in the dictionary with 'name' as the key
                    seqs[name] = ''
                else:  # if the line is not a header
                    # remove whitespaces
                    seq = rmv_new_lines_spaces(line)
                    # appends the line on the equivalent dictionary entry
                    seqs[name] += seq
    except IOError:
        print("File not found!!")
        return

    return seqs


def split_by_n_letters(seq, word_size=10):
    """This function splits the sequence to smaller 'words' separated by spaces"""
    chunks = [seq[i: i + word_size] for i in range(0, len(seq), word_size)]
    return " ".join(chunks)


def split_to_lines_of_length_n(seq, n=50):
    """This function splits the DNA sequence in 50 nts on each line."""
    lines = []
    for i in range(0, len(seq), n):
        lines.append(seq[i: i + n])
    return lines


def save_fasta(headers, seqs, file_name):
    """Function to save a list sequences in a .fasta file."""
    if file_name.endswith(".fasta") or file_name.endswith("faa"):
        file_name = file_name
    else:
        file_name = f"{file_name}.fasta"
    with open(file_name, "w") as handle:
        for header, seq in zip(headers, seqs):
            handle.write(f"{header}\n")
            lines = split_to_lines_of_length_n(seq)
            final_seq = ""
            for line in lines:
                final_seq += f"{split_by_n_letters(line)}\n"
            handle.write(final_seq)
            handle.write("\n")
