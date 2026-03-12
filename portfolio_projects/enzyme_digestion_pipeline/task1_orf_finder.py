#script to input a fasta file of multiple sequences or a genome and output a fasta file
#containing multiple open reading frames (of a user defined minimum size) with the header line containing
#the orf name, the frame it's in, the length of the orf, and the start position in the sequence
#created by Hannah Ehrlich (14226716)
import re, sys, os
import argparse
#function to read in fasta file and create dictionary of sequences and list of names
#first reads all lines in the file into a list, with debugging statement if the file isn't found
#parses the lines into a dictionary containing the name as a key and a list of the names in case this
#is needed for downstream analysis
#debugging statement for if fasta file is not formatted correctly and reassurance statement informing user
#of how many sequences were found in the input file

def input_fasta(fasta_file):
    try:
        with open(fasta_file, "r") as file:
            lines = file.readlines()
        print(f"File '{fasta_file}' successfully read")
    except FileNotFoundError:
        print(f"Error: File '{fasta_file}' not found.")
        return [], {}

    name = []
    seq = {}
    current = None
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            current = line[1:].strip()
            name.append(current)
            seq[current] = ""
        else:
            if current is None:
                raise ValueError("FASTA format error: please check file before continuing.")
            seq[current] += line
    count = len(name)
    word = "sequence" if count == 1 else "sequences"
    print(f"{count} {word} found")
    return name, seq
#create reverse complementary sequence dictionary and use it to iterate over the sequence dictionary
#for every sequence found, reverse it and use the complementary dictionary to obtain the revcomp seq
def revcomplement(seq, name):
    compdict = {
        'A':'T','C':'G','G':'C','T':'A','N':'N',
        'a':'t','c':'g','g':'c','t':'a','n':'n'
    }
    revcomplements = {}
    for id in name:
        seq = seq[id]
        reversed_seq = seq[::-1]
        revcomplements[id] = ''.join(compdict[base] for base in reversed_seq)
    return revcomplements
#translate the reverse complement sequence to amino acids using dictionary of pre-defined amino acids and their codons
#iterates over the reverse complement sequences, stopping 2 bp before the end to avoid incomplete codons and using X for unknown codons
def translate(revcomplement, name):
    codons = { "AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N", "ACA": "T", "ACC": "T",
    "ACG": "T", "ACT": "T", "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
    "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I", "CAA": "Q", "CAC": "H",
    "CAG": "Q", "CAT": "H", "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R", "CTA": "L", "CTC": "L",
    "CTG": "L", "CTT": "L", "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A", "GGA": "G", "GGC": "G",
    "GGG": "G", "GGT": "G", "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "TAA": "*", "TAC": "Y", "TAG": "*", "TAT": "Y", "TCA": "S", "TCC": "S",
    "TCG": "S", "TCT": "S", "TGA": "*", "TGC": "C", "TGG": "W", "TGT": "C",
    "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F" }
    proteins = {}
    for id in name:
        dna = revcomplement[id]
        protein = ""
        for i in range(0, len(dna) - 2, 3):
            codon = dna[i:i+3].upper()
            protein += codons.get(codon, "X")  #"X" for unknown codons
        proteins[id] = protein
    return proteins
#naming function for orfs where the first 5 letters of the organism name is taken, along with the frame the orf is in and the index
def make_orf_name(org_name, frame, orf_index):
    org_code = org_name.replace(".", "").upper()[:5]
    return f"{org_code}_F{frame}_{orf_index:04d}"
#function to find orfs of a minimum size
#assigns * for stop codons and M for start codons
#for each frame, iterate over AA sequence, if theres a start codon start recording the orf, if there's a stop codon within current orf,
#add it to list
#makes dictionary of orfs with the frame as the key
def find_orfs(proteins, min_size=50, frame_offset=0):
    results = {}
    for org_name, seq in proteins.items():
        all_orfs = {}
        stop = "*"
        start = "M"
        orf_counter = {}
        for frame in range(1, 4):
            true_frame = frame + frame_offset
            orfs = []
            current_orf = ""
            in_orf = False
            start_pos = None
            orf_counter[true_frame] = 1
            for i in range(frame, len(seq)):
                aa = seq[i]
                if aa == start and not in_orf:
                    current_orf = "M"
                    in_orf = True
                    start_pos = i
                elif aa == stop and in_orf:
                    if len(current_orf) >= min_size:
                        orf_name = make_orf_name(org_name, true_frame, orf_counter[true_frame])
                        orfs.append({
                            "orf_name": orf_name,
                            "sequence": current_orf,
                            "frame": true_frame,
                            "length": len(current_orf),
                            "start": start_pos
                        })
                        orf_counter[true_frame] += 1
                    current_orf = ""
                    in_orf = False
                    start_pos = None
                elif in_orf:
                    current_orf += aa
            all_orfs[true_frame] = orfs
        results[org_name] = all_orfs
    return results
def merge_orfs(forward, reverse):
    merged = {}
    for org in forward:
        merged[org] = {}
        for frame in forward[org]:
            merged[org][frame] = forward[org][frame]
        for frame in reverse[org]:
            merged[org][frame] = reverse[org][frame]
    return merged
#gets rid of identical sequences across different frames to ensure only unique ORFs appear
#iterates through each name and frame
#uses the sequence as the unique identifier and if it has already been seen, it skips adding them to the output
#and increases the total_removed counter
def deduplicate_orfs(results):
    dedup_results = {}
    total_removed = 0
    for org_name, frames in results.items():
        seen = set()
        new_frames = {i: [] for i in frames}
        total_removed = 0
        for frame, orfs in frames.items():
            for orf in orfs:
                seq = orf["sequence"]
                if seq in seen:
                    total_removed += 1
                    continue
                seen.add(seq)
                new_frames[frame].append(orf)
        dedup_results[org_name] = new_frames
    print(f"Total duplicates removed across all sequences = {total_removed}")
    return dedup_results
#output to fasta file by iterating over each orf in dict
#format the fasta file with the header line containing the name, frame, length, and start position
def output_fasta(dedup_results, filename):
    output = []
    with open(filename, "w") as f:
        for org_name, frames in dedup_results.items():
            orf_count = 1
            for frame, orfs in frames.items():
                for orf in orfs:
                    header = f">{orf['orf_name']} frame={orf['frame']} length={orf['length']} start={orf['start']}"
                    f.write(header + "\n")
                    f.write(orf["sequence"] + "\n")
                    orf_count += 1
                    output.append(header + "\n" + orf["sequence"])
    return output
#command line arguments
def main():
    parser = argparse.ArgumentParser(
        description="Find ORFs in reverse-complemented DNA sequences from a FASTA file."
    )
    parser.add_argument("input_fasta", help="Input FASTA file with DNA sequences or genomes")
    parser.add_argument("output_fasta", help="Output FASTA file for ORFs")
    parser.add_argument(
        "--min_size", type=int, default=50,
        help="Minimum ORF length in amino acids (default: 50)"
    )
    args = parser.parse_args()
    #read fasta
    names, seqs = input_fasta(args.input_fasta)
    if not names:   #if no sequences were read, stop here
        print(f'No sequences found - nothing written to {args.output_fasta}.')
        sys.exit(1)
    #forward strand translated
    forward_proteins = translate(seqs, names)  # frames 1–3
    #reverse strand translated
    revcomps = revcomplement(seqs, names)
    reverse_proteins = translate(revcomps, names)  # frames 4–6
    #find ORFs for both strands
    forward_orfs = find_orfs(forward_proteins, min_size=args.min_size, frame_offset=0)
    reverse_orfs = find_orfs(reverse_proteins, min_size=args.min_size, frame_offset=3)
    #merge and deduplicate
    results = merge_orfs(forward_orfs, reverse_orfs)
    results = deduplicate_orfs(results)
    #output to fasta
    output_fasta(results, args.output_fasta)
    print(f'ORFs written to {args.output_fasta}')
if __name__ == "__main__":
    main()
