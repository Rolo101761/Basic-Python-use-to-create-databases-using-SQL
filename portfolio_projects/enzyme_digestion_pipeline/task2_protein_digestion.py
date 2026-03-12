# alllows up to 1 missed cleavage as process isnt perfect. replicates real biological experiment.
# works with all 4 protease enzymes breaking long protien sequences for mass spectrophamaty. writes
# minimum peptide length filtering and writes peptide output to output file.
# amin miah
import argparse
import sys
# fasta reader
def fastaread(filename):
    order = []
    seqs = {}
    name = None
    seq = ""
  # trys opens the file gives error code if not found.
    try:
        f = open(filename, "r")
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)

    # read line by line
    with f:
        for line in f:
            line = line.strip()


            if line.startswith(">"):
                if name is not None:
                    order.append(name)
                    seqs[name] = seq
                name = line[1:].split()[0]
                seq = ""
            else:
                seq += line

    # after loop ends save final sequence
    if name is not None:
        order.append(name)
        seqs[name] = seq

    return order, seqs
# digestion function
def digest(sequence, enzyme, missed=1):
    # dictinary of cleavage rules for each enzyme tryspin cut after k or R unless aa= P
    rules = {
        "t":["K","R"],
        "lys-c":["K"],    # cuts only after k
        "arg-c":["R"],    # cuts only after r
        "glu-c":["E"]     # cuts after e
    }


    cut_sites = rules[enzyme]
    peptides = []
    final_out = []
    current = ""

    # primary cleavage loop

    for i in range(len(sequence)-1):
        aa = sequence[i]
        nxt = sequence[i+1]
        current += aa
        # checks if aa is in cleavage site and next aa is not p
        if aa in cut_sites and nxt != "P":
            peptides.append(current)
            current = ""  # starts new peptide

    # add final amino acid
    current += sequence[-1]
    peptides.append(current)


    # standard peptides

    for p in peptides:
        final_out.append((p,0))

    # missed cleavage max 1

    if missed == 1:
        for i in range(len(peptides) -1):
            combined = peptides[i] + peptides[i + 1]
            final_out.append((combined, 1))

    return final_out
# main function
if __name__=="__main__":

    parser = argparse.ArgumentParser(
    description="Task 2 - Protein digestion tool"

    )
    # positonal argument input fasta file
    parser.add_argument("fileName", help="Input FASTA file")

    #enzyme selection
    parser.add_argument(
    "-e", "--enzyme",
    choices=["t", "lys-c", "arg-c", "glu-c"],
    default="t",
    help="Enzyme: t (Tryspin), lys-c, arg-c, glu-c"

    )

    # missed cleavages (0 or 1)
    parser.add_argument(
    "-m", "--missed",
    type=int,
    default=1,
    help="Missed cleavages (0 or 1)"

    )

    # output file name
    parser.add_argument(
    "-o", "--output",
    default="digest.out",
    help="Output file name"

    )
    # minimum peptide length 5 is optimal.
    parser.add_argument(
    "-min", "--minLen",
    type=int,
    default=5,
    help="Minimum peptide length to output"
    )

    args = parser.parse_args()

    order, seqs = fastaread(args.fileName)

    # quality check
    print('%d sequences read in' % len(order))


    with open(args.output, "w") as out:
        for prot in order:
            peptides = digest(seqs[prot], args.enzyme, args.missed)

            count = 1
            for pep, miss in peptides:

                if len(pep) >= args.minLen:


                    out.write(f">{prot} peptide {count} missed={miss}\n{pep}\n")
                    count += 1


             # writes peptide which pass filter of length

        print("\nDigestion complete >", args.output)
