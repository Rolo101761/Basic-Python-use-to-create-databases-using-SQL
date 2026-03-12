import pandas as pd
import os, sys, argparse
# Amino acid masses dictionary, from here - https://www.matrixscience.com/help/aa_help.html
# Stores both the Monoisotopic (MI) and Average Mass (AM)
PROTON_MASS_CHARGE_NUM = 1.007
MASS_TO_CHARGE_NUMS = {
    "A":{"MI":71.037114,"AM":71.0779,},
    "R":{"MI":156.101111,"AM":156.1857,},
    "N":{"MI":114.042927,"AM":114.1026,},
    "D":{"MI":115.026943,"AM":115.0874,},
    "B":{"MI":0,"AM":0,},
    "C":{"MI":103.009185,"AM":103.1429,},
    "E":{"MI":129.042593,"AM":129.114,},
    "Q":{"MI":128.058578,"AM":128.1292,},
    "Z":{"MI":0,"AM":0,},
    "G":{"MI":57.021464,"AM":57.0513,},
    "H":{"MI":137.058912,"AM":137.1393,},
    "I":{"MI":113.084064,"AM":113.1576,},
    "L":{"MI":113.084064,"AM":113.1576,},
    "K":{"MI":128.094963,"AM":128.1723,},
    "M":{"MI":131.040485,"AM":131.1961,},
    "F":{"MI":147.068414,"AM":147.1739,},
    "P":{"MI":97.052764,"AM":97.1152,},
    "S":{"MI":87.032028,"AM":87.0773,},
    "T":{"MI":101.047679,"AM":101.1039,},
    "U":{"MI":150.95363,"AM":150.0379,},
    "W":{"MI":186.079313,"AM":186.2099,},
    "Y":{"MI":163.06332,"AM":163.1733,},
    "X":{"MI":0,"AM":0,},
    "V":{"MI":99.068414,"AM":99.1311,},
    "*":{"MI":0,"AM":0,},
    "WATER":{"MI":18.0106,"AM":18.0153,},
}
def read_peptide_fasta_file(file_name):
    if not os.path.exists(file_name):
        print(f"{file_name}: The file path does not exist.")
        sys.exit(0)

    def parse_peptide_seq(header, sequence):
        # get protein name, peptide number, and cleavage from the sequences header
        parts = header.split()
        protein_name = parts[0]
        peptide_num = int(parts[2])

        # handle cases where cleavage is in the format missed=X, X or is missing completely
        try:
            cleavage_num = parts[3]
            if "missed" in cleavage_num:
                cleavage_num = int(cleavage_num.split("=")[1])
            else:
                cleavage_num = int(cleavage_num)
        except (IndexError, ValueError):
            cleavage_num = 0

        # all values stored in a dictionary for each peptide sequence
        return {
            "Prot_name": protein_name,
            "peptide": peptide_num,
            "p": cleavage_num,
            "sequence": sequence,
        }

    peptide_data = []
    print(f"Reading FASTA file: {file_name}.")

    # opens the fasta file in read mode
    # header lines starting with '>' are stored
    # sequence lines are stored until the next header
    # when next header reached and we have a header and sequence pair,
    # get the sequence, protien name, peptide number and cleavage info
    with open(file_name, "r") as f:
        header = None
        sequence = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header and sequence:
                    pep = parse_peptide_seq(header, sequence)
                    peptide_data.append(pep)
                header = line[1:] # remove '>'
                sequence = None
            else:
                sequence = line
        # for final line/sequence in fasta file
        if header and sequence:
            pep = parse_peptide_seq(header, sequence)
            peptide_data.append(pep)

    print("File successfully read.")

    return pd.DataFrame(peptide_data)
def calculate_mass_to_charge(peptide_df, mass_type, charge):
    if charge < 1:
        print("Charge must be at least 1 for m/z calculation.")
        sys.exit(0)

    def get_seq_mass_to_charge(seq):
        # turn sequence into list of amino acids
        amino_acids = list(seq)
        # sum all the masses for each amino acid to get total mass of peptide
        total_mass = sum(MASS_TO_CHARGE_NUMS.get(aa, {}).get(mass_type, 0) for aa in amino_acids)
        # add water mass and proton charge mass
        total_mass += MASS_TO_CHARGE_NUMS.get("WATER", {}).get(mass_type, 0) + (charge * PROTON_MASS_CHARGE_NUM)
        # calculate final mass-to-charge value by dividing by charge
        m_z = total_mass / charge
        return m_z

    print(f"Calculating mass-to-charge for {len(peptide_df)} peptide sequences...")

    # for each sequence calculate its mass-to-charge (m/z), m/z = (peptide_mass + (z * proton_mass)) / z
    peptide_df["mass-to-charge"] = peptide_df["sequence"].apply(get_seq_mass_to_charge)

    print("All mass-to-charge values calculated successfully.")

    return peptide_df
def write_results(mass_results_df, file_name, charge):
    print(f"Writing mass-to-charge results to file: {file_name}.")

    # add charge column, and rearrange columns to fit specifications
    mass_results_df["z"] = charge
    mass_results_df = mass_results_df[["Prot_name", "peptide", "mass-to-charge", "z", "p", "sequence"]]
    # round all mass-to-charge values to 4 decimal places and then save
    mass_results_df["mass-to-charge"] = mass_results_df["mass-to-charge"].apply(lambda x: f"{x:.4f}")
    mass_results_df.to_csv(file_name, sep="\t", index=False, header=False)
    print("Mass-to-charge results successfully saved.")
def runner():
    task3_argparser = argparse.ArgumentParser(
        description='Calculates mass-to-charge values for protein peptide sequences in a FASTA formated file. '
    )
    task3_argparser.add_argument(
        "--mode",
        choices=["m", "a"],
        dest="mode",
        required=True,
        help="Mass type for mass-to-charge calculation: monoisotopic (m) or average masses (a).",
    )

    # only optional argument, defaults to 1
    task3_argparser.add_argument(
        "--charge",
        type=int,
        default=1,
        help="Charge state of the ion (default: +1)."
    )

    task3_argparser.add_argument(
        "-i",
        "--input_file",
        dest="input_file",
        required=True,
        help="Input FASTA file location containing the peptide sequence information.",
    )

    task3_argparser.add_argument(
        "-o",
        "--output_file",
        dest="output_file",
        required=True,
        help="Output file location for the mass-to-charge data.",
    )

    args = task3_argparser.parse_args()
    peptide_df = pd.DataFrame()
    mass_results_df = pd.DataFrame()

    # store all the parsed arguments
    mode = args.mode
    input_file = args.input_file
    output_file = args.output_file
    charge = args.charge

    # for selected mode get mass type key (m=MI=Monoisotopic, a=AM=Average Mass)
    mass_type = "MI" if mode == "m" else "AM"

    # read peptide sequences into a dataframe -> calculate masses for each peptide sequence -> save results
    peptide_df = read_peptide_fasta_file(input_file)
    mass_results_df = calculate_mass_to_charge(peptide_df, mass_type, charge)
    write_results(mass_results_df, output_file, charge)


if __name__ == "__main__":
    runner()
