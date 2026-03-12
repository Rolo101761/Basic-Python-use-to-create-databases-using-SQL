import argparse
import csv
from collections import defaultdict
def readData(file, lower, upper):
    """
    Read data file, keep only peptides within the specified m/z range
    Returns dictionary: {pepid: mass}
    pepid is composed as protein_pepnum_pepseq to ensure uniqueness
    """
    peps = {}
    with open(file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 6:
                continue
            protein = cols[0]
            pepnum = cols[1]
            mass   = float(cols[2])
            pepseq = cols[5]
            if mass < lower or mass > upper:
                continue
            pepid = f"{protein}_{pepnum}_{pepseq}"
            peps[pepid] = mass
    return peps
# Mode 1: Count peptides within range
def mode1(pepmass):
    print(f"Total peptides in range: {len(pepmass)}")
# Mode 2: Fixed bin histogram
def mode2(pepmass, start, end, binsize, outfile):
    numbins = int((end - start) / binsize)
    results = []
    for i in range(numbins):
        lower = start + i * binsize
        upper = lower + binsize
        count = sum(1 for mz in pepmass.values() if lower <= mz < upper)
        results.append([lower, upper, count])
    # Output CSV
    with open(outfile, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["bin_start", "bin_end", "count"])
        writer.writerows(results)
    print(f"Mode2 histogram saved to {outfile}")
# Mode 3: Sliding window
def mode3(pepmass, start, end, window, step, outfile):
    windows = []
    cur = start
    while cur + window <= end:
        lower = cur
        upper = cur + window
        count = sum(1 for mz in pepmass.values() if lower <= mz < upper)
        windows.append([lower, upper, count])
        cur += step
    with open(outfile, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["window_start", "window_end", "count"])
        w.writerows(windows)
    print(f"Mode3 sliding window saved to {outfile}")
# Mode 4: Count uniquely identified proteins
def mode4(pepmass, binsize):
    bin_map = defaultdict(list)
    # Assign peptides to "tolerance bins" based on binsize
    for pepid, mz in pepmass.items():
        bin_index = round(mz / binsize)
        bin_map[bin_index].append(pepid)
    # Find peptides with unique m/z values
    unique_peptides = set()
    for bin_id, peps in bin_map.items():
        if len(peps) == 1:
            unique_peptides.add(peps[0])
    # Extract protein names from unique peptides
    proteins = set()
    for pepid in unique_peptides:
        proteins.add(pepid)
    print(f"Proteins uniquely identified: {len(proteins)}")
# Main program entry
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Peptide analysis tool")
    parser.add_argument("fileName", help="input data file")
    parser.add_argument("--mode", type=int, required=True,
                        help="1=range count, 2=histogram, 3=sliding window, 4=unique proteins")
    parser.add_argument("-s", "--start", default=1000.0, type=float)
    parser.add_argument("-e", "--end",   default=1500.0, type=float)
    parser.add_argument("-b", "--binsize", default=1.0, type=float)
    parser.add_argument("--window", type=float, default=10.0)
    parser.add_argument("--step",   type=float, default=1.0)
    parser.add_argument("-o", "--output", default="output.tsv")
    args = parser.parse_args()
    # Read data
    pepmass = readData(args.fileName, args.start, args.end)
    # Execute different modes
    if args.mode == 1:
        mode1(pepmass)
    elif args.mode == 2:
        mode2(pepmass, args.start, args.end, args.binsize, args.output)
    elif args.mode == 3:
        mode3(pepmass, args.start, args.end, args.window, args.step, args.output)
    elif args.mode == 4:
        mode4(pepmass, args.binsize)
    else:
        print("Unknown mode. Use --mode 1/2/3/4.")
