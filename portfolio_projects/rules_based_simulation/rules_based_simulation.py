"""
Rules-Based Epigenome Simulation
==================================
Investigates how well a simple writer/eraser model explains epigenome-wide
patterns of successive histone methylation (H3K4, H3K9, H3K27, H3K36).

Model:
  - Nucleosomes sit on a 1D chromatin fibre
  - Each nucleosome carries a methylation state: {0 = unmodified, 1 = me1, 2 = me2, 3 = me3}
  - WRITER enzymes add a methyl group to a nucleosome with probability P_write
    if a neighbouring nucleosome is already methylated (spreading rule)
  - ERASER enzymes remove a methyl group with probability P_erase
  - Simulation runs for N Monte Carlo steps; statistics are recorded every epoch

Usage:
    python rules_based_simulation.py --nucleosomes 200 --steps 50000 --p_write 0.08 --p_erase 0.02
    python rules_based_simulation.py --sweep           # run parameter sweep
"""

import argparse
import random
import math
import csv
import os
from collections import defaultdict


# ---------------------------------------------------------------------------
# Chromatin fibre model
# ---------------------------------------------------------------------------
class Nucleosome:
    def __init__(self, position, state=0):
        self.position = position
        self.state    = state  # 0–3 (me0, me1, me2, me3)

    def __repr__(self):
        return f"Nuc({self.position}, me{self.state})"


class ChromatinFibre:
    """
    1D array of nucleosomes modelling a stretch of chromatin.
    Boundary conditions: periodic (circular chromosome).
    """
    def __init__(self, n_nucleosomes=100, seed_positions=None, seed_state=3):
        self.n      = n_nucleosomes
        self.nucs   = [Nucleosome(i) for i in range(n_nucleosomes)]
        # Seed initial methylation at specified positions
        if seed_positions:
            for pos in seed_positions:
                if 0 <= pos < n_nucleosomes:
                    self.nucs[pos].state = seed_state
        else:
            # Default: seed the centre nucleosome
            centre = n_nucleosomes // 2
            self.nucs[centre].state = seed_state

    def neighbours(self, idx):
        """Return left and right neighbour states (periodic boundary)."""
        left  = self.nucs[(idx - 1) % self.n].state
        right = self.nucs[(idx + 1) % self.n].state
        return left, right

    def max_neighbour_state(self, idx):
        left, right = self.neighbours(idx)
        return max(left, right)

    def methylation_profile(self):
        return [nuc.state for nuc in self.nucs]

    def fraction_methylated(self):
        return sum(1 for nuc in self.nucs if nuc.state > 0) / self.n

    def mean_state(self):
        return sum(nuc.state for nuc in self.nucs) / self.n

    def state_counts(self):
        counts = defaultdict(int)
        for nuc in self.nucs:
            counts[nuc.state] += 1
        return dict(counts)


# ---------------------------------------------------------------------------
# Monte Carlo simulation engine
# ---------------------------------------------------------------------------
class WriterEraserModel:
    """
    Simple writer/eraser model for histone methylation spreading.

    Rules:
      WRITE: A nucleosome at me_k gains +1 methylation (→ me_{k+1}) with
             probability p_write × (neighbour_state / 3).
             This encodes cooperative spreading: higher neighbour state = more likely to spread.

      ERASE: A nucleosome at me_k loses −1 methylation (→ me_{k-1}) with
             probability p_erase (constitutive erasure).

      AUTONOMOUS WRITE: A nucleosome gains +1 without neighbour dependence with
             probability p_auto (basal writer activity, noise term).
    """

    def __init__(self, fibre, p_write=0.08, p_erase=0.02, p_auto=0.001, rng_seed=42):
        self.fibre   = fibre
        self.p_write = p_write
        self.p_erase = p_erase
        self.p_auto  = p_auto
        random.seed(rng_seed)

        self.step        = 0
        self.history     = []   # list of (step, mean_state, fraction_methylated, state_counts)

    def _attempt_write(self, idx):
        nuc            = self.fibre.nucs[idx]
        max_nb_state   = self.fibre.max_neighbour_state(idx)
        if nuc.state < 3 and max_nb_state > 0:
            p = self.p_write * (max_nb_state / 3.0)
            if random.random() < p:
                nuc.state += 1

    def _attempt_erase(self, idx):
        nuc = self.fibre.nucs[idx]
        if nuc.state > 0:
            if random.random() < self.p_erase:
                nuc.state -= 1

    def _attempt_auto_write(self, idx):
        nuc = self.fibre.nucs[idx]
        if nuc.state < 3:
            if random.random() < self.p_auto:
                nuc.state += 1

    def mc_step(self):
        """One Monte Carlo step: attempt one update per nucleosome (random order)."""
        indices = list(range(self.fibre.n))
        random.shuffle(indices)
        for idx in indices:
            self._attempt_write(idx)
            self._attempt_erase(idx)
            self._attempt_auto_write(idx)
        self.step += 1

    def run(self, n_steps, record_every=500):
        """Run N Monte Carlo steps, recording statistics every `record_every` steps."""
        print(f"  Running {n_steps:,} MC steps on {self.fibre.n} nucleosomes...")
        for _ in range(n_steps):
            self.mc_step()
            if self.step % record_every == 0:
                self._record()

        print(f"  Simulation complete. {len(self.history)} snapshots recorded.")

    def _record(self):
        counts = self.fibre.state_counts()
        self.history.append({
            'step':               self.step,
            'mean_state':         round(self.fibre.mean_state(), 4),
            'fraction_methylated': round(self.fibre.fraction_methylated(), 4),
            'me0': counts.get(0, 0),
            'me1': counts.get(1, 0),
            'me2': counts.get(2, 0),
            'me3': counts.get(3, 0),
        })


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------
def compute_spreading_length(profile):
    """
    Compute the methylation spreading length as the half-width at half-maximum
    of the methylation profile around the seed nucleosome.
    Returns the number of nucleosomes over which spreading is detectable.
    """
    if not any(s > 0 for s in profile):
        return 0
    max_state = max(profile)
    half_max  = max_state / 2.0
    above_half = [i for i, s in enumerate(profile) if s >= half_max]
    if not above_half:
        return 0
    return above_half[-1] - above_half[0] + 1


def analyse_final_state(fibre, model_params, n_steps):
    """Print final state analysis."""
    profile   = fibre.methylation_profile()
    counts    = fibre.state_counts()
    spreading = compute_spreading_length(profile)
    frac_met  = fibre.fraction_methylated()
    mean_st   = fibre.mean_state()

    print(f"\n{'='*55}")
    print(f"Final State Analysis")
    print(f"{'='*55}")
    print(f"  Nucleosomes:        {fibre.n}")
    print(f"  MC steps:           {n_steps:,}")
    print(f"  P_write:            {model_params['p_write']}")
    print(f"  P_erase:            {model_params['p_erase']}")
    print(f"  P_auto:             {model_params['p_auto']}")
    print(f"\n  Final nucleosome states:")
    for state, count in sorted(counts.items()):
        bar = '#' * (count * 40 // fibre.n)
        print(f"    me{state}: {count:>4} ({count/fibre.n*100:5.1f}%)  {bar}")
    print(f"\n  Fraction methylated (me1+me2+me3): {frac_met:.3f}")
    print(f"  Mean methylation state:            {mean_st:.3f}")
    print(f"  Spreading length (HWHM):           {spreading} nucleosomes")

    # Pattern classification
    if frac_met > 0.8:
        pattern = "HYPERMETHYLATED (heterochromatin-like)"
    elif frac_met < 0.1:
        pattern = "HYPOMETHYLATED (euchromatin-like)"
    elif spreading < fibre.n // 5:
        pattern = "FOCAL domain (sharp boundary)"
    else:
        pattern = "BROAD spreading domain"
    print(f"  Pattern classification:            {pattern}")
    print(f"{'='*55}")


def save_timeseries(history, output_path):
    """Save simulation timeseries to CSV."""
    if not history:
        return
    fieldnames = list(history[0].keys())
    with open(output_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(history)
    print(f"\n  Timeseries saved to: {output_path}")


def save_final_profile(fibre, output_path):
    """Save the final methylation profile to CSV."""
    with open(output_path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['position', 'methylation_state'])
        for nuc in fibre.nucs:
            writer.writerow([nuc.position, nuc.state])
    print(f"  Final profile saved to: {output_path}")


# ---------------------------------------------------------------------------
# Parameter sweep
# ---------------------------------------------------------------------------
def parameter_sweep(n_nucleosomes=100, n_steps=20000):
    """
    Sweep over p_write and p_erase values to explore the phase space.
    Reports fraction methylated and spreading length at steady state.
    """
    p_write_values = [0.02, 0.05, 0.08, 0.12, 0.20]
    p_erase_values = [0.01, 0.02, 0.05, 0.10, 0.20]

    print(f"\n{'='*65}")
    print(f"Parameter Sweep: p_write vs p_erase  ({n_nucleosomes} nucleosomes, {n_steps} steps)")
    print(f"{'='*65}")
    print(f"{'p_write':>10} {'p_erase':>10} {'frac_met':>10} {'mean_state':>12} {'spreading':>12}")
    print(f"{'-'*10} {'-'*10} {'-'*10} {'-'*12} {'-'*12}")

    results = []
    for pw in p_write_values:
        for pe in p_erase_values:
            fibre = ChromatinFibre(n_nucleosomes)
            model = WriterEraserModel(fibre, p_write=pw, p_erase=pe, p_auto=0.001, rng_seed=42)
            # Suppress per-step output for sweep
            for _ in range(n_steps):
                model.mc_step()

            profile   = fibre.methylation_profile()
            frac_met  = fibre.fraction_methylated()
            mean_st   = fibre.mean_state()
            spreading = compute_spreading_length(profile)

            print(f"{pw:>10.3f} {pe:>10.3f} {frac_met:>10.3f} {mean_st:>12.3f} {spreading:>12}")
            results.append({
                'p_write': pw, 'p_erase': pe,
                'frac_methylated': round(frac_met, 4),
                'mean_state': round(mean_st, 4),
                'spreading_length': spreading,
            })

    # Save sweep results
    with open('parameter_sweep_results.csv', 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    print(f"\n  Sweep results saved to: parameter_sweep_results.csv")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description='Rules-Based Histone Methylation Simulation (Writer/Eraser Model)'
    )
    parser.add_argument('--nucleosomes', type=int,   default=200,    help='Number of nucleosomes')
    parser.add_argument('--steps',       type=int,   default=50000,  help='Monte Carlo steps')
    parser.add_argument('--p_write',     type=float, default=0.08,   help='Writer probability')
    parser.add_argument('--p_erase',     type=float, default=0.02,   help='Eraser probability')
    parser.add_argument('--p_auto',      type=float, default=0.001,  help='Autonomous write probability')
    parser.add_argument('--seed',        type=int,   default=42,     help='Random seed')
    parser.add_argument('--record',      type=int,   default=500,    help='Record every N steps')
    parser.add_argument('--sweep',       action='store_true',        help='Run parameter sweep')
    parser.add_argument('--output_dir',  default='.',                help='Output directory')
    args = parser.parse_args()

    if args.sweep:
        parameter_sweep()
        return

    os.makedirs(args.output_dir, exist_ok=True)

    # Build model
    fibre = ChromatinFibre(n_nucleosomes=args.nucleosomes)
    model = WriterEraserModel(
        fibre,
        p_write=args.p_write,
        p_erase=args.p_erase,
        p_auto=args.p_auto,
        rng_seed=args.seed,
    )

    params = {
        'p_write': args.p_write,
        'p_erase': args.p_erase,
        'p_auto':  args.p_auto,
    }

    print(f"Writer/Eraser Simulation")
    print(f"  Nucleosomes:  {args.nucleosomes}")
    print(f"  MC steps:     {args.steps:,}")
    print(f"  P_write:      {args.p_write}")
    print(f"  P_erase:      {args.p_erase}")
    print(f"  P_auto:       {args.p_auto}")

    model.run(args.steps, record_every=args.record)

    analyse_final_state(fibre, params, args.steps)

    save_timeseries(
        model.history,
        os.path.join(args.output_dir, 'simulation_timeseries.csv')
    )
    save_final_profile(
        fibre,
        os.path.join(args.output_dir, 'final_methylation_profile.csv')
    )


if __name__ == '__main__':
    main()
