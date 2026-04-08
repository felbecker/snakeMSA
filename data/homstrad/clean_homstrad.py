#!/usr/bin/env python3
"""
Based on:
https://github.com/steineggerlab/foldmason-analysis/blob/main/homstrad/clean_homstrad.py

Scans raw HOMSTRAD entries, keeps only those with at least --min-sequences
sequences, and for each qualifying family:
  - Splits *-sup.pdb into individual PDB files under pdb/<family>/
  - Writes an aligned FASTA file under aligned/<family>.fasta
  - Writes a raw (gap-free) FASTA file under unaligned/<family>.fasta

Usage:
    python3 clean_homstrad.py homstrad_with_PDB/ aligned/ pdbs/ unaligned/ [--min-sequences 5]
"""

import argparse
import re
from pathlib import Path
from collections import defaultdict


# Where Homstrad PDB files do not match sequence in Homstrad MSA
# 1ka5a (hpr)    Remove final K
# 1pi2  (bowman) SR -> AA end of sequence
# 1m85a (cat)    Remove final E
replacements = {
    "1ka5a": "MEQNSYVIIDETGIHARPATMLVQTASKFDSDIQLEYNGKKVNLKSIMGVMSLGVGKDAEITIYADGSDESDAIQAISDVLSKEGLT--",
    "1pi2": "---YSKPCCDLCMCTRSMPPQCSCED-RINSCHSDCKSCMCTRSQPGQCRCLDTNDFCYKPCKAA--------",
    "1m85a": "---------------------------------------------------KKLTTAAGAPVVDNNNVITAGPRGPMLLQDVWFLEKLAHFDREVIPERR-HAKGSGAFGTFTVTHDITKYTRAKIFSEVGKKTEMFARFSTVAGERGAADAERDIRGFALKFYTEEGNWDMVGNNTPVFYLRDPLKFPDLNHIVKRDPRTNM----RNMAYKWDFFSHLPESLHQLTIDMSDRGLPLSYRFVHGFGSHTYSFINKDNERFWVKFHFRCQQGIKNLMDDEAEALVGKDRESSQRDLFEAIERGDYPRWKLQIQIMPEKEASTVPYNPFDLTKVWPHADYPLMDVGYFELNRNPDNYFSDVEQAAFSPANIVPGISFSPDKMLQGRLFSYGDAHRYRL-GVNHHQIPVNAPK-CPFHNYHRDGAMRVDGNSGNGITYEPNS--GGVFQEQPDF------KEPPLSIEGAADHWNHREDEDYFSQPRALYE-LLSDDEHQRMFARIAGELSQA-SKETQQRQIDLFTKVHPEYGAGVEKAIKVL--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
}

amino_acid_map = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
    'B': 'ASX', 'Z': 'GLX', 'X': 'XAA', '-': '---'
}


def parse_ali(ali_file):
    """Parse a .ali file and return (members, chain_names, sequences)."""
    with ali_file.open() as fp:
        text = fp.read()
    members = []
    sequences = []
    chain_names = []
    for name, chain, sequence in re.findall(
        r"^>.+?;(?P<name>.+?)$\n^.+?:.+?:.+?:(?P<chain>.+?):.+?$\n(?P<sequence>.+?)\*",
        text,
        re.MULTILINE | re.DOTALL
    ):
        members.append(name)
        chain_names.append(chain if chain.isalnum() else 'A')
        if name in replacements:
            sequences.append(replacements[name])
        else:
            sequences.append(str(sequence).replace('\n', '').replace('/', '-'))
    return members, chain_names, sequences


def process_family(folder, aligned_dir, pdb_dir, unaligned_dir):
    """Process a single HOMSTRAD family folder into the three output directories."""
    folder = Path(folder)
    alignment = folder / f"{folder.stem}.ali"
    superposition = folder / f"{folder.stem}-sup.pdb"

    output_pdbs = Path(pdb_dir) / folder.name
    output_pdbs.mkdir(parents=True, exist_ok=True)

    members, chain_names, sequences = parse_ali(alignment)

    # Read superposition PDB
    with superposition.open() as fp:
        text = fp.read()

    # Find chain-to-name mapping from REMARK section, e.g.:
    # REMARK   1bif   chain A
    remark = {
        chain: name
        for (name, chain) in
        re.findall(r"^\s*?REMARK\s*(?P<name>\w+?)\s*chain\s*(?P<chain>[A-Z0-9]*)\s*?$", text, re.MULTILINE)
    }

    # Split superposed PDB into per-member files
    single_pdbs = defaultdict(str)
    last_chain = None
    last_residue = None
    chain_index = 0
    chain_ = chain_names[0]
    residue_index = 1
    name = ""
    for line in text.split('\n'):
        if line.startswith("ATOM"):
            chain = line[21:22]
            residue = line[22:27]
            if residue != last_residue:
                residue_index += 1
                last_residue = residue
            if chain != last_chain:
                name = remark.get(chain, members[chain_index])
                chain_ = chain_names[chain_index]
                chain_index += 1
                residue_index = 1
                last_chain = chain
            line = line[0:21] + chain_ + f" {residue_index:-3}" + line[26:] + '\n'
            single_pdbs[name] += line

    # Prepend SEQRES records
    for name, chain, sequence in zip(members, chain_names, sequences):
        three = [amino_acid_map.get(r, 'UNK') for r in sequence.replace('-', '')]
        lines = [three[i:i + 13] for i in range(0, len(three), 13)]
        seqres = [
            f"SEQRES {idx:>3}{chain:>2} {len(three):>4}  " + " ".join(residues)
            for idx, residues in enumerate(lines, start=1)
        ]
        single_pdbs[name] = "\n".join(seqres) + "\n" + single_pdbs[name]

    # Write individual PDB files
    for member in members:
        assert member in single_pdbs, f"{member} not parsed from {folder}"
        pdb_path = output_pdbs / f"{member}.pdb"
        with pdb_path.open('w') as fp:
            fp.write(single_pdbs[member] + '\n')

    # Write aligned MSA as FASTA
    Path(aligned_dir).mkdir(parents=True, exist_ok=True)
    msa_output = Path(aligned_dir) / f"{folder.name}.fasta"
    with msa_output.open('w') as fp:
        fp.write("".join(f">{name}\n{seq}\n" for name, seq in zip(members, sequences)))

    # Write gap-free sequence FASTA
    Path(unaligned_dir).mkdir(parents=True, exist_ok=True)
    aa_output = Path(unaligned_dir) / f"{folder.name}.fasta"
    with aa_output.open('w') as fp:
        fp.write("".join(f">{name}\n{seq.replace('-', '')}\n" for name, seq in zip(members, sequences)))


def main(homstrad_dir, aligned_dir, pdb_dir, unaligned_dir, min_sequences):
    homstrad = Path(homstrad_dir)
    assert homstrad.is_dir(), f"{homstrad} is not a directory"

    families = sorted(d for d in homstrad.iterdir() if d.is_dir())
    kept = 0
    skipped_missing = 0
    skipped_small = 0

    for family in families:
        ali = family / f"{family.stem}.ali"
        sup = family / f"{family.stem}-sup.pdb"

        if not ali.exists() or not sup.exists():
            print(f"Skipping {family.name}: missing .ali or -sup.pdb")
            skipped_missing += 1
            continue

        members, _, _ = parse_ali(ali)
        n = len(members)
        if n < min_sequences:
            print(f"Skipping {family.name}: {n} sequence(s) (< {min_sequences})")
            skipped_small += 1
            continue

        print(f"Processing {family.name}: {n} sequences")
        try:
            process_family(family, aligned_dir, pdb_dir, unaligned_dir)
            kept += 1
        except Exception as e:
            print(f"  ERROR processing {family.name}: {e}")
            skipped_missing += 1

    print(f"\nDone: {kept} families processed, "
          f"{skipped_small} skipped (< {min_sequences} sequences), "
          f"{skipped_missing} skipped (missing files / errors)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter and process HOMSTRAD families with at least N sequences"
    )
    parser.add_argument(
        "homstrad_dir",
        help="Directory containing raw HOMSTRAD family folders"
    )
    parser.add_argument(
        "--aligned_dir",
        type=str,
        help="Output folder for aligned FASTA files",
        default="aligned"
    )
    parser.add_argument(
        "--pdb_dir",
        type=str,
        help="Output folder for per-family PDB subfolders",
        default="pdbs"
    )
    parser.add_argument(
        "--unaligned_dir",
        type=str,
        help="Output folder for unaligned FASTA files",
        default="unaligned"
    )
    parser.add_argument(
        "--min-sequences",
        type=int,
        default=5,
        help="Minimum number of sequences to keep a family (default: 5)"
    )
    args = parser.parse_args()
    main(args.homstrad_dir, args.aligned_dir, args.pdb_dir, args.unaligned_dir, args.min_sequences)
