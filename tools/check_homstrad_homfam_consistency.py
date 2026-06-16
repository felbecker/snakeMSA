## Checks if the reference sequences in the homfam collection have exact
## matches to reference sequences in the homstrad collection.


import argparse
import os

from pathlib import Path

from learnMSA.util import SequenceDataset


argparser = argparse.ArgumentParser(
    description="Checks if the reference sequences in the homfam collection " \
    "have exact matches to reference sequences in the homstrad collection."
)
argparser.add_argument(
    "--homfam", type=str,
    default="data/homfam/aligned",
    help="Path to the homfam collection."
)
argparser.add_argument(
    "--homstrad", type=str,
    default="data/homstrad/aligned",
    help="Path to the homstrad collection."
)
args = argparser.parse_args()


if __name__ == '__main__':

    homfam_path = Path(args.homfam)
    homstrad_path = Path(args.homstrad)

    for homfam_file in homfam_path.glob("*.vie"):
        # Try to find the corresponding homstrad file by name matching
        homstrad_file = homstrad_path / (homfam_file.stem + ".fasta")
        if not homstrad_file.exists():
            print(f"Could not find homstrad file for {homfam_file.name}.")
            continue
        homfam_data = SequenceDataset(homfam_file)
        homstrad_data = SequenceDataset(homstrad_file)

        for sid in homfam_data.seq_ids:
            if sid not in homstrad_data.seq_ids:
                print(
                    f"Sequence ID {sid} in {homfam_file.name} not found in "\
                    f"{homstrad_file.name}."
                )
                continue
            homfam_i = homfam_data.index(sid)
            homstrad_i = homstrad_data.index(sid)
            homfam_s = homfam_data.get_standardized_seq(homfam_i)
            homstrad_s = homstrad_data.get_standardized_seq(homstrad_i)
            if homfam_s != homstrad_s:
                print(
                    f"Sequence ID {sid} in {homfam_file.name} does not match "
                    f"the corresponding sequence in {homstrad_file.name}:"
                    f"\n  homfam:  \t{homfam_s}"
                    f"\n  homstrad: \t{homstrad_s}"
                )