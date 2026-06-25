#!/usr/bin/env python3
"""
sync_matched.py

Copy files / pdb-subfolders from `source` into `target`, using `reference`
to decide what counts as a "match".

Expected layout for both `reference` and `source`:
    <path>/aligned     (files)
    <path>/unaligned   (files)
    <path>/pdb         (subfolders)

The reference/aligned subdirectory is treated as the single source of truth
for which name stems are valid (e.g. "1abc" counts as a real entity if and
only if it shows up in reference/aligned) - reference/unaligned and
reference/pdb are not consulted for matching at all.

On the source side, aligned/unaligned/pdb are assumed to share the same set
of name stems, so matching and reporting are done once per stem, not once
per subfolder - a stem is either matched against the reference/aligned truth
list (and copied wherever it appears in source: aligned, unaligned, and/or
pdb) or not, and is reported by its stem name only (never a full filepath).

  - Stems present in source but not in reference/aligned are skipped
    (reported under "Skipped - no match in reference").
  - Stems present in reference/aligned but not in source are left untouched
    (reported under "Untouched in reference - no match in source").

By default, a "stem" is the filename without extension (for pdb subfolders,
which normally have no extension, stem == folder name). Pass --exact to
require the full name including extension to match instead.

Usage:
    python sync_matched.py REFERENCE SOURCE TARGET [--exact] [--dry-run]
"""

import argparse
import shutil
import sys
from pathlib import Path

# subfolder name -> kind of entries it holds ("file" or "dir")
SUBFOLDERS = {
    "aligned": "file",
    "unaligned": "file",
    "pdb": "dir",
}


def match_key(entry: Path, exact: bool) -> str:
    """Key (stem) used to decide whether entries 'match' across roots."""
    if exact:
        return entry.name
    return entry.stem  # for directories with no '.', stem == name anyway


def index_entries(folder: Path, kind: str, exact: bool):
    """Index the direct children of `folder` that match `kind` ('file' or 'dir').

    Returns:
        index:      dict[stem -> Path]
        wrong_kind: list[str] stems of entries that exist but aren't the expected kind
        collisions: list[str] stems that appeared more than once (first one wins)
    """
    index = {}
    wrong_kind = []
    collisions = []

    if not folder.is_dir():
        return index, wrong_kind, collisions

    for entry in sorted(folder.iterdir()):
        if entry.name.startswith("."):
            continue  # ignore hidden files like .DS_Store
        is_right_kind = (kind == "file" and entry.is_file()) or (kind == "dir" and entry.is_dir())
        key = match_key(entry, exact)
        if not is_right_kind:
            wrong_kind.append(key)
            continue
        if key in index:
            collisions.append(key)
            continue
        index[key] = entry

    return index, wrong_kind, collisions


def index_root(root: Path, exact: bool):
    """Index all three subfolders of `root`.

    Returns:
        per_folder: dict[subfolder_name -> dict[stem -> Path]]
        all_stems:  set of every stem seen across the three subfolders
        wrong_kind: set of stems that showed up as the wrong kind somewhere
        collisions: set of stems that collided somewhere
    """
    per_folder = {}
    wrong_kind = set()
    collisions = set()

    for name, kind in SUBFOLDERS.items():
        idx, wk, col = index_entries(root / name, kind, exact)
        per_folder[name] = idx
        wrong_kind.update(wk)
        collisions.update(col)

    all_stems = set()
    for idx in per_folder.values():
        all_stems.update(idx.keys())

    return per_folder, all_stems, wrong_kind, collisions


def copy_entry(entry: Path, dest_folder: Path, kind: str, dry_run: bool) -> None:
    if dry_run:
        return
    dest_folder.mkdir(parents=True, exist_ok=True)
    dest_path = dest_folder / entry.name
    if kind == "file":
        shutil.copy2(entry, dest_path)
    else:
        if dest_path.exists():
            shutil.rmtree(dest_path)
        shutil.copytree(entry, dest_path)


def copy_matched_stem(stem: str, src_per_folder, target_root: Path, dry_run: bool) -> None:
    """Copy whichever subfolders (aligned/unaligned/pdb) contain this stem in source."""
    for name, kind in SUBFOLDERS.items():
        entry = src_per_folder[name].get(stem)
        if entry is not None:
            copy_entry(entry, target_root / name, kind, dry_run)


def fmt_list(stems):
    return ", ".join(sorted(stems)) if stems else "(none)"


def print_report(matched, source_unmatched, reference_untouched,
                  ref_wrong_kind, src_wrong_kind, ref_collisions, src_collisions,
                  dry_run, exact):
    prefix = "[DRY RUN] " if dry_run else ""
    match_mode = "exact filename" if exact else "filename stem (extension ignored)"
    print(f"Matching mode: {match_mode}\n")

    print(f"{prefix}Matched and copied ({len(matched)}): {fmt_list(matched)}")
    print(f"Skipped - no match in reference ({len(source_unmatched)}): {fmt_list(source_unmatched)}")
    print(f"Untouched in reference - no match in source ({len(reference_untouched)}): {fmt_list(reference_untouched)}")

    if src_wrong_kind:
        print(f"Ignored in source - unexpected type for that subfolder ({len(src_wrong_kind)}): {fmt_list(src_wrong_kind)}")
    if ref_wrong_kind:
        print(f"Ignored in reference - unexpected type for that subfolder ({len(ref_wrong_kind)}): {fmt_list(ref_wrong_kind)}")
    if src_collisions:
        print(f"Warning - duplicate stems within a source subfolder, first one wins ({len(src_collisions)}): {fmt_list(src_collisions)}")
    if ref_collisions:
        print(f"Warning - duplicate stems within a reference subfolder, first one wins ({len(ref_collisions)}): {fmt_list(ref_collisions)}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("reference", type=Path, help="Path to the reference root folder")
    parser.add_argument("source", type=Path, help="Path to the source root folder")
    parser.add_argument("target", type=Path, help="Path to the target root folder (created if missing)")
    parser.add_argument(
        "--exact", action="store_true",
        help="Match on full filename incl. extension instead of stem only",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Only report what would happen; copy nothing",
    )
    args = parser.parse_args()

    reference_root = args.reference.resolve()
    source_root = args.source.resolve()
    target_root = args.target.resolve()

    for label, root in (("reference", reference_root), ("source", source_root)):
        if not root.is_dir():
            sys.exit(f"Error: {label} path does not exist or is not a directory: {root}")

    if target_root.exists() and any(target_root.iterdir()):
        print(f"Warning: target directory {target_root} is not empty - continuing anyway.\n")

    # Reference truth list comes from reference/aligned only.
    ref_index, ref_wrong_kind, ref_collisions = index_entries(reference_root / "aligned", "file", args.exact)
    ref_stems = set(ref_index.keys())

    # Source side still scans all three subfolders, since that's what gets copied.
    src_per_folder, src_stems, src_wrong_kind, src_collisions = index_root(source_root, args.exact)

    matched = sorted(ref_stems & src_stems)
    source_unmatched = sorted(src_stems - ref_stems)
    reference_untouched = sorted(ref_stems - src_stems)

    for stem in matched:
        copy_matched_stem(stem, src_per_folder, target_root, args.dry_run)

    print_report(
        matched, source_unmatched, reference_untouched,
        ref_wrong_kind, src_wrong_kind, ref_collisions, src_collisions,
        args.dry_run, args.exact,
    )


if __name__ == "__main__":
    main()