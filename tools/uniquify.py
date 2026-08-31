import argparse
import os


def find_duplicates(path):
    """Returns the lines of the file and a dict mapping duplicate identifiers to their line indices."""
    with open(path, "r") as file:
        lines = list(file.readlines())
    ids = {}
    for i,line in enumerate(lines):
        line = line.strip()
        if line == "":
            continue
        if line[0] == ">":
            if line[1:] in ids:
                ids[line[1:]].append(i)
            else:
                ids[line[1:]] = [i]
    duplicates = {key : indices for key, indices in ids.items() if len(indices) > 1}
    return lines, duplicates


## finds sequences with identical identifiers and makes their identifiers unique
def uniquify(path, dry=False):
    lines, duplicates = find_duplicates(path)
    if dry:
        if duplicates:
            print(path + ": " + str(len(duplicates)) + " duplicate identifier(s)")
            for key, indices in duplicates.items():
                print("  " + key + " appears " + str(len(indices)) + " times.")
        return bool(duplicates)
    print("Uniquifying identifiers in " + path + "...")
    for key, indices in duplicates.items():
        print("Identifier " + key + " appears " + str(len(indices)) + " times.")
        for k,i in enumerate(indices):
            lines[i] = lines[i].strip() + "_" + str(k) + "\n"
    with open(path + ".uniquified", "w") as override_file:
        for line in lines:
            override_file.write(line)
    return bool(duplicates)


DEFAULT_EXTENSIONS = [".fa", ".fasta", ".fas", ".faa", ".fna"]


def collect_files(path, extensions):
    if os.path.isdir(path):
        return sorted(os.path.join(path, f) for f in os.listdir(path)
                      if os.path.isfile(os.path.join(path, f))
                      and any(f.endswith(ext) for ext in extensions))
    return [path]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Makes fasta identifiers unique by appending a running index.")
    parser.add_argument("path", help="a fasta file or a directory of fasta files")
    parser.add_argument("-d", "--dry", action="store_true",
                        help="only list files with duplicate identifiers, write nothing")
    parser.add_argument("-e", "--extensions", nargs="+", default=DEFAULT_EXTENSIONS,
                        help="file extensions to consider when path is a directory (default: %(default)s)")
    args = parser.parse_args()

    files = collect_files(args.path, args.extensions)
    affected = [f for f in files if uniquify(f, dry=args.dry)]
    if args.dry:
        print(str(len(affected)) + " of " + str(len(files)) + " file(s) contain duplicate identifiers.")
