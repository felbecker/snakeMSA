import os
import sys


to_replace = "UOBZJ"
replace_with = "X"

def replace(text):
    for aa in to_replace:
        text = text.replace(aa, replace_with)
    for aa in to_replace.lower():
        text = text.replace(aa, replace_with.lower())
    return text

def replace_in_file(filename):
    with open(filename, 'r') as file:
        for line in file:
            if not line.startswith(">"):
                print(replace(line.strip()))
            else:
                print(line.strip())

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print("Usage: python remove_gaps.py <file_name>")
		sys.exit(1)
	replace_in_file(sys.argv[1])