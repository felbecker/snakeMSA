import argparse

# a dummy aligner that does nothing but writing a trivial alignment to the 
# outputs

def main():
    parser = argparse.ArgumentParser(
        description="Dummy aligner that writes a trivial alignment to the output."
    )
    parser.add_argument(
        "-i", 
        "--input", 
        type=str, 
        required=True, 
        help="Input file."
    )
    parser.add_argument(
        "-o", 
        "--output", 
        type=str, 
        required=True, 
        help="Output file to write the trivial alignment."
    )
    
    args = parser.parse_args()

    with open(args.input, "r") as infile:
        lines = infile.readlines()

    # fill sequences with gaps up the length of the longest sequence
    max_length = max(
        len(line.strip()) 
        for line in lines if not line.startswith(">")
    )

    with open(args.output, "w") as outfile:
        for line in lines:
            line = line.strip()
            if line[0] != ">" and len(line) < max_length:
                line += "-" * (max_length - len(line))
            outfile.write(line + "\n")


if __name__ == "__main__":  
    main()