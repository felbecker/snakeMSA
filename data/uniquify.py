from sys import argv

## finds sequences with identical identifiers and makes their identifiers unique
with open(argv[1], "r") as file:
    lines = list(file.readlines())
    ids = {}
    for i,line in enumerate(lines):
        line = line.strip()
        if line[0] == ">":
            if line[1:] in ids:
                ids[line[1:]].append(i)
            else:
                ids[line[1:]] = [i]
for key, indices in ids.items():
    if len(indices) > 1:
        for k,i in enumerate(indices):
            lines[i] = lines[i].strip() + "_" + str(k) + "\n"
with open(argv[1], "w") as override_file:
    for line in lines:
        override_file.write(line)