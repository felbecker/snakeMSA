#parse tabular files and fix timestamps
#find entries like 1 day, 7:38:15 and convert to 1:07:38:15

import sys
import re

def fix_timestamps(filename):
    with open(filename, 'r') as f:
        for line in f:
            #find timestamps like 1 day, 7:38:15
            #and convert to 1:07:38:15
            #keep the rest of line before and after
            match = re.search(r"(\d+) day, (\d+):(\d+):(\d+)", line)
            if match:
                day = match.group(1)
                hour = match.group(2)
                minute = match.group(3)
                second = match.group(4)
                print(line.replace(match.group(0), day + ':' + hour + ':' + minute + ':' + second), end="")
            else:
                print(line, end="")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python fix_timestamps.py <filename>")
        sys.exit(1)
    fix_timestamps(sys.argv[1])