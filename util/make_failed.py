import os
import sys

filepath = sys.argv[1]
tool = sys.argv[2]
with open(filepath, 'r') as file:
	for line in file:
		if len(line.strip()) == 0:
			continue
		split = line.strip().split("/")
		family = split[-1]
		dataset = split[1]
		outfile = "../outputs/"+tool+"/alignments/"+dataset+"/"+family
		if os.path.exists(outfile):
			print(outfile + " exists")
		else:
			with open(outfile, 'w') as out_file:
				out_file.write("FAILED")
			with open("../outputs/"+tool+"/benchmarks/"+dataset+"/"+family+".txt", 'w') as log_file:
				log_file.write("""s	h:m:s	max_rss	max_vms	max_uss	max_pss	io_in	io_out	mean_load	cpu_time\n
259200	72:0:0	0.  0.  0.  0.  0.  0.  0.  0. \n
""")