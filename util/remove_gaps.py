import os
import sys

def remove_gaps(dir_name, outdir_name):
	for filename in os.listdir(dir_name):
		with open(dir_name+'/'+filename, 'r') as file:
			content = file.read().replace('-', '').replace('.', '')
			with open(outdir_name+'/'+filename, 'w') as out_file:
				out_file.write(content)

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("Usage: python remove_gaps.py <dir_name> <outdir_name>")
		sys.exit(1)
	remove_gaps(sys.argv[1], sys.argv[2])