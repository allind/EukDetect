#! /usr/bin/env python
#usage: script.py -fp fwd_prefix -rp rev_prefix -outdir [dir] -fd fastq_dir -t taxonomy_dir -o [config dest file]
#this will estimate your readlen

import sys,subprocess, shlex, os, argparse
def main(argv):
	
	parser = argparse.ArgumentParser(description = "Create a eukdetect config file")
	parser.add_argument("-outdir", type=str, action="store", dest="outdir", default=os.getcwd(), help="output directory")
	#will need to add a conditional about whether its pe or se
	parser.add_argument("-fs", type=str, action="store", dest="fs", help="suffix of forward files")
	parser.add_argument("-rs", type=str, action="store", dest="rs", help="suffix of reverse files")
	parser.add_argument("-fd", type=str, action="store", dest="fqdir", default=os.getcwd(), help="directory of fastq files")
	parser.add_argument("-t", type=str, action="store", dest="taxdir",  help="taxdir")
	parser.add_argument("-o", type=str, action="store", dest="outdest", default="config.yml", help="output config file")
	parser.add_argument("-readmin", type=str, action="store", dest="readmin", help="minimum read length")
	args = parser.parse_args()
	#note: will need to standardize names of things in taxonomy directory


	#print(args.outdir, args.fs, args.rs, args.taxdir, args.fqdir)
	#iterate over the fq dir and get sample names
	samples = []
	if os.path.isdir(args.fqdir):
		for f in os.listdir(args.fqdir):
			if args.fs in f:
				sn = f.replace(args.fs, "")
				#check if rs exist
				if not os.path.isfile(args.fqdir + "/" + sn + args.rs):
					print("no reverse reads for " + sn)
					sys.exit()
				samples.append(sn)
	else:
		print('fastq dir provided does not exist')
		sys.exit()
	
	dest = open(args.outdest, 'w')
	dest.write('rootdir: "' + args.outdir.rstrip('/') + '"\n')
	dest.write("paired_end: true\n")
	dest.write('fwd_suffix: "' + args.fs + '"\n')
	dest.write('rev_suffix: "' + args.rs + '"\n')
	dest.write('fq_dir: "' + args.fqdir.rstrip('/') + '"\n')
	dest.write('taxonomy_dir: "' + args.taxdir.rstrip('/') + '"\n')
	dest.write('readmin: ' + args.readmin + '\n')
	dest.write('\n')
	dest.write("samples:\n")
	for sn in samples:
		dest.write("  " + sn + ":\n")
	dest.close()

if __name__ == "__main__":
  main(sys.argv)
