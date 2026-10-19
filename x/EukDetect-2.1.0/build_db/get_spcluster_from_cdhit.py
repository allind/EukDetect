#! /usr/bin/env python
import sys
def main(argv):
	
	toprint = []
	shouldprint = False
	for line in open(sys.argv[1]):
		line = line.strip('\n')
		if ">" in line[0]:
			if shouldprint:
				for e in toprint:
					print(e)
			shouldprint = False
			toprint = []
			toprint.append(line)
		else:
			if line[0] == "0":
				toprint.append(line)
			elif line[0] == "1":
				toprint.append(line)
				shouldprint = True
			else:
				toprint.append(line)
	if shouldprint:
		for e in toprint:
			print(e)
if __name__ == "__main__":
  main(sys.argv)
