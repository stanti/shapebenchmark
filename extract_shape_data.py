#!/usr/bin/env python

import os
import sys
import xlrd

workbook = xlrd.open_workbook(sys.argv[1])
outpath = sys.argv[2]

#extract sequences
sheet = workbook.sheet_by_index(1)
sequences = {}

for r in range(1, sheet.nrows):
	row = sheet.row(r)

	name = row[0].value.strip()
	sequence = row[1].value.strip()
	sequences[name] = sequence

	with open(os.path.join(outpath, name + ".fa"), "w") as f:
		f.write(">" + name + "\n")
		f.write(sequence + "\n")


#extract SHAPE reactivities
sheet = workbook.sheet_by_index(2)
for c in range(0, sheet.ncols):
	column = sheet.col(c)

	name = column[0].value.strip()
	sequence = sequences[name]
	with open(os.path.join(outpath, name + ".shape"), "w") as f:
		for i in range(1,sheet.nrows):
			if column[i].ctype == 0:
				break

			assert(column[i].ctype == 2)
			v = column[i].value

			if v < -500: #values lower than -500 are treated as missing by RNAstructure
				continue
			if i-1 >= len(sequence): #unfortunately one sequence in the SHAPE data set has more reactivities than nucleotides :-/
				continue

			f.write("%d %s %g\n" % (i, sequences[name][i-1], v if v > 0 else 0)) #values between -500 and 0 will be mapped to 0 (same handling as RNAstructure)

