#!/usr/bin/env python

import os
import sys
import xlrd

workbook = xlrd.open_workbook(sys.argv[1])
outpath = sys.argv[2]

#extract sequences
sheet = workbook.sheet_by_index(1)

for r in range(1, sheet.nrows):
	row = sheet.row(r)

	with open(os.path.join(outpath, row[0].value.strip() + ".fa"), "w") as f:
		f.write(">" + row[0].value + "\n")
		f.write(row[1].value + "\n")


#extract SHAPE reactivities
sheet = workbook.sheet_by_index(2)
for c in range(0, sheet.ncols):
	column = sheet.col(c)

	with open(os.path.join(outpath, column[0].value.strip() + ".shape"), "w") as f:
		for i in range(1,sheet.nrows):
			if column[i].ctype == 0:
				break

			assert(column[i].ctype == 2)
			f.write("%d %f\n" % (i, column[i].value))
