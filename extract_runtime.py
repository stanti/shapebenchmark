#!/usr/bin/env python

import os
import sys

runtimefolder = sys.argv[1]
reffolder = sys.argv[2]

data = {}

for file in os.listdir(runtimefolder):
	path = os.path.join(runtimefolder, file)

	components = file.split('.')
	method = components[-2]
	name = '.'.join(components[:-2])

	if not name in data:
		data[name] = {}
	if not method in data[name]:
		data[name][method] = []

	with open(path, "r") as f:
		data[name][method].append(float(f.readline().strip()))

print '\t'.join(['Name', 'Length', 'RNAfold', 'Deigan', 'Zarringhalam', 'Washietl', 'Washietl perturbation vector calculation'])

for name,namedata in data.iteritems():
	out = [name]

	with open(os.path.join(reffolder, name + ".fa"), "r") as f:
		f.readline()
		out.append(str(len(f.readline().strip())))

	for m in ['R', 'D', 'Z', 'W', 'P']:
		avgtime = reduce(lambda x, y: x + y, namedata[m]) / len(namedata[m])
		out.append(str(avgtime))

	print '\t'.join(out)

