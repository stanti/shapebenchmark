#!/usr/bin/env python

import os
import re
import sys

resultdir = sys.argv[1]
refdir = sys.argv[2]
outdir = sys.argv[3]

def dotBracket2Pairs(string):
  pairs = {}
  stack = []
  for i in range(len(string)):
    c = string[i]
    if c == '(':
      stack.append(i)
    elif c == ')':
      pairs[stack.pop()+1] = i+1
    else:
      assert(c == '.')

  assert(len(stack) == 0)
  return pairs

def countCorrectPredictions(validated, predicted):
  ret = 0
  for (i, j) in predicted.iteritems():
    if i in validated and validated[i] == j:
      ret += 1

  return float(ret)

def calculateSensitivity(validated, predicted):
  return countCorrectPredictions(validated, predicted) / len(validated) if len(validated) else 0

def calculatePPV(validated, predicted):
  return countCorrectPredictions(validated, predicted) / len(predicted) if len(predicted) else 0


#resultfiles = os.listdir(refdir)
#names = sorted(set([x[0:-5] for x in resultfiles]))

#print resultfiles
names = [x[0:-3] for x in filter(lambda x: x.endswith('.fa'), os.listdir(refdir))]
allnames = {}
for n in names:
	allnames[n] = n

'''
allnames = {"PreQ1_riboswitch_B._subtilis": "Pre-Q1 riboswitch, B. subtilis *",
            "Telomerase_pseudoknot_human": "Telomerase pseudoknot, human *",
            "Fluoride_riboswitch_P.syringae": "Fluoride riboswitch, P. syringae *",
            "Adenine_riboswitch_V.vulnificus": "Adenine riboswitch, V. vulnificus",
            "tRNA_asp_yeast": "tRNA(aps), yeast",
            "tRNA_phe_E._coli": "tRNA(phe), E. coli",
            "TPP_riboswitch_E.coli": "TPP riboswitch, E. coli",
            "SARS_corona_virus_pseudoknot": "SARS corona virus pseudoknot *",
            "cyclic-di-GMP_riboswitch_V.cholerae": "cyclic-di-GMP riboswitch, V. cholerae",
            "SAM_I_riboswitch_T.tengcongensis": "SAM I riboswitch, T. tengcongensis *",
            "5S_rRNA_E.coli": "5S rRNA, E. coli",
            "Mbox_riboswitch_B._subtilis": "M-Box riboswitch, B. subtilis",
            "P546_domain_bI3_group_I_intron": "p546 domain, bI3 group I intron",
            "Lysine_riboswitch_T._martima": "Lysine riboswitch, T. martimia *",
            "Group_I_intron_Azoarcus_sp": "Group I intron, Azoarcus sp. *",
            "Signal_recognition_particle_RNA_human": "Signal recognition particle, human", #XX
            "Hepatitis_C_virus_IRES_domain": "Hepatitis C virus, IRES domain *",
            "RNase_P_B.subtilis": "RNase P, B. subtilis *", #XX
            "Group_II_intron_O.iheyensis": "Group II intron, O. iheyensis *",
            "Group_I_Intron_T.thermophila": "Group I intron, T. thermophilia *",
            "5domain16S_rRNA_H.volcanii": "5' domain of 23S rRNA, H. volcanii",
            "HIV-1_5prime_pseudoknot_domain": "HIV-1 5' pseudoknot domain *",
            "5domain23S_rRNA_E.coli": "5' domain of 23S rRNA, E. coli",
            "5domain16S_rRNA_E.coli": "5' domain of 16S rRNA, E. coli"}
'''

outMfeSens = []
outMfePpv = []
outMeaSens = []
outMeaPpv = []
outProb = []
outEnsembleDiv = []
outStructureDiv = []

for name,description in allnames.iteritems():
	sys.stderr.write(name+ "\n")

	reference = {}

	with open(os.path.join(refdir, name + '.ct')) as f:
		f.readline()
		for line in f:
			components = [x for x in line.strip().split() if x]
			if components:
				i = int(components[0])
				j = int(components[4])
				if j != 0 and not j in reference:
					reference[i] = j

	dataMfeSens = []
	dataMfePpv = []
	dataMeaSens = []
	dataMeaPpv = []
	dataProb = []
	dataEnsembleDiv = []
	dataStructureDiv = []

	for type in ['R', 'D', 'Z', 'W']:
		f = open(os.path.join(resultdir, name + '.' + type + '.out'))
		lines = f.readlines()
		f.close()

		mfe = lines[2].split(' ')[0]
		result = dotBracket2Pairs(mfe)
		dataMfeSens.append(calculateSensitivity(reference, result))
		dataMfePpv.append(calculatePPV(reference, result))

		assert('MEA' in lines[-2])
		mea = lines[-2].split(' ')[0]
		result = dotBracket2Pairs(mea)
		dataMeaSens.append(calculateSensitivity(reference, result))
		dataMeaPpv.append(calculatePPV(reference, result))

		length = len(mfe)
		dataEnsembleDiv.append(float(lines[-1].strip().split(' ')[-1]) / length)


		probs = {}
		f = open(os.path.join(resultdir, name + '.' + type + '.ps'))
		for line in f:
			if not line.endswith('ubox\n'):
				continue

			components = line.split(' ')
			if len(components) != 4:
				continue

			i = int(components[0])
			j = int(components[1])
			prob = float(components[2])

			if i in probs:
				probs[i][j] = prob
			else:
				probs[i] = {j: prob}
		f.close()

		values = []
		for (i, j) in reference.iteritems():
			if i in probs and j in probs[i]:
				values.append(probs[i][j])
			else:
				values.append(0)

		dataProb.append(sum(values)/len(values))

		structureDiversity = 0.
		for i in range(1, len(mfe) + 1):
			for j in range(i+1, len(mfe) + 1):
				prob = probs[i][j] if i in probs and j in probs[i] else 0
				if i in reference and reference[i] == j:
					structureDiversity += 1-prob
				else:
					structureDiversity += prob

		dataStructureDiv.append(structureDiversity/length)

	outMfeSens.append([description, length] + dataMfeSens)
	outMfePpv.append([description, length] + dataMfePpv)
	outMeaSens.append([description, length] + dataMeaSens)
	outMeaPpv.append([description, length] + dataMeaPpv)
	outProb.append([description, length] + dataProb)
	outEnsembleDiv.append([description, length] + dataEnsembleDiv)
	outStructureDiv.append([description, length] + dataStructureDiv)

out = {"mfesens": outMfeSens, "mfeppv": outMfePpv, "measens": outMeaSens, "meappv": outMeaPpv, "prob": outProb, "ensemblediv": outEnsembleDiv, "structurediv": outStructureDiv}

for name,data in out.iteritems():
	f = open(os.path.join(outdir, name + ".csv"), "w")
	for o in sorted(data, key = lambda x: x[1]):
		f.write("\t".join([str(x) for x in o]) + "\n")
	f.close()


