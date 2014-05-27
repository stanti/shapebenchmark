#!/usr/bin/env bash

set -e

rnafold=~/ViennaRNA/src/bin/RNAfold
rnapvmin=~/ViennaRNA/src/bin/RNApvmin

#optional data extraction step (preprocessing of shapeknots data)
if [[ "$*" == *--extract* ]]
then
	mkdir -p benchmarkdata

	cp 3rdparty/shapeknots/sequences_ct_files/*.ct benchmarkdata/

	#rename some files (sequences in the xls document sometimes have slightly different names)
	mv "benchmarkdata/5domain16S rRNA, E. coli.ct"       "benchmarkdata/5' domain of 16S rRNA, E. coli.ct"
	mv "benchmarkdata/5domain16S rRNA, H. volcanii.ct"   "benchmarkdata/5' domain of 16S rRNA, H. volcanii.ct"
	mv "benchmarkdata/5domain23S rRNA, E. coli.ct"       "benchmarkdata/5' domain of 23S rRNA, E. coli.ct"
	mv "benchmarkdata/Group I intron, Azoarcus sp.ct"    "benchmarkdata/Group I intron, Azoarcus sp..ct"
	mv "benchmarkdata/Lysine_riboswitch, T. martima.ct"  "benchmarkdata/Lysine riboswitch, T. maritime.ct"
	mv "benchmarkdata/Mbox_riboswitch, B. subtilis.ct"   "benchmarkdata/M-Box riboswitch, B. subtilis.ct"
	mv "benchmarkdata/PreQ1_riboswitch, B. subtilis.ct"  "benchmarkdata/Pre-Q1 riboswitch, B. subtilis.ct"

	./extract_shape_data.py 3rdparty/shapeknots/ShapeKnots_SNRNASM.xlsx benchmarkdata

	#remove last reactivity value since the number of reactivity exceeds the sequence length
	head -n -1 "benchmarkdata/5' domain of 16S rRNA, H. volcanii.shape" > temp.txt ; mv temp.txt "benchmarkdata/5' domain of 16S rRNA, H. volcanii.shape"
fi


#predict mfe structure and pairing probability matrices for all sequences using RNAfold and various SHAPE methods
if [[ "$*" != *--skipfolding* ]]
then
	mkdir -p predictions
	mkdir -p runtime

	for rna in benchmarkdata/*.fa
	do
		name=`echo $rna | sed 's/.[^.]*$//' | sed 's|benchmarkdata/||'`

		echo "$name - RNAfold"
		/usr/bin/time -f %e -o "runtime/$name.R.time" $rnafold -p --bppmThreshold=1e-15 --MEA < "benchmarkdata/$name.fa" > "predictions/$name.R.out"
		mv *_dp.ps "predictions/$name.R.ps"

		echo "$name - Deigan"
		/usr/bin/time -f %e -o "runtime/$name.D.time" $rnafold -p --bppmThreshold=1e-15 --MEA "--shape=benchmarkdata/$name.shape" --shapeMethod=D < "benchmarkdata/$name.fa" > "predictions/$name.D.out"
		mv *_dp.ps "predictions/$name.D.ps"

		echo "$name - Zarringhalam"
		/usr/bin/time -f %e -o "runtime/$name.Z.time" $rnafold -p --bppmThreshold=1e-15 --MEA "--shape=benchmarkdata/$name.shape" --shapeMethod=Z < "benchmarkdata/$name.fa" > "predictions/$name.Z.out"
		mv *_dp.ps "predictions/$name.Z.ps"

		echo "$name - Washietl"
		/usr/bin/time -f %e -o "runtime/$name.P.time" $rnapvmin "benchmarkdata/$name.shape" < "benchmarkdata/$name.fa" > "predictions/$name.W.pv"
		/usr/bin/time -f %e -o "runtime/$name.W.time" $rnafold -p --bppmThreshold=1e-15 --MEA "--shape=predictions/$name.W.pv" --shapeMethod=W < "benchmarkdata/$name.fa" > "predictions/$name.W.out"
		mv *_dp.ps "predictions/$name.W.ps"

		if [[ "$*" == *--foldingperformance* ]]
		then
			for n in 0 1 2 3 4 5 6 7 8 9
			do
				/usr/bin/time -f %e -o "runtime/$name.R.time$n" $rnafold -p --bppmThreshold=1e-15 --MEA < "benchmarkdata/$name.fa" > /dev/null
				/usr/bin/time -f %e -o "runtime/$name.D.time$n" $rnafold -p --bppmThreshold=1e-15 --MEA "--shape=benchmarkdata/$name.shape" --shapeMethod=D < "benchmarkdata/$name.fa" > /dev/null
				/usr/bin/time -f %e -o "runtime/$name.Z.time$n" $rnafold -p --bppmThreshold=1e-15 --MEA "--shape=benchmarkdata/$name.shape" --shapeMethod=Z < "benchmarkdata/$name.fa" > /dev/null
				/usr/bin/time -f %e -o "runtime/$name.W.time$n" $rnafold -p --bppmThreshold=1e-15 --MEA "--shape=predictions/$name.W.pv" --shapeMethod=W < "benchmarkdata/$name.fa" > /dev/null
			done

			/usr/bin/time -f %e -o "runtime/$name.S.time" $rnapvmin "benchmarkdata/$name.shape" --sampleSize=10000 < "benchmarkdata/$name.fa" > /dev/null
		fi

		rm *.ps
	done
fi


#rate predictions and create plots
mkdir -p results

./compare_sequences.py predictions benchmarkdata results

for n in mfesens mfeppv measens meappv prob ensemblediv structurediv
do
	awk -F'\t' -f transpose.awk < results/$n.csv > results/t$n.csv
	./plot.R results/t$n.csv results/$n.svg results/${n}_diff.svg $n
done


