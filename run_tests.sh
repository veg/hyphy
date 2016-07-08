echo "ALIGNMENTS";
`pwd`/HYPHYMP LIBPATH=`pwd`/res/ `pwd`/tests/hbltests/libv3/alignments.bf; cat errors.log;date

echo "MODELS";

echo "GENERIC CODON"
`pwd`/HYPHYMP LIBPATH=`pwd`/res/ `pwd`/tests/hbltests/libv3/models-codon.bf; cat errors.log;date

echo "HKY85";
`pwd`/HYPHYMP LIBPATH=`pwd`/res/ `pwd`/tests/hbltests/libv3/models/hky85.bf; cat errors.log;date

#echo "ABSREL";
#`pwd`/HYPHYMP LIBPATH=`pwd`/res/ `pwd`/tests/hbltests/libv3/models/absrel.bf; cat errors.log;date

echo "TREES";
`pwd`/HYPHYMP LIBPATH=`pwd`/res/ `pwd`/tests/hbltests/libv3/trees.bf; cat errors.log;date
