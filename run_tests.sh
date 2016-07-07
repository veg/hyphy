echo "ALIGNMENTS";
/Users/sweaver/Programming/hyphy/HYPHYMP LIBPATH=/Users/sweaver/Programming/hyphy/res/ /Users/sweaver/Programming/hyphy/tests/hbltests/libv3/alignments.bf; cat errors.log;date

echo "MODELS";

echo "GENERIC CODON"
/Users/sweaver/Programming/hyphy/HYPHYMP LIBPATH=/Users/sweaver/Programming/hyphy/res/ /Users/sweaver/Programming/hyphy/tests/hbltests/libv3/models-codon.bf; cat errors.log;date

echo "HKY85";
/Users/sweaver/Programming/hyphy/HYPHYMP LIBPATH=/Users/sweaver/Programming/hyphy/res/ /Users/sweaver/Programming/hyphy/tests/hbltests/libv3/models/hky85.bf; cat errors.log;date

echo "ABSREL";
/Users/sweaver/Programming/hyphy/HYPHYMP LIBPATH=/Users/sweaver/Programming/hyphy/res/ /Users/sweaver/Programming/hyphy/tests/hbltests/libv3/models/absrel.bf; cat errors.log;date

echo "TREES";
/Users/sweaver/Programming/hyphy/HYPHYMP LIBPATH=/Users/sweaver/Programming/hyphy/res/ /Users/sweaver/Programming/hyphy/tests/hbltests/libv3/trees.bf; cat errors.log;date
