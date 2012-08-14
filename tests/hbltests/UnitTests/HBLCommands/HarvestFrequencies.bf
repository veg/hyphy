ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
	return "HarvestFrequencies";
}		

function getTestedFunctions ()
{
	return {{" _ElementaryCommand::HandleHarvestFrequencies","_DataSet::HarvestFrequencies", "_DataSet::constructFreq"}};
}	


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult  	   = 0;
    
    assert (runCommandWithSoftErrors ("HarvestFrequencies (corn/holio,filter,2)","Incorrect number of arguments"), "Failed error checking for an invalid syntax");
    assert (runCommandWithSoftErrors ("HarvestFrequencies (corn/holio,filter,1,1,1)","is not a valid variable identifier in call to HarvestFrequencies"), "Failed error checking for an invalid receptacle");
    assert (runCommandWithSoftErrors ("HarvestFrequencies (cornholio,filter,1,1,1)","is neither a DataSet nor a DataSetFilter"), "Failed error checking for an invalid dataset/dataset filter parameter");

    ExecuteAFile (PATH_TO_CURRENT_BF  + "res" + DIRECTORY_SEPARATOR + "test_likefunc.nex");

    HarvestFrequencies (nucFreqs, filteredData, 1, 1, 1);
    assert (Abs (nucFreqs-{
{    0.404450757576}
{    0.166287878788}
{    0.209564393939}
{    0.219696969697}
}) < 1e-8, "Checking nucleotide frequency counts");

    assert (runCommandWithSoftErrors ("HarvestFrequencies (cornholio,filteredData,3,2,1)","Atom should divide unit"), "Failed error checking for an invalid unit/atom specification");
    
    HarvestFrequencies (nucFreqsFromDS, ds, 1, 1, 1);
    assert (Abs (nucFreqs-nucFreqsFromDS) < 1e-8, "Checking nucleotide frequency counts collected from the DataSet object");
    
    HarvestFrequencies (nucFreqsFromDSExplicit, ds, 1, 1, 1,"0");
    assert (Abs (nucFreqsFromDSExplicit-{{    0}
{    1}
{    0}
{    0}
}) < 1e-8, "Checking nucleotide frequency counts collected from the DataSet object with a site partition string");

    HarvestFrequencies (nucFreqsFromDSExplicit, ds, 1, 1, 1,"0-4",speciesIndex==1); // CCCAT
    assert (Abs (nucFreqsFromDSExplicit-{{    .2}
{    .6}
{    0}
{    .2}
}) < 1e-8, "Checking nucleotide frequency counts collected from the DataSet object with a site/species partition");


    simpleDataSet = ">1\nARGT--\n>2\nAGGYCC";
    DataSet simpleTest = ReadFromString (simpleDataSet);
    
    COUNT_GAPS_IN_FREQUENCIES = 1;
    HarvestFrequencies (count1, simpleTest, 1, 1, 1);
    assert (Abs (count1-{{    3}
{    3}
{    4}
{    2} 
}* (1/12)) < 1e-8, "Checking simple frequency counts with - mapped to N-way ambiguities");
    
    COUNT_GAPS_IN_FREQUENCIES = 0;
    HarvestFrequencies (count2, simpleTest, 1, 1, 1);
    assert (Abs (count2-{{    2.5}
{    2.5}
{    3.5}
{    1.5} 
}* (1/10)) < 1e-8, "Checking simple frequency counts with - mapped to nil");
    
    HarvestFrequencies (count12_positional, simpleTest, 2, 1, 1);
    assert (Abs (count12_positional-{{    2,0.5}
{    1,1.5}
{    2,1.5}
{    0,1.5} 
}* (1/5)) < 1e-8, "Checking positional dinucleotide frequency counts");
    
    
    
    HarvestFrequencies (count12_overall, simpleTest, 2, 1, 0);
    assert (Abs (count2-count12_overall) < 1e-8, "Checking the equivalence of counting nucleotides directly or counting doublets while ignoring position");
    
    HarvestFrequencies (count22, simpleTest, 2, 2, 1);
    assert (Abs (count22-{16,1,{0,0,0.5},{2,0,1.5},{5,0,1},{9,0,.5},{11,0,1.5}}* (1/5)) < 1e-8, "Checking dinucleotide frequency counts");
    
    HarvestFrequencies (count3, simpleTest, 4, 1, 1); // check the case when the number of sites is not a multiple of the unit

	testResult = 1;
		
	return testResult;
}