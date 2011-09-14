/* test preamble */

	_testDescription 		= " fit the non-reversible nucleotide model with parametric base frequencies to an alignment of 13 HIV env V3 sequences";
	_expectedLL 			= -1139.17487295371;
	ExecuteAFile 			("../Shared/TestInstrumentation.bf");
	startTestTimer 			(_testDescription);

/* end test preamble */

DataSet 		ds = ReadDataFile("../data/HIVenvSweden.seq");
DataSetFilter  filteredData = CreateFilter (ds,1);
ExecuteAFile   ( HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR +  "TemplateModels" + DIRECTORY_SEPARATOR + "NRM+Freqs.mdl");
Tree NREVTree = DATAFILE_TREE;

LikelihoodFunction lf = ( filteredData, NREVTree );	
Optimize ( res, lf );

/* test epilogue */
	timeMatrix = endTestTimer 				  (_testDescription);
	if (logTestResult (Abs (res[1][0] - _expectedLL) < 0.01))
	{
		return timeMatrix;
	}
	return 0;
/* end test epilogue */