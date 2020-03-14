ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "GetInformation";
}	




function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  
  // Simple example likelihood function for use later.
  DataSet         nucleotideSequences = ReadDataFile (PATH_TO_CURRENT_BF + "../../data/CD2_reduced.fna");
  DataSetFilter   filteredData = CreateFilter (nucleotideSequences,1);
  HarvestFrequencies (observedFreqs, filteredData, 1, 1, 1);
  F81RateMatrix = 
        {{*,mu,mu,mu}
         {mu,*,mu,mu}
         {mu,mu,*,mu}
         {mu,mu,mu,*}};
  Model   F81 = (F81RateMatrix, observedFreqs);
  Tree    givenTree = DATAFILE_TREE;
  LikelihoodFunction  LF = (filteredData, givenTree);
  Optimize (paramValues, LF);

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Returns a matrix whose contents vary upon the parameter value.

  // 	For a Number: The value of the variable, and the lower and upper bounds on the range as a 1x3 matrix(eg. \{\{var, lower, upper\}\})
  exampleNumber = 3221.174;
  GetInformation(exampleNumberInfo, exampleNumber);
  assert(exampleNumberInfo[0] == exampleNumber, "Failed to return the variable number as the first element in the info matrix");
  assert(exampleNumberInfo[1] == -1e+26, "Failed to return the expected lower bound of a float");
  assert(exampleNumberInfo[2] == 1e+26, "Failed to return the expected upper bound of a float");

  // For a Category: A two column matrix with a row per rate class, with column one containing the rates and column 2 containing the probabilities for the rates.
  category category1 = (4, EQUAL, MEAN, GammaDist(_x_,shapeParameter,shapeParameter), CGammaDist(_x_,shapeParameter,shapeParameter), 0 , 
			  							1e25,CGammaDist(_x_,shapeParameter+1,shapeParameter));
  GetInformation(category1Info, category1);
  assert(category1Info[1][0] == 0.25, "Failed to get information from a basic category variable");
  assert(category1Info[0][0] == 0, "Failed to get information from a basic category variable");

  // For a DataSetFilter: A string column, where each string is a sequence, ordered the same way they are in the filter.
  DataSet cd2nex = ReadDataFile (PATH_TO_CURRENT_BF + '../../data/CD2.nex');
  DataSetFilter onlyFiveThroughEleven = CreateFilter (cd2nex,1,"4-10");
  DataSetFilter onlyFiveThroughTen = CreateFilter (cd2nex,1,"4-9");
  GetInformation(onlyFiveThroughElevenInfo, onlyFiveThroughEleven);
  GetInformation(onlyFiveThroughTenInfo, onlyFiveThroughTen);
  assert((onlyFiveThroughTenInfo[0]+'C') == onlyFiveThroughElevenInfo[0], "Failed to output a matrix of sequence strings with GetInformation on a dataSetFilter");


  // TODO... GetInfo doesn't seem to be working for Trees, Likelihood Functions or Strings (regex)...
  
  //For a tree node: The rate matrix at that node (numeric). If the node has no associated rate matrix, the return value will be a 1x1 matrix whose value is not meaningful.
  /*
  Tree TT = ((a:0.1,b:0.2):0.4,c:0.15,d:0.33);
  fprintf (stdout, 'TT: ', TT, '\n');
  fprintf (stdout, 'TT.1: ', TT.2, '\n');
  GetInformation(node1Info, TT.2);
  fprintf (stdout, 'node1Info: ', node1Info, '\n');
  */

  // For a likelihood function: A string column vector containing the names of category variables the function depends on. The order is the same as that used for building marginal likelihood matrices.
  /*
  GetInformation(lfInfo, LF);
  fprintf (stdout, 'lfInfo: ', lfInfo, '\n');
  */

  //For a string: Finds all defined variable names which match the regular expression defined in the string.
  /*
  variableOne = 1;
  variableTwo = 2;
  variableThirtyTwo = 32;
  GetInformation(varabilesWithT, "Two");
  fprintf (stdout, 'variablesWithT: ', variablesWithT, '\n');
  */


  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  assert (runCommandWithSoftErrors ('GetInformation()', 'Incorrect number of arguments'), "Failed error checking for calling GetInformation without any arguments");
  assert (runCommandWithSoftErrors ('GetInformation(variable1)', 'Incorrect number of arguments'), "Failed error checking for calling GetInformation with only one argument");


  testResult = 1;

  return testResult;
}
