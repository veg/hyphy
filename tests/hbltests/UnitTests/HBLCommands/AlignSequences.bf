ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "AlignSequences";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Syntax:  AlignSequences (result, input_string_matrix,  options_matrix);
  // Parameters:         result: result is an empty container which will hold the sequence alignments.
  //                     input_string_matrix: A vector string matrix that contains the sequences to be aligned.
  //                     options_matrix: A series of alignment options to be configured.
  
  // Some stock parameters for the options_matrix from `res/TemplateBatchFiles/Utility/HXB2Mapper.bf`
  _hxb_alignOptions_nuc = {};
  _hxb_alignOptions_nuc ["SEQ_ALIGN_CHARACTER_MAP"]="ACGT";

  // Using a 4x4 matrix throws the following error: "The dimension of the scoring matrix (SEQ_ALIGN_SCORE_MATRIX) was not the expected dimension: 5x5"
  //_hxb_alignOptions_nuc ["SEQ_ALIGN_SCORE_MATRIX"] = 	{{5,-4,-4,-4}{-4,5,-4,-4}{-4,-4,5,-4}{-4,-4,-4,5}};
  _hxb_alignOptions_nuc ["SEQ_ALIGN_SCORE_MATRIX"] = 	{{5,-4,-4,-4,-4}{-4,5,-4,-4,-4}{-4,-4,5,-4,-4}{-4,-4,-4,5,-4}{-4,-4,-4,-4,5}};

  _hxb_alignOptions_nuc ["SEQ_ALIGN_GAP_OPEN"]	= 	50;
  _hxb_alignOptions_nuc ["SEQ_ALIGN_GAP_OPEN2"]	= 	50;
  _hxb_alignOptions_nuc ["SEQ_ALIGN_GAP_EXTEND"]	= 	1;
  _hxb_alignOptions_nuc ["SEQ_ALIGN_GAP_EXTEND2"]	= 	1;
  _hxb_alignOptions_nuc ["SEQ_ALIGN_AFFINE"]		=   1;
  _hxb_alignOptions_nuc ["SEQ_ALIGN_NO_TP"]		=   1;

  simpleSeqMatrix1 = {{"ACACACCCTTTTACACACACAC"}{"ACACATTTTAGAGAGAGAG"}};
  simpleSeqMatrix2 = {{"ACACACCCTTTTACACACACAC","ACACATTTTAGAGAGAGAG"}};

  
  // TODO: AlignSequences no longer seg faults but it just returns and empty array.
  AlignSequences(alignedSimple1, simpleSeqMatrix1, _hxb_alignOptions_nuc);
  AlignSequences(alignedSimple2, simpleSeqMatrix2, _hxb_alignOptions_nuc);

  fprintf (stdout, 'alignedSimple1: ', alignedSimple1, '\n');
  fprintf (stdout, 'alignedSimple2: ', alignedSimple2, '\n');





  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  list1 = {'key1':'val1', 'key2':'val2'};
  matrix1 = {{0.3,0.4}};
  Topology T1 = ((1,4),(3,4),5);
  Tree TT1 = ((1,2),(3,4),5);
  
  assert(runCommandWithSoftErrors("AlignSequences(placeholder, list1, _hxb_alignOptions_nuc);", "did not evaluate to a Matrix in call to AlignSequences"), "Failed error checking for trying to align sequences on list");
  assert(runCommandWithSoftErrors("AlignSequences(placeholder, matrix1, _hxb_alignOptions_nuc);", "did not evaluate to a dense string vector with"), "Failed error checking for trying to align sequences on matrix of numbers");
  assert(runCommandWithSoftErrors("AlignSequences(placeholder, T1, _hxb_alignOptions_nuc);", "did not evaluate to a Matrix in call to AlignSequences"), "Failed error checking for trying to align sequences on Topology");
  assert(runCommandWithSoftErrors("AlignSequences(placeholder, TT1, _hxb_alignOptions_nuc);", "did not evaluate to a Matrix in call to AlignSequences"), "Failed error checking for trying to align sequences on Tree");
  
  testResult = 1;

  return testResult;
}
