ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName ()
{
	return "Branch Length";
}		

function getTestedFunctions ()
{
	return {{"_TreeTopology::BranchLength","_CalcNode::BranchLength","_Matrix::ExpNumberOfSubs"}};
}	

function runTest ()
{
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult  	   = 0;
	
	Topology T 			   = (a,b,c);
	assert (BranchLength (T,0) == (-1), "Retrieve a branch length by a valid index when no branch length has been defined");

	Topology T 			   = ((a:0.1,b:0.2):0.4,c:0.15,d:0.33);

	assert (BranchLength (T,0) == 0.1, "Retrieve a branch length by a valid index");
	assert (Abs(BranchLength (T,-1) - {{0.1,0.2,0.4,0.15,0.33,0.0}}) == 0.0, "Retrieve all branch lengths");
	assert (BranchLength (T,"b") == 0.2, "Retrieve a branch length by a valid name");
	assert (BranchLength (T,"a;c") == 0.65, "Retrieve a valid path length; the first node is deeper than the second in the tree");
	assert (BranchLength (T,"d;d") == 0.0, "Retrieve a valid trivial path length; checking boundary conditions");
	assert (BranchLength (T,"d;b") == BranchLength (T,1)+BranchLength(T,2)+BranchLength(T,4), "Retrieve a valid path length; the second node is deeper than the first in the tree");
	
	/* invalid parameters */
	
	assert (Type (BranchLength (T,5)) == "Unknown", "Trying to retrieve a branch length using an invalid index");
	assert (Type (BranchLength (T,"bb;d")) == "Unknown", "Trying to retrieve a path length using an 'invalid;valid' path specification");
	assert (Type (BranchLength (T,"b;ddd")) == "Unknown", "Trying to retrieve a path length using an 'valid;invalid' path specification");
	assert (Type (BranchLength (T,{{1,2}})) == "Unknown", "Trying to retrieve a branch length using an invalid argument type (matrix)");
	assert (Type (BranchLength (T,"nosuchnode")) == "Unknown","Trying to retrieve a branch length using an invalid node name");
	assert (Type(BranchLength (T,"a;EXPECTED_NUMBER_OF_SUBSTITUTIONS")) == "Unknown", "Trying to retrieve a branch length from a topology object");
	
	global 	  kappa = 4;
	
	Q_HKY85 = {{*,t,kappa*t,t}
			   {t,*,t,t*kappa}
			   {t*kappa,t,*,t}
			   {t,t*kappa,t,*}};
			   
	freqs	= {{0.4}{0.3}{0.2}{0.1}};
	
	Model HKY85 = (Q_HKY85, freqs, 1);
	
	/* now the default branch lengths are
	   measured in expected numbers of substitutions per site */
	   
	BRANCH_LENGTH_STENCIL = 0;
	LARGE_MATRIX_BRANCH_LENGTH_MODIFIER_DIMENSION = 21;
	
	Tree  T = ((a:0.1,b:0.2):0.4,c:0.15,d:0.33);
	
	assert 	   (BranchLength (T,0) == 0.136, "Retrieve a branch length by a valid index, when there is a valid (HKY85) model attached to the node");
	assert 	   (BranchLength (T,"a;EXPECTED_NUMBER_OF_SUBSTITUTIONS") == "0.4800000000000001*t+0.22*kappa*t", 
								  "Retrieve the branch length expression for a valid node name; no BRANCH_LENGTH_STENCIL");
	
	/* using BRANCH_LENGTH_STENCIL we can specify
	WHICH substitutions to count in the branch length;
	in the case below -- transitions only first and then 
	transversions only */
	
	BRANCH_LENGTH_STENCIL = {{0,0,1,0}
							 {0,0,0,1}
							 {1,0,0,0}
							 {0,1,0,0}};
	
	assert 	   (BranchLength (T,"a;EXPECTED_NUMBER_OF_SUBSTITUTIONS") == "0.22*kappa*t", "Retrieve the branch length expression for a valid node name subject to BRANCH_LENGTH_STENCIL");
	assert 	   (BranchLength (T,0) == 0.088, "Retrieve a branch length by a valid index, conditioned by BRANCH_LENGTH_STENCIL (transitions)");
	
	BRANCH_LENGTH_STENCIL = BRANCH_LENGTH_STENCIL ["_MATRIX_ELEMENT_VALUE_==0"];
			   
	assert 	   (BranchLength (T,0) == 0.048, "Retrieve a branch length by a valid index, conditioned by BRANCH_LENGTH_STENCIL (transversions)");
	
	/* 
	override model definition and retrieve only the 
	'original' branch length supplied during tree construction 
	*/
	
	BRANCH_LENGTH_STENCIL = "STRING_SUPPLIED_LENGTHS";
	assert 	   (BranchLength (T,0) == 0.1, "Retrieve a branch length by a valid index, conditioned by BRANCH_LENGTH_STENCIL = 'STRING_SUPPLIED_LENGTHS'");
	
	/* reset the branch length stencil */
	
	BRANCH_LENGTH_STENCIL = 0;
	assert 	   (BranchLength (T,0) == 0.136, "Testing BRANCH_LENGTH_STENCIL reset");
	
	LARGE_MATRIX_BRANCH_LENGTH_MODIFIER = 2;
	assert 	   (BranchLength (T,0) == 0.136, "Testing LARGE_MATRIX_BRANCH_LENGTH_MODIFIER; should not apply because of LARGE_MATRIX_BRANCH_LENGTH_MODIFIER_DIMENSION");
	LARGE_MATRIX_BRANCH_LENGTH_MODIFIER_DIMENSION = 2;
	assert 	   (BranchLength (T,0) == 0.136/2, "Testing LARGE_MATRIX_BRANCH_LENGTH_MODIFIER; should apply ");
	
	testResult = 1;
		
	return testResult;
}

/* execution stub */

