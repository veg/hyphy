ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();

function getTestName () {
	return "ReplicateConstraint";
}

function getTestedFunctions () {
	return {{"_ExecutionList::ReplicateConstraint"}};
}

lfunction checkTree (tree_id, callback) {
	tree_avl = (^tree_id)^0;
	
	result = TRUE;
	
	for (i = 1; i < Abs (tree_avl) - 1; i+=1) {
		result = result && Call (callback, tree_id + "." + (tree_avl[i])["Name"]);
	}
	
	return result;
}

lfunction check_syn_rate (name) {
	GetString (res, ^(name+".synRate"),-1);
	if (Type (res) == "AssociativeList") {
		if (Columns (res["Local"]) == 1) {
			return TRUE;
		}
	}
	
	return FALSE;
	
}

lfunction check_syn_rate0 (name) {
	GetString (res, ^(name+".synRate"),-1);
	if (Type (res) == "AssociativeList") {
		if (Columns (res["Local"]) == 0) {
			return TRUE;
		}
	}
	
	return FALSE;
	
}

lfunction check_syn_rate2 (name) {
	GetString (res, ^(name+".synRate"),-1);
	if (Type (res) == "AssociativeList") {
		if (Columns (res["Local"]) == 2) {
			return TRUE;
		}
	}
	
	return FALSE;
	
}

lfunction check_both (name) {
	GetString (res, ^(name+".nonSynRate"),-1);
	if (Type (res) == "AssociativeList") {
		if (Columns (res["Local"]) == 1) {
			return check_syn_rate (name);
		}
	}
	
	return FALSE;
	
}

lfunction check_both0 (name) {
	GetString (res, ^(name+".nonSynRate"),-1);
	if (Type (res) == "AssociativeList") {
		if (Columns (res["Local"]) == 0) {
			return check_syn_rate0 (name);
		}
	}
	
	return FALSE;
	
}
	

function runTest () {

	/** load test data **/
	
	ExecuteAFile ("res/replicate_constraint.nex");



	Tree second_tree =((D_CD_83_ELI_ACC_K03454,D_UG_94_94UG114_ACC_U88824)Node2,D_CD_84_84ZR085_ACC_U88822,(B_US_83_RF_ACC_M17451,((B_FR_83_HXB2_ACC_K03455,B_US_86_JRFL_ACC_U63632)Node11,B_US_90_WEAU160_ACC_U21135)Node10)Node8);
	Tree third_tree  = givenTree;
	ACCEPT_ROOTED_TREES = 1;
	Tree fourth_tree = (B_US_83_RF_ACC_M17451,((B_FR_83_HXB2_ACC_K03455,B_US_86_JRFL_ACC_U63632)Node11,B_US_90_WEAU160_ACC_U21135)Node10);
	
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult  	   = 0;


	//-----------------------------------------------------------------------------------------------------------------
	// ERROR HANDLING
	//-----------------------------------------------------------------------------------------------------------------

    assert (runCommandWithSoftErrors (
    	"ReplicateConstraint ('/*', givenTree)", 
    "The template for ReplicateConstraint must be an assignment"), "Failed error checking for an invalid receptacle");
    
    assert (runCommandWithSoftErrors (
    	"ReplicateConstraint ('this1[0][2]=this1[0][1]', givenTree)", 
    "The template for ReplicateConstraint must be an assignment"), "Failed error checking for an invalid receptacle");

    assert (runCommandWithSoftErrors (
    	"ReplicateConstraint ('this1.?.synRate', givenTree, givenTree)", 
    "'givenTree' \\(argument 2\\) did not appear in the contstraint expression"), "Failed error checking for an invalid receptacle");

     assert (runCommandWithSoftErrors (
    	"ReplicateConstraint ('this1.?.synRate := this2.?.nonSynRate', givenTree)", 
    "'this2.\\?.nonSynRate' does not have a matched positional argument"), "Failed error checking for an invalid receptacle");

     assert (runCommandWithSoftErrors (
    	"ReplicateConstraint ('this1.?.synRate := this2.?.nonSynRate', givenTree, second_tree)", 
    "'second_tree' \(argument 2\) is topologically incompatible with the reference argument"), "Failed error checking for an invalid receptacle");

    
 	//-----------------------------------------------------------------------------------------------------------------
	// SIMPLE CONSTRAINTS
	//-----------------------------------------------------------------------------------------------------------------


	global dNdS = 1;
	ReplicateConstraint ('this1.?.synRate := dNdS * this2.?.nonSynRate', givenTree, third_tree);
	// simple set of constraints; check to see that 
	//fprintf (stdout, LAST_SET_OF_CONSTRAINTS, "\n\n");
	
	assert (checkTree ("givenTree", "check_syn_rate"), 
			"Incorrectly generated constraints on synRate from ReplicateConstraint ('this1.?.synRate := this2.?.nonSynRate', givenTree, third_tree)"
		   );
	
	ReplicateConstraint ('this1.?.synRate := this2.?.nonSynRate', givenTree, third_tree);
	// should generate nothing

	assert (Abs (LAST_SET_OF_CONSTRAINTS) == 0, 
			"Second application of ReplicateConstraint ('this1.?.synRate := this2.?.nonSynRate', givenTree, third_tree) generated new constraints"
		   );


	ReplicateConstraint ('this1.?.? := this2.?.?', givenTree, third_tree);
	assert (checkTree ("givenTree", "check_both"), 
			"Incorrectly generated constraints on nonSynRate from ReplicateConstraint ('this1.?.synRate := this2.?.nonSynRate', givenTree, third_tree)"
		   );
	
	ClearConstraints (givenTree);
	ReplicateConstraint ('this1.?.? := this1.?.synRate__', givenTree);
	//fprintf (stdout, LAST_SET_OF_CONSTRAINTS, "\n\n");
	assert (checkTree ("givenTree", "check_syn_rate0"), 
			"Incorrectly generated constraints on synRate from ReplicateConstraint ('this1.?.? := this1.?.synRate__', givenTree)"
		   );
	   
	ReplicateConstraint ('this1.?.? := this1.?.?__', givenTree);
	assert (checkTree ("givenTree", "check_both0"), 
			"Incorrectly generated constraints on nonSynRate from ReplicateConstraint ('this1.?.? := this1.?.?__', givenTree)"
		   );

	ClearConstraints (givenTree);
	ReplicateConstraint ('this1.?.synRate := this2.?.synRate / this3.?.nonSynRate', fourth_tree, second_tree.Node8, givenTree.Node8);

	assert (checkTree ("fourth_tree", "check_syn_rate2"), 
			"Incorrectly generated constraints on synRate from ReplicateConstraint ('this1.?.synRate := this2.?.synRate / this3.?.nonSynRate', fourth_tree, second_tree.Node8, givenTree.Node8)"
		   );


	testResult = 1;

	return testResult;
}
