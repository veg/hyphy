ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "Model";
}		


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = TRUE;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Model id=(inst_transition_matrix_ident, equilibrium_frequencies_ident, <multiply by frequencies>);

  // Make sure you can create a simple model.
  // Testing of model usage is handled by the functions that consume or use models.
  global kappa = 4;
  Q_HKY85 = {{*,t,kappa*t,t}
           {t,*,t,t*kappa}
           {t*kappa,t,*,t}
           {t,t*kappa,t,*}};
  freqs = {{0.4}{0.3}{0.2}{0.1}};
  Model HKYd85 = (Q_HKY85, freqs, 1);


  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------
  assert (runCommandWithSoftErrors ('Model oneArg = (Q_HKY85);', "missing in Model definition. Must have a matrix and a compatible eqiulibrium frequencies vector"), "Failed error checking for trying to create a model with only one parameter");

  freqs2 = {{0.25}{0.25}{0.25}};
  assert (runCommandWithSoftErrors ('Model oneArg = (Q_HKY85, freqs2);', "must be a column vector of the same dimension as the model matrix in the call to Model"), "Failed error checking for trying to create a model with a frequencies vector of different dimensions to the transition matrix");
  
  rectangularTransitionMatrix = {{*,t,kappa*t,t}
                                {t,*,t,t*kappa}
                                {t*kappa,t,*,t}
                                {t,t*kappa,t,*}
                                {t,t*kappa,t,*}};
  assert (runCommandWithSoftErrors ('Model rectTransMatrix = (rectangularTransitionMatrix, freqs);', "must be a square matrix of dimension>=2 in the call to Model"), "Failed error checking for trying to create a model with a rectangular transition matrix");

  testResult = 1;

  return testResult;
}
