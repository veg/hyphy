ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
return runATest ();


function getTestName () {
  return "Polynomial";
}		

function test_a_polynomial (S) {
    test_r = 0;
    ExecuteCommands ("_LAMBDA_:=" + S);
    ExecuteCommands ("_POLY_ = Polynomial ('`S`');");
    ExecuteCommands ("_POLY_2_ = Polynomial ('`S`');");
    GetString (vars, _LAMBDA_, -1);
    local_vars = vars["Local"];
    
    assert (_POLY_ == _POLY_2_, "Polynomial `S` equal to self");
    
    _POLY_SQ_ = _POLY_ * _POLY_;
    
    for (i = 0; i < 100; i+=1) {
        for (vn; in; local_vars) {
            ^vn = Random (-1, 1);
        }
        _ep = Eval (_POLY_);
        assert (Abs (_ep-_LAMBDA_) < 1e-10, "Polynomial `S` evaluates correctly");
        assert (Abs(_ep^2 - Eval (_POLY_SQ_)) < 1e-10, "Poly == Poly^2");
     }
    
    assert (_POLY_ == _POLY_2_, "Polynomial `S` equal to self following round-trip conversion");
    test_r = 1;
    return test_r;
}



function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult = 0;
  

  //---------------------------------------------------------------------------------------------------------
  // SIMPLE FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
    
  
 


  // test that polynomials generate without errors
  PN=Polynomial (""); // null (0) polynomial
  P0=Polynomial (0);
  P1=Polynomial (Exp(1));
  P2=Polynomial ("x+y+z^3");
  P3=Polynomial ("y+z^3+x");
  P4=Polynomial ("2a+3b+4c");
  P23=Polynomial ("2z^3+2x+2y");
  assert (Type(PN)=='Polynomial', "Creating a null polynomial");
  assert (Type(P0)=='Polynomial', "Creating a 0 polynomial");
  assert (Type(P1)=='Polynomial', "Creating an 'e' polynomial");
  assert (Type(P2)=='Polynomial', "Creating an 'x+y+z' polynomial");
  assert (Type(P3)=='Polynomial', "Creating an 'y+x+z' polynomial");
  assert (Type(P4)=='Polynomial', "Creating an '2a+3b+4b' polynomial");  
  P2x=Polynomial(P2);
  assert (P2x==P2, "P2==P2"); 
  assert (P0==P0, "0==0"); 
  assert (PN+Exp(1)==P1, "0+Exp(1)==Exp(1)"); 
  assert (P2==P2+PN, "x==x+0"); 
  assert (P1==P1, "Exp(1)==Exp(1)"); 
  assert (P2==P2, "x+y+z==x+y+z"); 
  assert (P2==P3, "x+y+z==y+z+x"); 
  assert (P23==P2+P3, "P23==P2+P3"); 
  assert (P2+P4==Polynomial ("x+y+z^3+2a+3b+4c"), "P2+P4"); 
  assert (P0+P4-P4==Polynomial (0), "P0+P4-P4=P0"); 
  assert (P0+P1==P1+P0, "P0+P1==P1+P0"); 
  assert (P2*P4==P4*P2, "P2*P4==P4*P2"); 
  assert (Polynomial ("x+x+x-3x")==Polynomaial (0), "x+x+x-3x==0"); 
  assert (Eval (P23) == Eval(P2) + Eval(P3), "Eval works for non-trivial polynomials");
  
  // checking evaluation blocks 
  
  x = 1;
  y = 2; 
  z = 3;
  w = 4;
  
  E = Polynomial (0);
  
  test_a_polynomial ("1+x+2x^2+3x^3");
  
  
  assert (P2+P3==P23, "x+y+z+y+z+x==z+2x+2y"); 
  assert (Eval (P0*P4) == 0, "0*(poly)==0");
  assert (Eval (P1) == Exp(1), "Eval works");
  x = 1; y = 2; z = 3;

  p_list = {{
        "0",
        "2.71",
        "x^2",
        "1+x+x^2+x^3",
        "1+x+4x^3+x^5",
        "x*y^2",
        "x*y+y^2+z^2+x*y*z",
        "x+2y+3z",
        "a*b*c+x^2+y^2"}};
  
  
   for (_p; in; p_list) {
    assert (test_a_polynomial (_p), _p);
   }


  //---------------------------------------------------------------------------------------------------------
  // ERROR HANDLING
  //---------------------------------------------------------------------------------------------------------

 
  assert (runCommandWithSoftErrors ('Polynomial ({{1,1}})', "Operation 'Polynomial' is not implemented/defined for a Matrix"), "Failed error checking for Polynomial ({{1,1}})");

  testResult = 1;

  return testResult;
}

