ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
  return "^";
}		

function power(x, n) {
  pow = 1;
  for (i=0; i<n; i=i+1) {
    pow = x*pow;
  }
  return pow;
}

function runTest () {
	testResult = 0;
  
  //---------------------------------------------------------------------------------------------------------
  // BASIC FUNCTIONALITY
  //---------------------------------------------------------------------------------------------------------
  // Integer powers
  x = {{ 1, 1.5, 2, 3, 4.71, 5.33 }};
  n = {{ 0, 1, 2, 3, 4, 5 }};
  for (i=0; i<6; i=i+1) {
    for (j=0; j<6; j=j+1) {
      caret = x[i]^n[j];
      comparison = power(x[i], n[j]);
      assert(Abs(caret - comparison) < 1e-12, "Fails to agree with a simple function that calculates integer powers");
    }
  }
  // Assorted fractional powers
  assert(4^(1/2) == 2, "Fails to compute square root");
  assert(8^(1/3) == 2, "Fails to compute cube root");
  assert(16^(1/4) == 2, "Fails to compute fourth root");
  assert(Abs(4^(3/2) - 8) < 1e-8, "Fails to compute fractional power");
  
  // Poisson log-likelihood for arrays
  m = {{1, 2, 3, 4}};
  assert(Abs(m^2 - (-6.73149)) < 1e-4, "Fails to compute Poisson log-likelihood of a vector");

  // String regexing
  assert("abc"^{{"b","d"}} == "adc", "Fails to handle a single character replacement regex");
  assert("abcd"^{{"b|d","e"}} == "aece", "Fails to handle a logical or character replacement regex");

  // TODO: avl_true comes from fprintf(stdout, avl_caret). Why do a,b and c all have the same parent?
  Topology T = ((a:.2,b:.4):.6,c:.8);
  avl_caret = T^1;
  avl_true = {
   "0":{
     "Name":"T",
     "Root":1
    },
   "1":{
     "Children":{
       "0":2,
       "1":3,
       "2":4
      },
     "Depth":0,
     "Length":0.6,
     "Name":"Node1"
    },
   "2":{
     "Depth":1,
     "Length":0.2,
     "Name":"a",
     "Parent":1
    },
   "3":{
     "Depth":1,
     "Length":0.4,
     "Name":"b",
     "Parent":1
    },
   "4":{
     "Depth":1,
     "Length":0.8,
     "Name":"c",
     "Parent":1
    }
  };

  //TODO: Operation '==' is not implemented/defined for a AssociativeList; would be nice instead of the following copy/paste
  //TODO: assert(avl_caret == avl_true, "failed test");
  assert((avl_caret["0"])["Name"] == (avl_true["0"])["Name"], "AVL tree does not match pre-computed value");
  assert((avl_caret["0"])["Root"] == (avl_true["0"])["Root"], "AVL tree does not match pre-computed value");
  assert((avl_caret["1"])["Depth"] == (avl_true["1"])["Depth"], "AVL tree does not match pre-computed value");
  assert((avl_caret["1"])["Length"] == (avl_true["1"])["Length"], "AVL tree does not match pre-computed value");
  assert((avl_caret["1"])["Name"] == (avl_true["1"])["Name"], "AVL tree does not match pre-computed value");
  assert(((avl_caret["1"])["Children"])["0"] == ((avl_true["1"])["Children"])["0"], "AVL tree does not match pre-computed value");
  assert(((avl_caret["1"])["Children"])["1"] == ((avl_true["1"])["Children"])["1"], "AVL tree does not match pre-computed value");
  assert(((avl_caret["1"])["Children"])["2"] == ((avl_true["1"])["Children"])["2"], "AVL tree does not match pre-computed value");
  assert((avl_caret["2"])["Depth"] == (avl_true["2"])["Depth"], "AVL tree does not match pre-computed value");
  assert((avl_caret["2"])["Length"] == (avl_true["2"])["Length"], "AVL tree does not match pre-computed value");
  assert((avl_caret["2"])["Name"] == (avl_true["2"])["Name"], "AVL tree does not match pre-computed value");
  assert((avl_caret["2"])["Parent"] == (avl_true["2"])["Parent"], "AVL tree does not match pre-computed value");
  assert((avl_caret["3"])["Depth"] == (avl_true["3"])["Depth"], "AVL tree does not match pre-computed value");
  assert((avl_caret["3"])["Length"] == (avl_true["3"])["Length"], "AVL tree does not match pre-computed value");
  assert((avl_caret["3"])["Name"] == (avl_true["3"])["Name"], "AVL tree does not match pre-computed value");
  assert((avl_caret["3"])["Parent"] == (avl_true["3"])["Parent"], "AVL tree does not match pre-computed value");
  assert((avl_caret["4"])["Depth"] == (avl_true["4"])["Depth"], "AVL tree does not match pre-computed value");
  assert((avl_caret["4"])["Length"] == (avl_true["4"])["Length"], "AVL tree does not match pre-computed value");
  assert((avl_caret["4"])["Name"] == (avl_true["4"])["Name"], "AVL tree does not match pre-computed value");
  assert((avl_caret["4"])["Parent"] == (avl_true["4"])["Parent"], "AVL tree does not match pre-computed value");

  testResult = 1;

  return testResult;
}
