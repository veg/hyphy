# import the HyPhy library
# and standard OS utilities

import os, HyPhy

# first, create a HyPhy interface instance (class _THyPhy)
# the first argument defines the root directory for HyPhy
# and the second - how many threads the computational core
# should spawn

hyphyInstance = HyPhy._THyPhy (os.getcwd(),2)

# the basic interface command is 'ExecuteBF' which
# executes HyPhy batch language commands in HyPhy
# and returns a string representation of the return value
# (if any) from HYPHY
# The returned object is of type _THyPhyString with 
# sData and sLength fields
# HyPhy will take care of disposing of the memory needed 
# to store the result

hyphyResult = hyphyInstance.ExecuteBF ("return 2+2;");
print ("Testing a trivial HyPhy command. 2+2 = ", hyphyResult.sData)

# an optional second argument to ExecuteBF
# can be used to "flush" the current state of the system

# this is the default option for the call of ExecuteBF
# passing the second argument of False or 0 will preserve
# the execution state 

print ("Consecutive command exection")
hyphyInstance.ExecuteBF ("z:=x+y;",False);
hyphyInstance.ExecuteBF ("x=3;",False);
hyphyInstance.ExecuteBF ("y=5;",False);
hyphyResult = hyphyInstance.ExecuteBF ("return z;",False);
print ("The value of z is ", hyphyResult.sData)

print ("Resetting the state of the execution erases the value of 'z'")
hyphyResult = hyphyInstance.ExecuteBF ("return z;");
print ("The value of z is ", hyphyResult.sData)

# the real utility of the interface is to be able 
# to execute prewritten analyses from HBL files

print ("Executing the example F81.bf file")
hyphyInstance.ExecuteBF ("ExecuteAFile(\"../HBL/F81.bf\")");

# retrive the standard output, error and runtime warnings

hyphyOut 		= hyphyInstance.GetStdout()
#errors will be empty UNLESS there was an exection error
hyphyErrors     = hyphyInstance.GetErrors()
hyphyWarnings	= hyphyInstance.GetWarnings()

print ("Standard out: \n", hyphyOut.sData)
print ("Errors: \n", hyphyErrors.sData)
print ("Warnings/Log messages: \n", hyphyWarnings.sData)

# these variables can be explicitly deleted when they are no longer needed
# python garbage collection should take care of disposing of disused variables

del hyphyOut
del hyphyErrors
del hyphyWarnings

# A tighter intergration can be achieved by defining a retrieval function 
# with the reserved name _THyPhyAskFor in the HBL file; it retrieves data
# by key and returns them in a internal format that can be converted 
# to one of the basic return types: number, string or matrix

def retrieveValueByKey (key, returnType,hyphyInstance):
	theResult = hyphyInstance.AskFor(key)
	# see if HyPhy can retrieve a value with the requested key
	if theResult:
   		canICast = hyphyInstance.CanCast(theResult,returnType)
   		# see if HyPhy can cast the value to the requested type
   		if canICast:
   			# do the casting
   			theResult = hyphyInstance.CastResult(theResult,returnType)
   			# the last step is to convert from the basic return type
   			# to a derived class that we can use in python directly
   			if (returnType == HyPhy.THYPHY_TYPE_NUMBER):
   				return theResult.castToNumber()
   			if (returnType == HyPhy.THYPHY_TYPE_STRING):
   				return theResult.castToString()
   			if (returnType == HyPhy.THYPHY_TYPE_MATRIX):
   				return theResult.castToMatrix()
	return null   		
   		

hyphyResult = hyphyInstance.ExecuteBF ("ExecuteAFile(\"../HBL/HKY85.bf\")");


print ("Log-L = ", retrieveValueByKey ("LogL", HyPhy.THYPHY_TYPE_NUMBER, hyphyInstance).nValue)
print ("kappa = ", retrieveValueByKey ("kappa", HyPhy.THYPHY_TYPE_NUMBER, hyphyInstance).nValue)
print ("tree string = ", retrieveValueByKey ("Tree", HyPhy.THYPHY_TYPE_STRING, hyphyInstance).sData)
bl = retrieveValueByKey ("Branch lengths", HyPhy.THYPHY_TYPE_MATRIX, hyphyInstance)
print ("retrieved ", bl.mCols-1, "branch lengths")
for i in range(0,bl.mCols-1):
	print ("Branch ", i+1, " has length ", bl.MatrixCell(0,i))












