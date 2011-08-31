# import the HyPhy library and R glue
# change the paths according to your distribution

dyn.load ("LibraryModules/R/HyPhy.so")
source  ("LibraryModules/R/HyPhy.R")

# first, create a HyPhy interface instance (class _THyPhy)
# the first argument defines the root directory for HyPhy
# and the second - how many threads the computational core
# should spawn

hyphyInstance<-`_THyPhy` (paste(getwd(),'/',sep=''),2)

# the basic interface command is 'ExecuteBF' which
# executes HyPhy batch language commands in HyPhy
# and returns a string representation of the return value
# (if any) from HYPHY
# The returned object is of type _THyPhyString with 
# sData and sLength fields
# HyPhy will take care of disposing of the memory needed 
# to store the result

hyphyResult<-hyphyInstance$ExecuteBF (hyphyInstance,"return 2+2;")
print(paste("Testing a trivial HyPhy command. 2+2 = ",hyphyResult$sData))

# an optional second argument to ExecuteBF
# can be used to "flush" the current state of the system

# this is the default option for the call of ExecuteBF
# passing the second argument of FALSE or 0 will preserve
# the execution state 

print("Consecutive command exection")
hyphyInstance$ExecuteBF (hyphyInstance,"z:=x+y;",FALSE)
hyphyInstance$ExecuteBF (hyphyInstance,"x=3;",FALSE)
hyphyInstance$ExecuteBF (hyphyInstance,"y=5;",FALSE)
hyphyResult = hyphyInstance$ExecuteBF (hyphyInstance,"return z;",FALSE)
print(paste("The value of z is ",hyphyResult$sData))

print("Resetting the state of the execution erases the value of 'z'")
hyphyResult<-hyphyInstance$ExecuteBF (hyphyInstance,"return z;");
print(paste("The value of z is ",hyphyResult$sData))

# the real utility of the interface is to be able 
# to execute prewritten analyses from HBL files

print("Executing the example F81.bf file")
hyphyResult<-hyphyInstance$ExecuteBF (hyphyInstance,"ExecuteAFile(\"Examples/HBL/F81.bf\")");

# retrive the standard output, error and runtime warnings

hyphyOut 		<- hyphyInstance$GetStdout(hyphyInstance)
#errors will be empty UNLESS there was an exection error
hyphyErrors     <- hyphyInstance$GetErrors(hyphyInstance)
hyphyWarnings	<- hyphyInstance$GetWarnings(hyphyInstance)

print (paste("Standard out: \n", hyphyOut$sData))
print (paste("Errors: \n", hyphyErrors$sData))
print (paste("Warnings/Log messages: \n", hyphyWarnings$sData))

# these variables can be explicitly deleted when they are no longer needed
# R garbage collection should take care of disposing of disused variables

rm('hyphyOut')
rm('hyphyErrors')
rm('hyphyWarnings')

# A tighter intergration can be achieved by defining a retrieval function 
# with the reserved name _THyPhyAskFor in the HBL file; it retrieves data
# by key and returns them in a internal format that can be converted 
# to one of the basic return types: number(0), string(1) or matrix(2)

retrieveValueByKey <- function(key, returnType, hyphyInstance){
	theResult <- hyphyInstance$AskFor(hyphyInstance,key);
	# see if HyPhy can retrieve a value with the requested key
	if (!is.null(theResult))
	{
  		canICast <- hyphyInstance$CanCast(hyphyInstance,theResult,returnType);
   		# see if HyPhy can cast the value to the requested type
   		if(canICast)
   		{
   			# do the casting
   			theResult <- hyphyInstance$CastResult(hyphyInstance,theResult,returnType);
   			# the last step is to convert from the basic return type
   			# to a derived class that we can use in python directly
   			if (returnType == 0)
   			{
   				return(theResult$castToNumber(theResult));
   			}
   			if (returnType == 1)
   			{
   				return(theResult$castToString(theResult));
   			}
   			if (returnType == 2)
   			{
   				return(theResult$castToMatrix(theResult));
   			}
  		}
  	}
	return(NULL);
}
 		
   		

hyphyResult<-hyphyInstance$ExecuteBF (hyphyInstance,"ExecuteAFile(\"Examples/HBL/HKY85.bf\")");

print(paste("Log-L = ", retrieveValueByKey ("LogL", 0, hyphyInstance)$nValue))
print(paste("kappa = ", retrieveValueByKey ("kappa", 0, hyphyInstance)$nValue))
print(paste("tree string = ", retrieveValueByKey ("Tree", 1, hyphyInstance)$sData))
bl<-retrieveValueByKey ("Branch lengths", 2, hyphyInstance)
print(paste("retrieved ", bl$mCols-1, "branch lengths"))
for(i in 0:(bl$mCols-2))
{
	print (paste("Branch ", i+1, " has length ", bl$MatrixCell(bl,0,i)))
}











