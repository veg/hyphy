/***** SJS Convenience batchfile for loading function libraries. *******/

RequireVersion("2.31"); // Should also be in the calling batchfile, but repeated here for additional sanity

LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf"); 
LoadFunctionLibrary("libv3/all-terms.bf");

LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");
LoadFunctionLibrary("libv3/tasks/ancestral.bf"); 

LoadFunctionLibrary("libv3/models/rate_variation.bf");

LoadFunctionLibrary("libv3/convenience/math.bf");



