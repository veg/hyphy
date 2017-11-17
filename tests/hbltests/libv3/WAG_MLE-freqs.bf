LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/models/protein/empirical.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");

file_name = "`PATH_TO_CURRENT_BF`data/CD2.prot";
wag_mlf.data = {};
wag_mlf.filter = {};

wag_mlf.alignment_info = alignments.ReadNucleotideAlignment(file_name, "wag_mlf.data", "wag_mlf.filter");

wag_mlf.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions(wag_mlf.alignment_info[terms.data.partitions], {});

wag_mlf.partitions = alignments.DefineFiltersForPartitions (wag_mlf.partitions_and_trees,  "wag_mlf.data" , "wag_mlf.filter.", wag_mlf.alignment_info);

trees = utility.Map (wag_mlf.partitions_and_trees, "_partition_", '_partition_[terms.data.tree]');
filter_names = utility.Map (wag_mlf.partitions, "_partition_", '_partition_[terms.data.name]');


                                         
wag_freq.results = estimators.FitSingleModel_Ext(filter_names,
                                         trees,
                                         "models.protein.WAGF.ModelDescription",
                                         None,None);

console.log ("WAG + F, Log(L) => " + wag_freq.results[terms.fit.log_likelihood]);

wag_mlf.results = estimators.FitSingleModel_Ext(filter_names,
                                         trees,
                                         "models.protein.WAGML.ModelDescription",
                                         None,None);

console.log ("WAG + ML, Log(L) => " + wag_mlf.results[terms.fit.log_likelihood]);
