RequireVersion("2.3");

/*------------------------------------------------------------------------------
    Load library files
*/

LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");
LoadFunctionLibrary("libv3/terms-json.bf");

LoadFunctionLibrary("libv3/tasks/ancestral.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");

LoadFunctionLibrary("libv3/models/DNA/GTR.bf");
LoadFunctionLibrary("libv3/models/rate_variation.bf");

LoadFunctionLibrary("libv3/convenience/math.bf");

/*------------------------------------------------------------------------------ 
Display analysis information
*/

io.DisplayAnalysisBanner({
    "info": "A simple demonstration of how rate variation can be added onto standard models",
    "version": "0.01",
    "authors": "Sergei L Kosakovsky Pond",
    "contact": "spond@temple.edu",
    "requirements": "a nucleotide alignment"
});


msa = alignments.ReadNucleotideAlignment (utility.GetEnvVariable ("PATH_TO_CURRENT_BF") + "/data/CD2.nex", "data", "filter");

name_mapping = msa[utility.getGlobalValue("terms.json.name_mapping")];
if (None == name_mapping) { /** create a 1-1 mapping if nothing was done */
	name_mapping = {};
	utility.ForEach (alignments.GetSequenceNames ("filter"), "_value_", "`&name_mapping`[_value_] = _value_");
}

partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (msa[utility.getGlobalValue("terms.json.partitions")], name_mapping);


trees = utility.Map (partitions_and_trees, "_partition_", '_partition_["tree"]');
io.ReportProgressMessageMD ("rv", "nuc-fit", "Obtaining branch lengths and nucleotide rates under the constant rate GTR model");
gtr_results = estimators.FitSingleModel_Ext({{"filter"}},
									 trees,
									 "models.DNA.GTR.ModelDescription",
									 parameters.helper.tree_lengths_to_initial_values (trees, None),
									 None);


io.ReportProgressMessageMD ("rv", "nuc-fit", "* Log(L) = " + Format (gtr_results["LogL"], 8, 2));

//------------------------------------------------------------------------------------------------------------------------
lfunction models.DNA.GTR.ModelDescription.withGamma (options) {
	def = models.DNA.GTR.ModelDescription (options);
	def [^"terms.rate_variation"] = rate_variation.types.Gamma.factory ({^"terms.rate_variation.bins" : 4});
	return def;
};

io.ReportProgressMessageMD ("rv", "nuc-fit", "Obtaining branch lengths and nucleotide rates under the GTR+G/4 model");
gtr_results_gamma = estimators.FitSingleModel_Ext({{"filter"}},
									 trees,
									 "models.DNA.GTR.ModelDescription.withGamma",
									 gtr_results,
									 None);
io.ReportProgressMessageMD ("rv", "nuc-fit", "* Log(L) = " + Format (gtr_results_gamma["LogL"], 8, 2));

//------------------------------------------------------------------------------------------------------------------------
lfunction models.DNA.GTR.ModelDescription.withGammaInv (options) {
	def = models.DNA.GTR.ModelDescription (options);
	def [^"terms.rate_variation"] = rate_variation.types.GammaInv.factory ({^"terms.rate_variation.bins" : 4});
	return def;
};

io.ReportProgressMessageMD ("rv", "nuc-fit", "Obtaining branch lengths and nucleotide rates under the GTR+G/4+I model");
gtr_results_gamma_i = estimators.FitSingleModel_Ext({{"filter"}},
									 trees,
									 "models.DNA.GTR.ModelDescription.withGammaInv",
									 gtr_results_gamma,
									 None);
io.ReportProgressMessageMD ("rv", "nuc-fit", "* Log(L) = " + Format (gtr_results_gamma_i["LogL"], 8, 2));

//------------------------------------------------------------------------------------------------------------------------

lfunction models.DNA.GTR.ModelDescription.withGDD2 (options) {
	def = models.DNA.GTR.ModelDescription (options);
	def [^"terms.rate_variation"] = rate_variation.types.GDD.factory ({^"terms.rate_variation.bins" : 2});
	return def;
};

io.ReportProgressMessageMD ("rv", "nuc-fit", "Obtaining branch lengths and nucleotide rates under the (2-bin) general discrete model");
gtr_results_gdd2 = estimators.FitSingleModel_Ext({{"filter"}},
									 trees,
									 "models.DNA.GTR.ModelDescription.withGDD2",
									 gtr_results_gamma_i,
									 None);
io.ReportProgressMessageMD ("rv", "nuc-fit", "* Log(L) = " + Format (gtr_results_gdd2["LogL"], 8, 2));

//------------------------------------------------------------------------------------------------------------------------

lfunction models.DNA.GTR.ModelDescription.withGDD4 (options) {
	def = models.DNA.GTR.ModelDescription (options);
	def [^"terms.rate_variation"] = rate_variation.types.GDD.factory ({^"terms.rate_variation.bins" : 4});
	return def;
};


io.ReportProgressMessageMD ("rv", "nuc-fit", "Obtaining branch lengths and nucleotide rates under the (4-bin) general discrete model");
gtr_results_gdd4 = estimators.FitSingleModel_Ext({{"filter"}},
									 trees,
									 "models.DNA.GTR.ModelDescription.withGDD4",
									 gtr_results_gdd2,
									 None);
io.ReportProgressMessageMD ("rv", "nuc-fit", "* Log(L) = " + Format (gtr_results_gdd4["LogL"], 8, 2));

