LoadFunctionLibrary ("chooseGeneticCode", {"0" : "Universal"});

function models.codon.map_code (genetic_code) {
	return {"sense" : Columns (ComputeCodonCodeToStringMap (genetic_code)), 
	        "stop"  : Columns (ComputeCodonCodeToStringMapStop (genetic_code))};
}



