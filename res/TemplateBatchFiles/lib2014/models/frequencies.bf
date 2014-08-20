LoadFunctionLibrary ("terms.bf");
LoadFunctionLibrary ("parameters.bf");
LoadFunctionLibrary("../UtilityFunctions.bf");

function frequencies.equal (model, namespace, datafilter) {
	__N = Abs (model["alphabet"]);
	model["EFV"] = {__N,1}["1/__N"];
	return model;
}

function frequencies.empirical.nucleotide (model, namespace, datafilter) {
	__N = Abs (model["alphabet"]);
	utility.toggleEnvVariable ("COUNT_GAPS_IN_FREQUENCIES", 0);
	HarvestFrequencies (__f, *datafilter, 1,1,1);
	utility.toggleEnvVariable ("COUNT_GAPS_IN_FREQUENCIES", None);
	model[terms.efv_estimate] = __f;
	model[terms.efv_estimate_name] = terms.freqs.4x1;
	return model;
}

function frequencies.mle (model, namespace, datafilter) {
	assert (0, "frequencies.mle is TBD");
}
