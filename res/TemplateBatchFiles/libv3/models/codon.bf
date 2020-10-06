LoadFunctionLibrary ("../tasks/genetic_code.bf");
LoadFunctionLibrary ("../UtilityFunctions.bf");

/** @module models.codon */

/**
 * @name models.codon.MapCode
 * @param {String} genetic_code
 * @returns {Dictionary} the sense, stop, and translation-table for the genetic code
 */
function models.codon.MapCode (genetic_code) {

	return {terms.sense_codons : utility.UniqueValues (genetic_code.ComputeCodonCodeToStringMap (genetic_code)),
	        terms.stop_codons  : utility.UniqueValues (genetic_code.ComputeCodonCodeToStringMapStop (genetic_code)),
	        terms.translation_table : genetic_code.DefineCodonToAAGivenCode (genetic_code) };
}

/**
 * @name models.codon.generic.DefineQMatrix
 * @param modelSpect
 * @param namespace
 */
function models.codon.generic.DefineQMatrix (modelSpec, namespace) {
//#profile START;

	__base_alphabet = modelSpec [terms.bases];
	assert (Type (__base_alphabet) == "Matrix" && Columns (__base_alphabet) == 4, "Unsupported codon bases '" + __base_alphabet + "'");

	__alphabet      = modelSpec [terms.alphabet];
	// need at least one codon per amino-acid

	__dimension     = model.Dimension (modelSpec);

	__stops         = modelSpec [terms.stop_codons]; // "stop"
	__table         = modelSpec [terms.translation_table];

	assert (Type (__alphabet) == "Matrix" && __dimension > 20 && __dimension + Columns (__stops) && Abs (__table) == 64, "Unsupported or missing alphabet '" + __alphabet + "' (stop codons :'" + __codons + "')");

	__modelType = modelSpec[terms.model.type];
	if (Type (__modelType) == "None" || Type (__modelType) == "Number") {
		__modelType = terms.global;
	}
	modelSpec[terms.model.type] = __modelType;
	assert (__modelType == terms.local || __modelType == terms.global, "Unsupported or missing model type '" + __modelType + "'");

	__rate_function = modelSpec [terms.model.q_ij];
	assert (utility.IsFunction (__rate_function), "Missing q_ij callback in model specification");

	__time_function = modelSpec [terms.model.time];
	assert (utility.IsFunction (__time_function), "Missing time callback in model specification");

	__rate_matrix = {__dimension,__dimension};
	__rate_matrix [0][0] = "";

	__rate_variation = model.generic.get_rate_variation (modelSpec);

	__global_cache = {};


	if (None != __rate_variation) {
		__rp = Call (__rate_variation[terms.rate_variation.distribution], __rate_variation[terms.rate_variation.options], namespace);
		__rate_variation [terms.id] = (__rp[terms.category])[terms.id];
		parameters.DeclareCategory   (__rp[terms.category]);
        parameters.helper.copy_definitions (modelSpec[terms.parameters], __rp);
	}


	for (_rowChar = 0; _rowChar < __dimension; _rowChar +=1 ){
		for (_colChar = _rowChar + 1; _colChar < __dimension; _colChar += 1) {


            __rp = Call (__rate_function, __alphabet[_rowChar],
                                          __alphabet[_colChar],
                                          namespace,
                                          __modelType,
                                          modelSpec);

            if (Abs (__rp) == 0) { // null rate
                continue;
            }

		 	if (None != __rate_variation) {
				__rp = Call (__rate_variation[terms.rate_variation.rate_modifier],
							 __rp,
							 __alphabet[_rowChar],
							 __alphabet[_colChar],
							 namespace,
							 __rate_variation [terms.id]);
 			}

            if (Abs (__rp[terms.model.rate_entry])) {
                parameters.DeclareGlobal (__rp[terms.global], __global_cache);
                parameters.helper.copy_definitions (modelSpec[terms.parameters], __rp);
                __rate_matrix [_rowChar][_colChar] = __rp[terms.model.rate_entry];
                __rate_matrix [_colChar][_rowChar] = __rp[terms.model.rate_entry];
            }
		}
	}

	if (__modelType == terms.global) {
	    __rp = Call (__time_function, __modelType);

        if (Type (__rp) == "String") {
            if (Abs (__rp)) {
                ((modelSpec[terms.parameters])[terms.local])[terms.parameters.synonymous_rate] = __rp;
                modelSpec [terms.model.rate_matrix] = parameters.AddMultiplicativeTerm (__rate_matrix, __rp, FALSE);
            } else {
                modelSpec [terms.model.rate_matrix] = __rate_matrix;
            }
        } else {
            if (Type (__rp) == "AssociativeList") {
                parameters.DeclareGlobal (__rp[terms.global], __global_cache);
                parameters.helper.copy_definitions (modelSpec[terms.parameters], __rp);
                modelSpec [terms.model.rate_matrix] = parameters.AddMultiplicativeTerm (__rate_matrix, __rp[terms.model.rate_entry], FALSE);
            }
        }
    } else {
        modelSpec [terms.model.rate_matrix] = __rate_matrix;
    }

/*
#profile _hyphy_profile_dump;



stats  			= _hyphy_profile_dump["STATS"];
_profile_summer = ({1,Rows(stats)}["1"]) * stats;
_instructions   = _hyphy_profile_dump["INSTRUCTION"];
_indices	    = _hyphy_profile_dump["INSTRUCTION INDEX"];

fprintf (stdout, "\nTotal run time (seconds)      : ", Format(_profile_summer[1],15,6),
                 "\nTotal number of steps         : ", Format(_profile_summer[0],15,0), "\n\n");

to_sort        =  stats["-_MATRIX_ELEMENT_VALUE_*_MATRIX_ELEMENT_COLUMN_+(_MATRIX_ELEMENT_COLUMN_==0)*_MATRIX_ELEMENT_ROW_"] % 1;

for (k=0; k<Columns(_instructions); k=k+1)
{
    k2 = to_sort[k][0];
    fprintf (stdout, Format (_indices[k2],6,0), " : ", _instructions[k2], "\n\tCall count: ", stats[k2][0],
                                                   "\n\tTime (seconds): ", stats[k2][1], "\n");
}
*/

    return modelSpec;

}

/**
 * Generates a branch length extraction stencil for synonymous and non-synonymous rates
 * @name models.codon.generate_stencil 
 * @param {String} type - terms.genetic_code.synonymous or terms.genetic_code.nonsynonymous or callback
 * @param {Dictionary} model - the model object
 * @returns {Matrix} - the appropriate NxN matrix stencil
*/

lfunction models.codon.generate_stencil._non_synonymous (from, to, model) {
    return (model[^"terms.translation_table"])[from] != (model[^"terms.translation_table"])[to];
}

lfunction models.codon.generate_stencil._synonymous (from, to, model) {
    return (model[^"terms.translation_table"])[from] == (model[^"terms.translation_table"])[to];
}


lfunction models.codon.generate_stencil (type, model) {
    
    if (type == utility.getGlobalValue("terms.genetic_code.synonymous")) {
        callback = "models.codon.generate_stencil._synonymous";
    } else {
        if (type == utility.getGlobalValue("terms.genetic_code.non_synonymous")) {
            callback = "models.codon.generate_stencil._non_synonymous";
        } else {
            callback = type;
        }
    }
 
	__dimension     = model.Dimension (model);
	__alphabet      = model [^"terms.alphabet"];
    stencil         = {__dimension,__dimension};
    
	for (_rowChar = 0; _rowChar < __dimension; _rowChar +=1 ){
		for (_colChar = _rowChar + 1; _colChar < __dimension; _colChar += 1) {
            stencil [_rowChar][_colChar] = Call (callback, __alphabet[_rowChar], __alphabet[_colChar], model);
            stencil [_colChar][_rowChar] = Call (callback, __alphabet[_colChar], __alphabet[_rowChar], model);
        }
    }
    
    return stencil;
}

lfunction models.codon.diff (a,b) {

    r = {};
    r [^"terms.diff.position"] = None;
    
    for (i = 0; i < 3; i += 1) {
        if (a[i] != b[i]) {
            if (r[^"terms.diff.position"] != None) {
                return None;
            }
            r[^("terms.diff.from")] = a[i];
            r[^("terms.diff.to")] = b[i];
            r[^("terms.diff.position")] = i;
        }
    }
    if (r[^"terms.diff.position"] == None)  {
        return None;
    }
    return r;
    
}


/** return complete differences between two codons **/

lfunction models.codon.diff.complete (a,b) {
    r = {};

    for (i = 0; i < 3; i += 1) {
        if (a[i] != b[i]) {
            c = {};
            c[^"terms.diff.from"] = a[i];
            c[^"terms.diff.to"] = b[i];
            c[^"terms.diff.position"] = i;
            r+c;
        }
    }

    return r;
}
