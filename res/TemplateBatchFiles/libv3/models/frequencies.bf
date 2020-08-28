LoadFunctionLibrary("../all-terms.bf");
LoadFunctionLibrary("parameters.bf");
LoadFunctionLibrary("../UtilityFunctions.bf");
LoadFunctionLibrary("model_functions.bf");

/** @module frequencies */

/**
 * Sets model's equilibrium frequency estimator to Binary chatacter 2x1 estimator
 * @name frequencies.empirical.binary
 * @param {Dictionary} model
 * @param {String} namespace - does nothing
 * @param {DataSetFilter} datafilter - does nothing
 * @returns {Dictionary} updated model
 */
function frequencies.empirical.binary(model, namespace, datafilter) {
    model = frequencies._aux.empirical.singlechar(model, namespace, datafilter);
    model[terms.model.efv_estimate_name] = terms.frequencies.binary;
    return model;
}

/**
 * Sets model's equilibrium frequency estimator to equal frequencies
 * @name frequencies.equal
 * @param {Dictionary} model
 * @param {String} namespace - does nothing
 * @param {DataSetFilter} datafilter - does nothing
 * @returns {Dictionary} updated model
 */
function frequencies.equal(model, namespace, datafilter) {


    __N = utility.Array1D(model[terms.alphabet]);
    model[terms.efv_estimate] = {
        __N,
        1
    }["1/__N"];
    model[terms.model.efv_estimate_name] = terms.frequencies.equal;
    (model[terms.parameters])[terms.model.empirical] = 0;
    return model;
}

/**
 * Sets model's equilibrium frequency estimator to Nucleotide 4x1 estimator
 * @name frequencies.empirical.nucleotide
 * @param {Dictionary} model
 * @param {String} namespace
 * @param {DataSetFilter} datafilter
 * @returns {Dictionary} updated model
 */
function frequencies.empirical.nucleotide(model, namespace, datafilter) {
    model = frequencies._aux.empirical.singlechar(model, namespace, datafilter);
    model[terms.model.efv_estimate_name] = terms.frequencies._4x1;
    (model[terms.parameters])[terms.model.empirical] = 3;
    return model;
}

/**
 * Compute equilibrium frequencies at run-time using Q inversion
 * @name frequencies.empirical.nucleotide
 * @param {Dictionary} model
 * @param {String} namespace
 * @param {DataSetFilter} datafilter
 * @returns {Dictionary} updated model
 */
 
lfunction frequencies.runtime.nucleotide(model, namespace, datafilter) {

    /*
        The run-time calculation depends on the current value of the rate matrix and
        proceeds by 
        
        (1) setting the LAST column of the Q matrix to all 1s
        (2) inverting the modified matrix 
        (3) reading off equilibrium frequencies from the last row of the inverse
        
        This function needs to have the rate matrix instantiated before it can set 
        up the dependancies.        
    */
    
    
    return model;
}

/**
 * This is an function that computes stationary frequencies of Markov process based on its  
 * rate matrix
 * @param {Matrix} Q matrix
*/

lfunction frequencies._aux.invert_model (Q) {
    dim = Rows(Q);
    for (i = 0; i < dim; i++) {
        Q[i][dim-1] = 1;
    }
    return (Inverse (Q))[{{dim-1,0}}][{{dim-1,dim-1}}];
}


/**
 * Sets model's equilibrium frequency estimator to ML for binary data
 * @name frequencies.empirical.binary
 * @param {Dictionary} model
 * @param {String} namespace
 * @param {DataSetFilter} datafilter
 * @returns {Dictionary} updated model
 */
function frequencies.ML.binary (model, namespace, datafilter) {
    model = frequencies._aux.empirical.singlechar(model, namespace, datafilter);
    // define 2 frequency parameters
    // initialize to empirical freqs
    // add to model parameter manifest
    frequencies.ML.binary.emp = model[terms.efv_estimate];
    model[terms.efv_estimate] = {2,1};
    frequencies.ML.binary.variables = {2,1};
    frequencies.ML.binary.scaler = namespace + ".frequency_scaler";

    utility.ForEachPair (model[terms.alphabet], "_index", "_letter",
                         '
                            _idx = _index[1];
                            _freq_parameter = namespace + ".equilibrium_frequency_of." + _letter;
                            frequencies.ML.binary.variables [_idx] = _freq_parameter;
                            (model[terms.efv_estimate]) [_idx] = _freq_parameter + "/" + frequencies.ML.binary.scaler;
                            parameters.DeclareGlobalWithRanges (_freq_parameter, frequencies.ML.binary.emp[_idx], 0, 1);
                            model.generic.AddGlobal (model, _freq_parameter, terms.characterFrequency (_letter));
                         '
                         );

    // constrain pi_0 + pi_1 = 1,

    parameters.SetConstraint ( frequencies.ML.binary.scaler, Join ("+", frequencies.ML.binary.variables), "global");

    model[terms.model.efv_estimate_name] = terms.frequencies.ml;
    (model[terms.parameters])[terms.model.empirical] = -1; // correct for the restricted sum
    return model;
}

/**
 * Sets model's equilibrium frequency estimator to protein 20x1 estimator
 * @name frequencies.empirical.protein
 * @param {Dictionary} model
 * @param {String} namespace
 * @param {DataSetFilter} datafilter
 * @returns {Dictionary} updated model
 */
function frequencies.empirical.protein (model, namespace, datafilter) {
    model = frequencies._aux.empirical.singlechar(model, namespace, datafilter);
    model[terms.model.efv_estimate_name] = terms.frequencies._20x1;
    (model[terms.parameters])[terms.model.empirical] = 19;
    return model;
}

/**
 * Sets model's equilibrium frequency estimator to ML for protein data
 * @name frequencies.empirical.protein
 * @param {Dictionary} model
 * @param {String} namespace
 * @param {DataSetFilter} datafilter
 * @returns {Dictionary} updated model
 */
function frequencies.ML.protein (model, namespace, datafilter) {
    return frequencies.mle (model, namespace, datafilter);
}


/**
 * Multiplies nucleotide frequencies into the empirical codon estimate
 * @name frequencies.empirical.protein
 * @param {Dictionary} model
 * @param {Matrix} __estimates 4x3 matrix of frequency estimates
 * @returns {Dictionary} updated model
 */
function frequencies.codon.multiply_in_frequencies (model, __estimates) {

    __dimension = model.Dimension(model);
    __alphabet = model[^"terms.alphabet"];

    if (Type (model[terms.model.rate_matrix]) == "AssociativeList") { // handle mixture models
        __components = Rows (model[terms.model.rate_matrix]);
        __component_count = utility.Array1D (__components);

        for (_rowChar = 0; _rowChar < __dimension; _rowChar += 1) {
            for (_colChar = _rowChar + 1; _colChar < __dimension; _colChar += 1) {

                //__diff = models.codon.diff(__alphabet[_rowChar], __alphabet[_colChar]);
                for (__component_id = 0; __component_id < __component_count; __component_id += 1) {
                    if (Abs (((model[terms.model.rate_matrix])[__components[__component_id]]) [_colChar][_rowChar])) {
                        __diff = models.codon.diff.complete (__alphabet[_rowChar], __alphabet[_colChar]);
                         for (__diff_id = 0; __diff_id < Abs (__diff); __diff_id += 1) {
                             ((model[terms.model.rate_matrix])[__components[__component_id]]) [_rowChar][_colChar] += "*" + (__estimates[(__diff[__diff_id])["to"]])[(__diff[__diff_id])["position"]];
                             ((model[terms.model.rate_matrix])[__components[__component_id]]) [_colChar][_rowChar] += "*" + (__estimates[(__diff[__diff_id])["from"]])[(__diff[__diff_id])["position"]];
                         }
                    }
                }
            }
        }
    } else {
        for (_rowChar = 0; _rowChar < __dimension; _rowChar += 1) {
            for (_colChar = _rowChar + 1; _colChar < __dimension; _colChar += 1) {

                if (Abs ((model[terms.model.rate_matrix])[_colChar][_rowChar])) {
                    __diff = models.codon.diff.complete (__alphabet[_rowChar], __alphabet[_colChar]);
                    for (__diff_id, __diff_value; in; __diff) {
                        (model[terms.model.rate_matrix])[_rowChar][_colChar] += "*" + (__estimates[__diff_value["to"]])[__diff_value["position"]];
                        (model[terms.model.rate_matrix])[_colChar][_rowChar] += "*" + (__estimates[__diff_value["from"]])[__diff_value["position"]];
                    }
                }
            }
        }
    }
    return model;
}



lfunction frequencies._aux.validate (model) {

    __dimension = model.Dimension(model);
    __alphabet = model[^"terms.alphabet"];

    if (Type (model[^"terms.model.rate_matrix"]) == "AssociativeList") {
        utility.ForEach (model[^"terms.model.rate_matrix"], "_q_matrix_", '
            assert(Type(_q_matrix_) == "Matrix" && Rows(_q_matrix_) == __dimension && Columns(_q_matrix_) == __dimension,
            "`^'terms.model.rate_matrix'` must be defined prior to calling frequency estimators")
        ');
    } else {
        assert(Type(model[^"terms.model.rate_matrix"]) == "Matrix" && Rows(model[^"terms.model.rate_matrix"]) == __dimension && Columns(model[^"terms.model.rate_matrix"]) == __dimension,
            "`^'terms.model.rate_matrix'` must be defined prior to calling frequency estimators");
    }

}

lfunction frequencies.empirical.codon_from_nuc (model, nuc_dict) {

 
    result = {model.Dimension (model),1};
    corrector = 1;
    utility.ForEach (model[^"terms.stop_codons"], "codon",
    '
        `&corrector` += - (`&nuc_dict`[codon[0]])[0] * (`&nuc_dict`[codon[1]])[1] * (`&nuc_dict`[codon[2]])[2];
    ');

    utility.ForEachPair (model[^"terms.alphabet"], "index", "codon",
    '
        (`&result`)[index[1]] =  (`&nuc_dict`[codon[0]])[0] * (`&nuc_dict`[codon[1]])[1] * (`&nuc_dict`[codon[2]])[2] / `&corrector`;
    ');
    
    return (result);
}
/**
 * @name frequencies.empirical.F3x4
 * @param {Dictionary} model
 * @param {String} namespace
 * @param {DataSetFilter} datafilter
 * @returns {Dictionary} updated model
 */
 
lfunction frequencies.empirical.F3x4(model, namespace, datafilter) {

    frequencies._aux.validate (model);

    if (None != datafilter) {
         utility.ToggleEnvVariable("COUNT_GAPS_IN_FREQUENCIES", FALSE);
        __f = frequencies._aux.empirical.collect_data(datafilter, 3, 1, 1);
        utility.ToggleEnvVariable("COUNT_GAPS_IN_FREQUENCIES", None);
    } else {
        __f = model [^"terms.efv_estimate"];
    }


    __alphabet = model[^"terms.bases"];
    nuc_dict = {};
    utility.ForEachPair (__alphabet, "index","base",
    '
        `&nuc_dict`[base] = `&__f`[index[1]][-1];
    ');    
    
    frequencies.codon.multiply_in_frequencies  (model, nuc_dict);
    model[^"terms.model.efv_estimate_name"] = ^"terms.frequencies.F3x4";
    model[^"terms.efv_estimate"] = frequencies.empirical.codon_from_nuc  (model, nuc_dict);
    (model[^"terms.parameters"])[^"terms.model.empirical"] = 9;
    
    return model;
}     

/**
 * @name frequencies.empirical.F1x4
 * @param {Dictionary} model
 * @param {String} namespace
 * @param {DataSetFilter} datafilter
 * @returns {Dictionary} updated model
 */

lfunction frequencies.empirical.F1x4(model, namespace, datafilter) {

    frequencies._aux.validate (model);

    if (None != datafilter) {
        utility.ToggleEnvVariable("COUNT_GAPS_IN_FREQUENCIES", FALSE);
        __f1 = frequencies._aux.empirical.collect_data(datafilter, 1, 1, 1);
        __f = {4,3} ["__f1[_MATRIX_ELEMENT_ROW_][0]"];
        utility.ToggleEnvVariable("COUNT_GAPS_IN_FREQUENCIES", None);
    } else {
        __f = model [^"terms.efv_estimate"];
    }

    __alphabet = model[^"terms.bases"];
    nuc_dict = {};
    utility.ForEachPair (__alphabet, "index","base",
    '
        `&nuc_dict`[base] = `&__f`[index[1]][-1];
    ');    
    
    frequencies.codon.multiply_in_frequencies  (model, nuc_dict);
    model[^"terms.model.efv_estimate_name"] = ^"terms.frequencies.F1x4";
    model[^"terms.efv_estimate"] = frequencies.empirical.codon_from_nuc  (model, nuc_dict);
    (model[^"terms.parameters"])[^"terms.model.empirical"] = 3;
    return model;
}       


/**
 * @name frequencies.empirical.corrected.CF3x4
 * @param {Dictionary} model
 * @param {String} namespace
 * @param {DataSetFilter} datafilter
 * @returns {Dictionary} updated model
 */
lfunction frequencies.empirical.corrected.CF3x4(model, namespace, datafilter) {

    frequencies._aux.validate (model);
    
    if (None != datafilter) {
        utility.ToggleEnvVariable("COUNT_GAPS_IN_FREQUENCIES", 0);
        __f = frequencies._aux.empirical.collect_data(datafilter, 3, 1, 1);
        utility.ToggleEnvVariable("COUNT_GAPS_IN_FREQUENCIES", None);
    } else {
        __f = model [^"terms.efv_estimate"];
    }
    
    //TODO
    __alphabet = model[^"terms.alphabet"];
    __estimates = frequencies._aux.CF3x4(__f, model[^"terms.bases"], __alphabet, model[^"terms.stop_codons"]);
    
    model[^"terms.efv_estimate"] = __estimates[^"terms.codons"];
        
    frequencies.codon.multiply_in_frequencies  (model, __estimates[^"terms.bases"]);

    model[^"terms.model.efv_estimate_name"] = ^"terms.frequencies.CF3x4";
    (model[^"terms.parameters"])[^"terms.model.empirical"] = 9;
    return model;
}

/**
 * To be implemented
 * @name frequencies.mle
 * @param model
 * @param namespace
 * @param datafilter
 */
lfunction frequencies.mle(model, namespace, datafilter) {
    frequencies._aux.empirical.singlechar(model, namespace, datafilter); // gather empirical frequencies

    alphabet       = model[utility.getGlobalValue ("terms.alphabet")];
    alphabet_size  = utility.Array1D (alphabet);
    empirical       = model[utility.getGlobalValue ("terms.efv_estimate")];

    model[utility.getGlobalValue ("terms.efv_estimate")] = {alphabet_size,1};

    ML.variables = {alphabet_size,1};
    ML.scaler    = namespace + ".frequency_scaler";

    utility.ForEachPair (model[utility.getGlobalValue ("terms.alphabet")], "_index", "_letter",
                         '
                            _idx = _index[1];
                            _freq_parameter = `&namespace` + ".equilibrium_frequency_of." + _letter;
                            (`&ML.variables`) [_idx] = _freq_parameter;
                            (`&model`[terms.efv_estimate]) [_idx] = _freq_parameter + "/" + `&ML.scaler`;
                            parameters.DeclareGlobalWithRanges (_freq_parameter, `&empirical`[_idx], 0, 1);
                            model.generic.AddGlobal (`&model`, _freq_parameter, terms.characterFrequency (_letter));
                         '
                         );

    // constrain pi_A + ... + pi_A = 1,

    parameters.SetConstraint ( ML.scaler, Join ("+", ML.variables), "global");

    model[utility.getGlobalValue("terms.model.efv_estimate_name")] = utility.getGlobalValue ("terms.frequencies.MLE");
    (model[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.model.empirical")] = -1; // correct for the restricted sum

    return model;

}

//--- AUX FUNCTIONS FROM THIS POINT ON ---//

/**
 * Returns entire character count (# of sites * # of species)
 * @name frequencies._aux.empirical.character_count
 * @private
 * @param {DataSetFilter} datafilter - the datafilter containing site and species information
 */
function frequencies._aux.empirical.character_count(datafilter) {
    return Eval("`datafilter`.sites * `datafilter`.species");
}

/**
 * Collects frequency data from dataset
 * @name frequencies._aux.empirical.collect_data
 * @private
 * @param {DataSetFilter} datafilter
 * @param {Number} unit
 * @param {Number} stride
 * @param {Number} position_specific
 * @returns {Matrix} frequency data
 */
lfunction frequencies._aux.empirical.collect_data(datafilter, unit, stride, position_specific) {

    assert(Type(datafilter) == "Matrix" || Type(datafilter) == "AssociativeList" || Type(datafilter) == "String",
        "datafilter must be a String/Matrix/Associative List in call to frequencies._aux.empirical.collect_data, instead had a " + Type(datafilter)
    );
    if (Type(datafilter) == "String") {
        HarvestFrequencies(__f, ^ datafilter, unit, stride, position_specific);
    } else {
        site_count = 0;
        dim = utility.Array1D(datafilter);

        for (i = 0; i < dim; i += 1) {

            HarvestFrequencies(__f, ^ (datafilter[i]), unit, stride, position_specific);
            local_sites = frequencies._aux.empirical.character_count(datafilter[i]);

            if (i) {
                __f_composite += __f * local_sites;
            } else {
                __f_composite = __f * local_sites;

            }
            site_count += local_sites;
        }
        return __f_composite * (1 / site_count);
    }
    

    return __f;
}

/**
 * Updates model's equilibrium frequency estimator to not count gaps in frequencies
 * @name frequencies._aux.empirical.singlechar
 * @private
 * @param {Dictionary} model - the model to update the
 * @param {String} namespace - does nothing
 * @param {DataSetFilter} datafilter - the datasetfilter to collect data from
 * @returns {Dictionary} updated model
 */
function frequencies._aux.empirical.singlechar(model, namespace, datafilter) {
    utility.ToggleEnvVariable("COUNT_GAPS_IN_FREQUENCIES", 0);
    __f = frequencies._aux.empirical.collect_data(datafilter, 1, 1, 1);
    utility.ToggleEnvVariable("COUNT_GAPS_IN_FREQUENCIES", None);
    model[terms.efv_estimate] = __f;
    return model;
}

/**
 * @name frequencies._aux.CF3x4
 * @private
 * @param observed_3x4
 * @param base_alphabet
 * @param sense_codons
 * @param stop_codons
 * @returns  {Dictionary} codons and bases
 */
function frequencies._aux.CF3x4(observed_3x4, base_alphabet, sense_codons, stop_codons) {

    frequencies._aux.CF3x4.p = {};

    frequencies._aux.CF3x4.args = {};

    for (frequencies._aux.CF3x4.k = 0; frequencies._aux.CF3x4.k < 3; frequencies._aux.CF3x4.k += 1) {
        frequencies._aux.CF3x4.p[frequencies._aux.CF3x4.k] = parameters.GenerateSequentialNames("frequencies._aux.CF3x4.p" + frequencies._aux.CF3x4.k, 3, None);
        parameters.DeclareGlobal(frequencies._aux.CF3x4.p[frequencies._aux.CF3x4.k], None);
        parameters.SetRange(frequencies._aux.CF3x4.p[frequencies._aux.CF3x4.k], terms.range01);
        frequencies._aux.CF3x4.args + (Join(",", frequencies._aux.CF3x4.p[frequencies._aux.CF3x4.k]));
    }

    frequencies._aux.CF3x4.args = Join(",", frequencies._aux.CF3x4.args);

    frequencies._aux.CF3x4.n = {};

    for (frequencies._aux.CF3x4.k = 0; frequencies._aux.CF3x4.k < 3; frequencies._aux.CF3x4.k += 1) {
        frequencies._aux.CF3x4.n[frequencies._aux.CF3x4.k] = parameters.GenerateAttributedNames("frequencies._aux.CF3x4.n" + frequencies._aux.CF3x4.k, base_alphabet, None);
        parameters.SetConstraint(frequencies._aux.CF3x4.n[frequencies._aux.CF3x4.k],
            parameters.helper.stick_breaking(frequencies._aux.CF3x4.p[frequencies._aux.CF3x4.k], observed_3x4[-1][frequencies._aux.CF3x4.k]),
            utility.getGlobalValue("terms.global"));
    }


    frequencies._aux.stop_count = Columns(stop_codons);

    frequencies._aux.CF3x4.stop_correction = {};

    for (frequencies._aux.CF3x4.i = 0; frequencies._aux.CF3x4.i < Columns(base_alphabet); frequencies._aux.CF3x4.i += 1) {
        frequencies._aux.CF3x4.stop_correction[base_alphabet[frequencies._aux.CF3x4.i]] = {
            {
                "",
                "",
                ""
            }
        };
    }

    frequencies._aux.CF3x4.denominator = "1";

    for (frequencies._aux.CF3x4.i = 0; frequencies._aux.CF3x4.i < frequencies._aux.stop_count; frequencies._aux.CF3x4.i += 1) {
        frequencies._aux.CF3x4.sc = stop_codons[frequencies._aux.CF3x4.i];

        (frequencies._aux.CF3x4.stop_correction[frequencies._aux.CF3x4.sc[0]])[0] +=
        "-frequencies._aux.CF3x4.n1_" + frequencies._aux.CF3x4.sc[1] +
            "*frequencies._aux.CF3x4.n2_" + frequencies._aux.CF3x4.sc[2];

        (frequencies._aux.CF3x4.stop_correction[frequencies._aux.CF3x4.sc[1]])[1] +=
        "-frequencies._aux.CF3x4.n0_" + frequencies._aux.CF3x4.sc[0] +
            "*frequencies._aux.CF3x4.n2_" + frequencies._aux.CF3x4.sc[2];

        (frequencies._aux.CF3x4.stop_correction[frequencies._aux.CF3x4.sc[2]])[2] +=
        "-frequencies._aux.CF3x4.n0_" + frequencies._aux.CF3x4.sc[0] +
            "*frequencies._aux.CF3x4.n1_" + frequencies._aux.CF3x4.sc[1];

        frequencies._aux.CF3x4.denominator += "-frequencies._aux.CF3x4.n0_" + frequencies._aux.CF3x4.sc[0] +
            "*frequencies._aux.CF3x4.n1_" + frequencies._aux.CF3x4.sc[1] +
            "*frequencies._aux.CF3x4.n2_" + frequencies._aux.CF3x4.sc[2];
    }


    parameters.SetConstraint("frequencies._aux.CF3x4.denominator", frequencies._aux.CF3x4.denominator, utility.getGlobalValue("terms.global"));

    frequencies._aux.N = {
        Columns(base_alphabet),
        3
    };
    frequencies._aux.res = {};
    frequencies._aux.codons = {
        Columns(sense_codons),
        1
    };

    for (frequencies._aux.CF3x4.i = 0; frequencies._aux.CF3x4.i < Columns(sense_codons); frequencies._aux.CF3x4.i += 1) {
        frequencies._aux.CF3x4.sc = {
            3,
            1
        };
        for (frequencies._aux.CF3x4.pos = 0; frequencies._aux.CF3x4.pos < 3; frequencies._aux.CF3x4.pos += 1) {
            frequencies._aux.CF3x4.sc[frequencies._aux.CF3x4.pos] = "frequencies._aux.CF3x4.n" + frequencies._aux.CF3x4.pos + "_" + (sense_codons[frequencies._aux.CF3x4.i])[frequencies._aux.CF3x4.pos];
        }
        ExecuteCommands("frequencies._aux.codons[frequencies._aux.CF3x4.i] := " + Join("*", frequencies._aux.CF3x4.sc) + "/frequencies._aux.CF3x4.denominator");
    }

    for (frequencies._aux.CF3x4.i = 0; frequencies._aux.CF3x4.i < Columns(base_alphabet); frequencies._aux.CF3x4.i += 1) {
        frequencies._aux.CF3x4.n = base_alphabet[frequencies._aux.CF3x4.i];
        frequencies._aux.res[frequencies._aux.CF3x4.n] = {
            3,
            1
        };
        for (frequencies._aux.CF3x4.pos = 0; frequencies._aux.CF3x4.pos < 3; frequencies._aux.CF3x4.pos += 1) {

            frequencies._aux.CF3x4.sc = (frequencies._aux.CF3x4.stop_correction[frequencies._aux.CF3x4.n])[frequencies._aux.CF3x4.pos];

            if (Abs(frequencies._aux.CF3x4.sc)) {
                frequencies._aux.CF3x4.sc = "*(1" + frequencies._aux.CF3x4.sc + ")";
            }

            ExecuteCommands("frequencies._aux.N[" + frequencies._aux.CF3x4.i + "][" + frequencies._aux.CF3x4.pos + "] := frequencies._aux.CF3x4.n" + frequencies._aux.CF3x4.pos + "_" + frequencies._aux.CF3x4.n + frequencies._aux.CF3x4.sc + "/frequencies._aux.CF3x4.denominator");

            ExecuteCommands("(frequencies._aux.res[frequencies._aux.CF3x4.n])[frequencies._aux.CF3x4.pos] := frequencies._aux.CF3x4.n" + frequencies._aux.CF3x4.pos + "_" + frequencies._aux.CF3x4.n);
        }
    }



    ExecuteCommands("Optimize (frequencies._aux.CF3x4.p, frequencies._aux._CF3x4_minimizer( " +
        frequencies._aux.CF3x4.args + "))");

    return {
        utility.getGlobalValue("terms.codons"): Eval("frequencies._aux.codons"),
        utility.getGlobalValue("terms.bases"):  utility.Map (frequencies._aux.res,"_value_", "Eval(_value_)") // this is an in-place shallow copy
    };
}

/**
 * @name frequencies._aux._CF3x4_minimizer
 * @private
 * @param p11
 * @param p12
 * @param p13
 * @param p21
 * @param p22
 * @param p23
 * @param p31
 * @param p32
 * @param p33
 */
function frequencies._aux._CF3x4_minimizer(p11, p12, p13, p21, p22, p23, p31, p32, p33) {
    frequencies._aux._CF3x4_minimizer.error = frequencies._aux.N - observed_3x4;
    return -(+frequencies._aux._CF3x4_minimizer.error$frequencies._aux._CF3x4_minimizer.error);
}
