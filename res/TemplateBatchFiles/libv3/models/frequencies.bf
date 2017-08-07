LoadFunctionLibrary("terms.bf");
LoadFunctionLibrary("parameters.bf");
LoadFunctionLibrary("../UtilityFunctions.bf");
LoadFunctionLibrary("model_functions.bf");


/** @module frequencies */

/**
 * Sets model's equilibrium frequency estimator to equal frequencies
 * @name frequencies.equal
 * @param {Dictionary} model
 * @param {String} namespace - does nothing
 * @param {DataSetFilter} datafilter - does nothing
 * @returns {Dictionary} updated model
 */
function frequencies.equal(model, namespace, datafilter) {
    __N = Abs(model["alphabet"]);
    model[terms.efv_estimate] = {
        __N,
        1
    }["1/__N"];
    model[terms.efv_estimate_name] = terms.freqs.equal;
    (model["parameters"])["empirical"] = 0;
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
    model[terms.efv_estimate_name] = terms.freqs.4x1;
    (model["parameters"])["empirical"] = 3;
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
    model[terms.efv_estimate_name] = terms.freqs.20x1;
    (model["parameters"])["empirical"] = 19;
    return model;
}

/**
 * @name frequencies.empirical.corrected.CF3x4
 * @param {Dictionary} model
 * @param {String} namespace
 * @param {DataSetFilter} datafilter
 * @returns {Dictionary} updated model
 */
function frequencies.empirical.corrected.CF3x4(model, namespace, datafilter) {


    __dimension = model.Dimension(model);
    __alphabet = model["alphabet"];

    if (Type (model[terms.rate_matrix]) == "AssociativeList") {
        utility.ForEach (model[terms.rate_matrix], "_q_matrix_", '
            assert(Type(_q_matrix_) == "Matrix" && Rows(_q_matrix_) == __dimension && Columns(_q_matrix_) == __dimension,
            "`terms.rate_matrix` must be defined prior to calling frequencies.empirical.corrected.CF3x4")
        ');
    } else {
        assert(Type(model[terms.rate_matrix]) == "Matrix" && Rows(model[terms.rate_matrix]) == __dimension && Columns(model[terms.rate_matrix]) == __dimension,
            "`terms.rate_matrix` must be defined prior to calling frequencies.empirical.corrected.CF3x4");
    }

    utility.ToggleEnvVariable("COUNT_GAPS_IN_FREQUENCIES", 0);
    __f = frequencies._aux.empirical.collect_data(datafilter, 3, 1, 1);
    utility.ToggleEnvVariable("COUNT_GAPS_IN_FREQUENCIES", None);

    __estimates = frequencies._aux.CF3x4(__f, model["bases"], __alphabet, model["stop"]);
    model[terms.efv_estimate] = __estimates["codons"];
    __estimates = __estimates["bases"];

    if (Type (model[terms.rate_matrix]) == "AssociativeList") { // handle mixture models
        __components = Rows (model[terms.rate_matrix]);
        __component_count = utility.Array1D (__components);

        for (_rowChar = 0; _rowChar < __dimension; _rowChar += 1) {
            for (_colChar = _rowChar + 1; _colChar < __dimension; _colChar += 1) {

                __diff = models.codon.diff(__alphabet[_rowChar], __alphabet[_colChar]);
                if (None != __diff) {
                    for (__component_id = 0; __component_id < __component_count; __component_id += 1) {
                         ((model[terms.rate_matrix])[__components[__component_id]]) [_rowChar][_colChar] += "*" + (__estimates[__diff["to"]])[__diff["position"]];
                         ((model[terms.rate_matrix])[__components[__component_id]]) [_colChar][_rowChar] += "*" + (__estimates[__diff["from"]])[__diff["position"]];
                    }
                 }
            }
        }
    } else {
        for (_rowChar = 0; _rowChar < __dimension; _rowChar += 1) {
            for (_colChar = _rowChar + 1; _colChar < __dimension; _colChar += 1) {

                __diff = models.codon.diff(__alphabet[_rowChar], __alphabet[_colChar]);
                if (None != __diff) {
                    (model[terms.rate_matrix])[_rowChar][_colChar] += "*" + (__estimates[__diff["to"]])[__diff["position"]];
                    (model[terms.rate_matrix])[_colChar][_rowChar] += "*" + (__estimates[__diff["from"]])[__diff["position"]];
                }
            }
        }
    }


    model[terms.efv_estimate_name] = terms.freqs.CF3x4;
    (model["parameters"])["empirical"] = 9;
    return model;
}

/**
 * To be implemented
 * @name frequencies.mle
 * @param model
 * @param namespace
 * @param datafilter
 */
function frequencies.mle(model, namespace, datafilter) {
    assert(0, "frequencies.mle is TBD");
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
            "global");
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


    parameters.SetConstraint("frequencies._aux.CF3x4.denominator", frequencies._aux.CF3x4.denominator, "global");

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
        "codons": Eval("frequencies._aux.codons"),
        "bases":  utility.Map (frequencies._aux.res,"_value_", "Eval(_value_)") // this is an in-place shallow copy
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
