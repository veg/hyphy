LoadFunctionLibrary("parameters.bf");

/** @module rate_variation */

rate_variation.types = {
    "Gamma+I": "rate_variation.types.GammaI"
};

/**
 * @name rate_variation.Add 
 * @param model
 * @param specs
 */
function rate_variation.Add(model, specs) {

    rate_variation.Add.spec = rate_variation.types[specs["type"]];
    assert(Type(rate_variation.Add.spec) == "String",
        "Unsupported rate variation type '`type`' in call to rate_variation. Use one of " + Rows(rate_variation.types));

    model[terms.rate_variation.bins] = specs["bins"];
    model["original q call-back"] = model["defineQ"];
    model["defineQ"] = rate_variation.Add.spec;
    return model;
}

/**
 * to be implemented
 * @name rate_variation.MultiplyIn
 * @param matrix
 * @param variable
 */
function rate_variation.MultiplyIn(matrix, variable) {

}

/**
 * @name rate_variation.types.GammaI
 * @param model
 * @param namespace
 * @returns rate_variation.types.GammaI.q
 */
function rate_variation.types.GammaI(model, namespace) {

    __global_cache = {};



    rate_variation.types.GammaI.alpha = parameters.ApplyNameSpace("alpha", namespace);
    parameters.DeclareGlobal(rate_variation.types.GammaI.alpha, __global_cache);

    rate_variation.types.GammaI.range = {};
    rate_variation.types.GammaI.range[terms.lower_bound] = 0.001;
    rate_variation.types.GammaI.range[terms.upper_bound] = 100;

    parameters.SetRange(rate_variation.types.GammaI.alpha,
        rate_variation.types.GammaI.range);


    rate_variation.types.GammaI.inv_p = parameters.ApplyNameSpace("inv_p", namespace);
    parameters.DeclareGlobal(rate_variation.types.GammaI.inv_p, __global_cache);

    parameters.SetValue(rate_variation.types.GammaI.inv_p, 0.1);

    parameters.SetRange(rate_variation.types.GammaI.inv_p, terms.range01);

    rate_variation.types.GammaI.beta = parameters.ApplyNameSpace("beta", namespace);
    parameters.SetConstraint(rate_variation.types.GammaI.beta,
        "(1-`rate_variation.types.GammaI.inv_p`)*`rate_variation.types.GammaI.alpha`",
        "global");

    rate_variation.types.GammaI.range[terms.upper_bound] = 200;
    parameters.SetRange(rate_variation.types.GammaI.beta,
        rate_variation.types.GammaI.range);

    ((model["parameters"])["global"])[terms.rate_variation.gamma_alpha] = rate_variation.types.GammaI.alpha;
    ((model["parameters"])["global"])[terms.rate_variation.gamma_beta] = rate_variation.types.GammaI.beta;
    ((model["parameters"])["global"])[terms.rate_variation.gamma_p_inv] = rate_variation.types.GammaI.inv_p;

    rate_variation.types.GammaI.bins = 0 + model[terms.rate_variation.bins];


    rate_variation.types.GammaI.main = parameters.ApplyNameSpace("gamma_i", namespace);
    rate_variation.types.GammaI.cats = parameters.ApplyNameSpace("gamma_i.freqs", namespace);

    ExecuteCommands(rate_variation.types.GammaI.cats + " = {1 + rate_variation.types.GammaI.bins, 1}");
    ExecuteCommands(rate_variation.types.GammaI.cats + "[0] := `rate_variation.types.GammaI.inv_p`");

    for (rate_variation.types.GammaI.i = 1; rate_variation.types.GammaI.i <= rate_variation.types.GammaI.bins; rate_variation.types.GammaI.i += 1) {
        ExecuteCommands(rate_variation.types.GammaI.cats +
            "[rate_variation.types.GammaI.i] := (1-`rate_variation.types.GammaI.inv_p`)/" +
            rate_variation.types.GammaI.bins);

    }


    ExecuteCommands("category `rate_variation.types.GammaI.main` = (" +
        (1 + rate_variation.types.GammaI.bins) + "," +
        (rate_variation.types.GammaI.cats) +
        ",MEAN, 
        (1 - `rate_variation.types.GammaI.inv_p`) * GammaDist(_x_, `rate_variation.types.GammaI.alpha`, `rate_variation.types.GammaI.beta`) * (_x_ > 0),
        (1 - `rate_variation.types.GammaI.inv_p`) * CGammaDist(_x_, `rate_variation.types.GammaI.alpha`, `rate_variation.types.GammaI.beta`) * (_x_ > 0) + `rate_variation.types.GammaI.inv_p`,
        0,
        1e25,
        (1 - `rate_variation.types.GammaI.inv_p`) * CGammaDist(_x_, `rate_variation.types.GammaI.alpha` + 1, `rate_variation.types.GammaI.beta`) * (`rate_variation.types.GammaI.alpha` / `rate_variation.types.GammaI.beta`) * (_x_ > 0)
    )
    ");

    ExecuteCommands("GetInformation (info, `rate_variation.types.GammaI.main`)");

    rate_variation.types.GammaI.q = Call(model["original q call-back"],
        model,
        namespace);

    rate_variation.types.GammaI.q[terms.rate_matrix] = parameters.AddMultiplicativeTerm(rate_variation.types.GammaI.q[terms.rate_matrix], rate_variation.types.GammaI.main, 0);
    return rate_variation.types.GammaI.q;
}
