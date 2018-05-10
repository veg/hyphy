LoadFunctionLibrary("../protein.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");
LoadFunctionLibrary("../../all-terms.bf");

models.protein.empirical.default_generators = {"binary": "models.binary.ModelDescription"};
                             

models.binary.empirical.mleF_generators = {"binary": "models.binaryML.ModelDescription"};

/**
 * @name models.binary.ModelDescription
 * @description Create the baseline schema (dictionary) for the binary model of character evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.binaryML.ModelDescription(type) {
    models.binaryML.ModelDescription.model_definition = models.binary.ModelDescription(type);
    models.binaryML.ModelDescription.model_definition [terms.model.frequency_estimator] = "frequencies.ML.binary";
    models.binaryML.ModelDescription.model_definition [terms.model.efv_estimate_name]   =  utility.getGlobalValue("terms.frequencies.MLE");
    return models.binaryML.ModelDescription.model_definition;
}
