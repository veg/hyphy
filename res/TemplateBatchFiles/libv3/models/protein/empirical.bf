LoadFunctionLibrary("../protein.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");

/** @module models.protein.empirical */

/**
 * @name models.protein.empirical.ModelDescription
 * @param {String} type
 */

function models.protein.empirical.ModelDescription(type) {
    return {
        "alphabet": models.protein.alphabet,
        "description": "General class of empirical substitution matrices for amino-acids",
        "canonical": 1, // is of the r_ij \times \pi_j form
        "reversible": 1,
        terms.efv_estimate_name: terms.freqs.predefined,
        "parameters": {
            "global": {},
            "local": {},
            "empirical": 0
        },
        "type": type,
        "get-branch-length": "",
        "set-branch-length": "models.generic.SetBranchLength",
        "constrain-branch-length": "models.generic.constrain_branch_length",
        "frequency-estimator": "frequencies.empirical.protein",
        "q_ij": "",
        "time": "models.protein.generic.Time",
        "defineQ": "models.protein.empirical._DefineQ",
        "post-definition": "models.generic.post.definition"
    };
}

function models.protein.WAG._GenerateRate (from,to,namespace,modelType) {
    return models.protein.empirical._GenerateRate (models.protein.WAG.empirical_Q, from,to,namespace,modelType);
}

function models.protein.WAG.ModelDescription(type) {
    models.protein.WAG.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.WAG.ModelDescription.model_definition ["empirical-rates"] = models.protein.WAG.empirical_Q;
    models.protein.WAG.ModelDescription.model_definition ["frequency-estimator"] = "models.protein.WAG.frequencies";
    models.protein.WAG.ModelDescription.model_definition ["q_ij"] = "models.protein.WAG._GenerateRate";
    return models.protein.WAG.ModelDescription.model_definition;
}

function models.protein.WAG.frequencies (model, namespace, datafilter) {
    model[terms.efv_estimate] =
        {{     0.0866279}
        {     0.0193078}
        {     0.0570451}
        {     0.0580589}
        {     0.0384319}
        {     0.0832518}
        {     0.0244313}
        {      0.048466}
        {     0.0620286}
        {      0.086209}
        {     0.0195027}
        {     0.0390894}
        {     0.0457631}
        {     0.0367281}
        {      0.043972}
        {     0.0695179}
        {     0.0610127}
        {     0.0708956}
        {     0.0143859}
        {     0.0352742}
        };

    model[terms.efv_estimate_name] = terms.freqs.predefined;
    (model["parameters"])["empirical"] = 0;
    return model;
}

models.protein.WAG.empirical_Q = {
 "A":{
   "C":3.26324,
   "D":2.34804,
   "E":5.02923,
   "F":0.668808,
   "G":4.50138,
   "H":1.00707,
   "I":0.6142879999999999,
   "K":2.8795,
   "L":1.26431,
   "M":2.83893,
   "N":1.61995,
   "P":4.57074,
   "Q":2.88691,
   "R":1.75252,
   "S":10.7101,
   "T":6.73946,
   "V":6.37375,
   "W":0.35946,
   "Y":0.764894
  },
 "C":{
   "D":0.0962568,
   "E":0.06784229999999999,
   "F":1.26464,
   "G":0.974403,
   "H":0.791065,
   "I":0.540574,
   "K":0.23523,
   "L":1.22101,
   "M":1.24069,
   "N":0.842805,
   "P":0.347612,
   "Q":0.313977,
   "R":1.67824,
   "S":4.4726,
   "T":1.62992,
   "V":3.18413,
   "W":2.27837,
   "Y":1.72794
  },
 "D":{
   "E":19.6173,
   "F":0.148478,
   "G":2.75024,
   "H":2.95706,
   "I":0.125304,
   "K":1.52466,
   "L":0.269452,
   "M":0.32966,
   "N":17.251,
   "P":1.34714,
   "Q":1.95972,
   "R":0.468033,
   "S":3.40533,
   "T":1.19107,
   "V":0.484018,
   "W":0.412312,
   "Y":1.03489
  },
 "E":{
   "F":0.257789,
   "G":1.80382,
   "H":1.81116,
   "I":0.404776,
   "K":8.21158,
   "L":0.490144,
   "M":1.00125,
   "N":3.00956,
   "P":2.16806,
   "Q":17.3783,
   "R":1.39535,
   "S":2.23982,
   "T":2.61419,
   "V":1.87059,
   "W":0.497433,
   "Y":0.623719
  },
 "F":{
   "G":0.158647,
   "H":2.15858,
   "I":3.36628,
   "K":0.282261,
   "L":6.72059,
   "M":3.78302,
   "N":0.305538,
   "P":0.51296,
   "Q":0.317481,
   "R":0.326346,
   "S":1.7346,
   "T":0.546192,
   "V":2.06492,
   "W":4.86017,
   "Y":20.5074
  },
 "G":{
   "H":0.792457,
   "I":0.0967499,
   "K":1.18692,
   "L":0.194782,
   "M":0.553173,
   "N":3.57627,
   "P":0.773901,
   "Q":1.04868,
   "R":1.85767,
   "S":4.2634,
   "T":0.717545,
   "V":0.5949449999999999,
   "W":1.07071,
   "Y":0.329184
  },
 "H":{
   "I":0.439075,
   "K":2.82919,
   "L":1.58695,
   "M":1.28409,
   "N":12.5704,
   "P":2.21205,
   "Q":13.6438,
   "R":6.79042,
   "S":2.35176,
   "T":1.50385,
   "V":0.376062,
   "W":0.834267,
   "Y":12.3072
  },
 "I":{
   "K":1.02892,
   "L":10.0752,
   "M":13.5273,
   "N":1.76099,
   "P":0.317506,
   "Q":0.361952,
   "R":0.594093,
   "S":1.01497,
   "T":4.63305,
   "V":24.8508,
   "W":0.675128,
   "Y":1.33502
  },
 "K":{
   "L":0.818336,
   "M":2.9685,
   "N":9.57014,
   "P":1.76944,
   "Q":12.3754,
   "R":17.0032,
   "S":3.07289,
   "T":4.40689,
   "V":0.970464,
   "W":0.436898,
   "Y":0.423423
  },
 "L":{
   "M":15.4228,
   "N":0.417907,
   "P":1.32127,
   "Q":2.76265,
   "R":1.58126,
   "S":1.09535,
   "T":1.03778,
   "V":5.72027,
   "W":2.1139,
   "Y":1.26654
  },
 "M":{
   "N":0.629813,
   "P":0.544368,
   "Q":4.9098,
   "R":2.17063,
   "S":1.5693,
   "T":4.81721,
   "V":6.54037,
   "W":1.63857,
   "Y":1.36128
  },
 "N":{
   "P":0.6198360000000001,
   "Q":4.90465,
   "R":2.0187,
   "S":12.6274,
   "T":6.45016,
   "V":0.623538,
   "W":0.228503,
   "Y":3.45058
  },
 "P":{
   "Q":2.96563,
   "R":2.15896,
   "S":5.12592,
   "T":2.52719,
   "V":1.0005,
   "W":0.442935,
   "Y":0.686449
  },
 "Q":{
   "R":9.644769999999999,
   "S":3.26906,
   "T":2.72592,
   "V":0.957268,
   "W":0.685467,
   "Y":0.723509
  },
 "R":{
   "S":3.88965,
   "T":1.76155,
   "V":0.800207,
   "W":3.69815,
   "Y":1.21225
  },
 "S":{
   "T":13.9104,
   "V":0.739488,
   "W":1.6641,
   "Y":2.50053
  },
 "T":{
   "V":4.41086,
   "W":0.352251,
   "Y":0.925072
  },
 "V":{
   "W":1.1609,
   "Y":1
  },
 "W":{
   "Y":7.8969
  }
};

/**
 * @name models.protein.Empirical._GenerateRate
 * return the r_ij component of
 * q_ij := r_ij * time * freq_j
 * @param {Dict} rateDict
 * @param {Number} fromChar
 * @param {Number} toChar
 * @param {String} namespace
 * @param {String} model_type
 * @return list of parameters
 */


function models.protein.empirical._GenerateRate(rateDict, fromChar, toChar, namespace, model_type) {

    models.protein.empirical._GenerateRate.p = {};
    models.protein.empirical._GenerateRate.p  [model_type]       = {};
    if (fromChar < toChar) {
        models.protein.empirical._GenerateRate.p  [terms.rate_entry] = "" + (rateDict[fromChar])[toChar];
    } else {
        models.protein.empirical._GenerateRate.p  [terms.rate_entry] = "" + (rateDict[toChar])[fromChar];
    }
    return models.protein.empirical._GenerateRate.p;
}

/**
 * @name models.protein.empirical._DefineQ
 * @param {Dictionary} model definition
 * @param {String} namespace
 */
function models.protein.empirical._DefineQ(model_dict, namespace) {
    models.protein.generic.DefineQMatrix (model_dict, namespace);
    return model_dict;
}
