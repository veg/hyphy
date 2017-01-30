LoadFunctionLibrary("../protein.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");
LoadFunctionLibrary("../terms.bf");

/** @module models.protein.empirical */

/**
 * @name models.protein.empirical.ModelDescription
 * @param {String} type
 * @returns {Dictionary} model description
 * @description Create the baseline schema (dictionary) for empirical protein model definitions
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


/**
 * @name models.protein.WAG._GenerateRate
 * @private
 * @description Generate rate matrix for the WAG model of protein evolution (http://mbe.oxfordjournals.org/content/18/5/691)
 * @param {String} from - source amino acid
 * @param {String} to   - target amino acid
 * @param {String} namespace
 * @param {String} modelType
 */
function models.protein.WAG._GenerateRate (from,to,namespace,modelType) {
    return models.protein.empirical._GenerateRate (models.protein.WAG.empirical_Q, from,to,namespace,modelType);
}

/**
 * @name models.protein.WAG.ModelDescription
 * @description Create the baseline schema (dictionary) for the WAG model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.WAG.ModelDescription(type) {
    models.protein.WAG.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.WAG.ModelDescription.model_definition ["empirical-rates"] = models.protein.WAG.empirical_Q;
    models.protein.WAG.ModelDescription.model_definition ["frequency-estimator"] = "models.protein.WAG.frequencies";
    models.protein.WAG.ModelDescription.model_definition ["q_ij"] = "models.protein.WAG._GenerateRate";
    return models.protein.WAG.ModelDescription.model_definition;
}

/**
 * @name models.protein.WAGF.ModelDescription
 * @description Create the baseline schema (dictionary) for the WAG+F model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.WAGF.ModelDescription(type) {
    models.protein.WAGF.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.WAGF.ModelDescription.model_definition ["empirical-rates"] = models.protein.WAG.empirical_Q;
    models.protein.WAGF.ModelDescription.model_definition ["frequency-estimator"] = "frequencies.empirical.protein";
    models.protein.WAGF.ModelDescription.model_definition ["parameters"] = {"global": {}, "local": {}, "empirical": 19};
    models.protein.WAGF.ModelDescription.model_definition ["q_ij"] = "models.protein.WAG._GenerateRate";
    return models.protein.WAGF.ModelDescription.model_definition;
}

/**
 * @name models.protein.WAG.frequencies
 * @param {Dictionary} Baseline WAG model
 * @returns {Dictionary} Updated WAG model with empirical frequencies
 * @description Define the empirical amino acid frequencies associated with the WAG model of protein evolution
 */
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


/* Define a dictionary of amino-acid exchangeability rates for the WAG model of protein evolution. */ 
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
 * @name models.protein.LG._GenerateRate
 * @private
 * @description Generate rate matrix for the LG model of protein evolution (doi: 10.1093/molbev/msn067)
 * @param {String} from - source amino acid
 * @param {String} to   - target amino acid
 * @param {String} namespace
 * @param {String} modelType
 */
function models.protein.LG._GenerateRate (from,to,namespace,modelType) {
    return models.protein.empirical._GenerateRate (models.protein.LG.empirical_Q, from,to,namespace,modelType);
}

/**
 * @name models.protein.LG.ModelDescription
 * @description Create the baseline schema (dictionary) for the LG model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
 function models.protein.LG.ModelDescription(type) {
    models.protein.LG.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.LG.ModelDescription.model_definition ["empirical-rates"] = models.protein.LG.empirical_Q;
    models.protein.LG.ModelDescription.model_definition ["frequency-estimator"] = "models.protein.LG.frequencies";
    models.protein.LG.ModelDescription.model_definition ["q_ij"] = "models.protein.LG._GenerateRate";
    return models.protein.LG.ModelDescription.model_definition;
}

/**
 * @name models.protein.LGF.ModelDescription
 * @description Create the baseline schema (dictionary) for the LG+F model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.LGF.ModelDescription(type) {
    models.protein.LGF.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.LGF.ModelDescription.model_definition ["empirical-rates"] = models.protein.LG.empirical_Q;
    models.protein.LGF.ModelDescription.model_definition ["frequency-estimator"] = "frequencies.empirical.protein";
    models.protein.LGF.ModelDescription.model_definition ["parameters"] = {"global": {}, "local": {}, "empirical": 19};
    models.protein.LGF.ModelDescription.model_definition ["q_ij"] = "models.protein.LG._GenerateRate";
    return models.protein.LGF.ModelDescription.model_definition;
}


/**
 * @name models.protein.LG.frequencies
 * @param {Dictionary} Baseline LG model
 * @returns {Dictionary} Updated LG model with empirical frequencies
 * @description Define the empirical amino acid frequencies associated with the LG model of protein evolution
 */
function models.protein.LG.frequencies (model, namespace, datafilter) {
    model[terms.efv_estimate] =
 
          {{     0.07906500000000008}
           {     0.012937}
           {     0.053052}
           {     0.071586}
           {     0.042302}
           {     0.057337}
           {     0.022355}
           {     0.062157}
           {     0.0646}
           {     0.099081}
           {     0.022951}
           {     0.041977}
           {     0.04404}
           {     0.040767}
           {     0.055941}
           {     0.061197}
           {     0.053287}
           {     0.069147}
           {     0.012066}
           {     0.034155}
          };
    model[terms.efv_estimate_name] = terms.freqs.predefined;
    (model["parameters"])["empirical"] = 0;
    return model;
}

/* Define a dictionary of amino-acid exchangeability rates for the LG model of protein evolution. */ 
models.protein.LG.empirical_Q = {
    "A":
        {"C":7.90863,
        "D":1.255501,
        "E":3.299796,
        "F":0.806091,
        "G":6.564482,
        "H":1.140209,
        "I":0.476059,
        "K":1.704692,
        "L":1.256114,
        "M":3.571425,
        "N":0.879541,
        "P":3.741781,
        "Q":3.081669,
        "R":1.350659,
        "S":15.019796,
        "T":6.797891,
        "V":8.095413000000001,
        "W":0.574197,
        "Y":0.695704},
    "C":
        {"D":0.198761,
        "E":0.011117,
        "F":3.511742,
        "G":1.80874,
        "H":2.035214,
        "I":1.018736,
        "K":0.04215,
        "L":1.887354,
        "M":2.839512,
        "N":1.680068,
        "P":0.239513,
        "Q":0.269463,
        "R":1.698443,
        "S":8.847193000000001,
        "T":3.633208,
        "V":6.225305,
        "W":2.129215,
        "Y":3.703275},
    "D":
        {"E":16.661482,
        "F":0.055336,
        "G":2.684605,
        "H":2.945743,
        "I":0.033966,
        "K":0.899053,
        "L":0.047901,
        "M":0.081174,
        "N":16.128578,
        "P":1.253315,
        "Q":1.662968,
        "R":0.393842,
        "S":3.940757,
        "T":1.353096,
        "V":0.120634,
        "W":0.09497,
        "Y":0.429279},
    "E":
        {"F":0.059769,
        "G":1.1084,
        "H":1.346808,
        "I":0.140644,
        "K":5.74199,
        "L":0.221374,
        "M":0.552013,
        "N":1.721195,
        "P":1.332599,
        "Q":13.117878,
        "R":1.156451,
        "S":1.944437,
        "T":1.920836,
        "V":0.7785530000000001,
        "W":0.247361,
        "Y":0.381397},
    "F":
        {"G":0.284644,
        "H":2.167378,
        "I":3.535496,
        "K":0.07599499999999999,
        "L":8.237826999999999,
        "M":5.715542,
        "N":0.28445,
        "P":0.300143,
        "Q":0.113923,
        "R":0.167515,
        "S":1.149617,
        "T":0.524262,
        "V":2.080141,
        "W":7.807073,
        "Y":24.795538},
    "G":
        {"H":0.989686,
        "I":0.027659,
        "K":0.942509,
        "L":0.140632,
        "M":0.443358,
        "N":4.567866,
        "P":0.6258089999999999,
        "Q":0.851393,
        "R":1.239767,
        "S":5.528515,
        "T":0.412531,
        "V":0.243704,
        "W":0.853083,
        "Y":0.173733},
    "H":
        {"I":0.345954,
        "K":2.215435,
        "L":1.163908,
        "M":1.405878,
        "N":14.327317,
        "P":1.616785,
        "Q":15.294073,
        "R":7.710101,
        "S":3.14559,
        "T":1.856391,
        "V":0.378143,
        "W":1.897035,
        "Y":16.861539},
    "I":
        {"K":0.505414,
        "L":13.170227,
        "M":13.578641,
        "N":0.608467,
        "P":0.248724,
        "Q":0.231481,
        "R":0.403492,
        "S":0.203682,
        "T":3.284525,
        "V":2.062424,
        "W":0.35478,
        "Y":0.738801},
    "K":
        {"L":0.436882,
        "M":2.086245,
        "N":6.815611,
        "P":1.24018,
        "Q":10.276405,
        "R":20.099975,
        "S":2.37881,
        "T":3.612184,
        "V":0.5884470000000001,
        "W":0.158568,
        "Y":0.419191},
    "L":
        {"M":20.056417,
        "N":0.217415,
        "P":0.791345,
        "Q":1.850656,
        "R":0.9590689999999999,
        "S":0.5791849999999999,
        "T":0.962526,
        "V":5.410175,
        "W":1.968773,
        "Y":0.952079},
    "M":
        {"N":1.178801,
        "P":0.317253,
        "Q":5.314296,
        "R":1.538248,
        "S":1.102405,
        "T":6.419361,
        "V":6.032845,
        "W":2.211975,
        "Y":1.529266},
    "N":
        {"P":0.51405,
        "Q":5.387956,
        "R":2.388961,
        "S":12.735858,
        "T":6.356809,
        "V":0.265904,
        "W":0.144174,
        "Y":1.944603},
    "P":
        {"Q":1.983585,
        "R":1.056566,
        "S":4.251681,
        "T":1.81574,
        "V":0.94208,
        "W":0.302262,
        "Y":0.28473},
    "Q":
        {"R":8.921638,
        "S":3.8885,
        "T":3.431944,
        "V":0.668293,
        "W":0.750481,
        "Y":0.81764},
    "R":
        {"S":2.726625,
        "T":1.83963,
        "V":0.542964,
        "W":1.886083,
        "Y":0.999078},
    "S":
        {"T":20.564538,
        "V":0.31255,
        "W":0.790716,
        "Y":1.272668},
    "T":
        {"V":6.95249,
        "W":0.447447,
        "Y":0.781117},
    "V":
        {"W":0.602135,
        "Y":0.792149},
    "W":
        {"Y":10.014342}
};
//==================================================================================================================//



/**
 * @name models.protein.JTT._GenerateRate
 * @private
 * @description Generate rate matrix for the JTT model of protein evolution (doi: 10.1093/bioinformatics/8.3.275)
 * @param {String} from - source amino acid
 * @param {String} to   - target amino acid
 * @param {String} namespace
 * @param {String} modelType
 */
function models.protein.JTT._GenerateRate (from,to,namespace,modelType) {
    return models.protein.empirical._GenerateRate (models.protein.JTT.empirical_Q, from,to,namespace,modelType);
}

/**
 * @name models.protein.JTT.ModelDescription
 * @description Create the baseline schema (dictionary) for the JTT model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
 function models.protein.JTT.ModelDescription(type) {
    models.protein.JTT.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.JTT.ModelDescription.model_definition ["empirical-rates"] = models.protein.JTT.empirical_Q;
    models.protein.JTT.ModelDescription.model_definition ["frequency-estimator"] = "models.protein.JTT.frequencies";
    models.protein.JTT.ModelDescription.model_definition ["q_ij"] = "models.protein.JTT._GenerateRate";
    return models.protein.JTT.ModelDescription.model_definition;
}

/**
 * @name models.protein.JTTF.ModelDescription
 * @description Create the baseline schema (dictionary) for the JTT+F model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.JTTF.ModelDescription(type) {
    models.protein.JTTF.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.JTTF.ModelDescription.model_definition ["empirical-rates"] = models.protein.JTT.empirical_Q;
    models.protein.JTTF.ModelDescription.model_definition ["frequency-estimator"] = "frequencies.empirical.protein";
    models.protein.JTTF.ModelDescription.model_definition ["parameters"] = {"global": {}, "local": {}, "empirical": 19};
    models.protein.JTTF.ModelDescription.model_definition ["q_ij"] = "models.protein.JTT._GenerateRate";
    return models.protein.JTTF.ModelDescription.model_definition;
}


/**
 * @name models.protein.JTT.frequencies
 * @param {Dictionary} Baseline JTT model
 * @returns {Dictionary} Updated JTT model with empirical frequencies
 * @description Define the empirical amino acid frequencies associated with the JTT model of protein evolution
 */
function models.protein.JTT.frequencies (model, namespace, datafilter) {
    model[terms.efv_estimate] =
      {{	0.07686099999999986}
        {	0.020279}
        {	0.051269}
        {	0.06182}
        {	0.04053}
        {	0.074714}
        {	0.022983}
        {	0.052569}
        {	0.059498}
        {	0.091111}
        {	0.023414}
        {	0.042546}
        {	0.050532}
        {	0.041061}
        {	0.051057}
        {	0.068225}
        {	0.058518}
        {	0.066374}
        {	0.014336}
        {   0.032303}
    };
    model[terms.efv_estimate_name] = terms.freqs.predefined;
    (model["parameters"])["empirical"] = 0;
    return model;
}

 /* Define a dictionary of amino-acid exchangeability rates for the JTT model of protein evolution. */ 
models.protein.JTT.empirical_Q = {
    "A":
        {"C":1.825304,
        "D":2.629062,
        "E":3.389193,
        "F":0.439402,
        "G":5.529052,
        "H":0.698916,
        "I":1.149188,
        "K":1.173822,
        "L":0.984993,
        "M":1.491421,
        "N":1.772843,
        "P":6.226284,
        "Q":1.768897,
        "R":1.689314,
        "S":12.350566,
        "T":14.560301,
        "V":9.291012,
        "W":0.267941,
        "Y":0.443212},
    "C":
        {"D":0.335605,
        "E":0.17128,
        "F":2.155291,
        "G":1.736056,
        "H":2.303555,
        "I":0.478375,
        "K":0.155718,
        "L":0.522966,
        "M":1.300168,
        "N":0.995491,
        "P":0.392886,
        "Q":0.290103,
        "R":3.240373,
        "S":6.848188,
        "T":1.492781,
        "V":1.974145,
        "W":3.508343,
        "Y":6.719573},
    "D":
        {"E":24.67688,
        "F":0.103333,
        "G":4.042937,
        "H":3.280087,
        "I":0.368468,
        "K":0.897486,
        "L":0.195361,
        "M":0.603695,
        "N":17.632664,
        "P":0.404041,
        "Q":1.657439,
        "R":0.492165,
        "S":1.872296,
        "T":1.350869,
        "V":1.001687,
        "W":0.182588,
        "Y":1.442353},
    "E":
        {"F":0.139259,
        "G":3.544726,
        "H":0.7745300000000001,
        "I":0.355139,
        "K":5.502124,
        "L":0.309742,
        "M":0.556299,
        "N":1.83686,
        "P":0.610027,
        "Q":10.859165,
        "R":1.011924,
        "S":0.992752,
        "T":1.05355,
        "V":1.478317,
        "W":0.363426,
        "Y":0.201608},
    "F":
        {"G":0.15954,
        "H":1.440688,
        "I":2.469068,
        "K":0.07791099999999999,
        "L":7.944248,
        "M":1.385889,
        "N":0.233473,
        "P":0.471779,
        "Q":0.14515,
        "R":0.207524,
        "S":2.999303,
        "T":0.441343,
        "V":1.885673,
        "W":1.709153,
        "Y":17.425203},
    "G":
        {"H":0.640854,
        "I":0.170842,
        "K":0.85737,
        "L":0.220799,
        "M":0.414257,
        "N":2.457067,
        "P":0.661141,
        "Q":0.734896,
        "R":4.320057,
        "S":5.955249,
        "T":1.006774,
        "V":1.493788,
        "W":1.729037,
        "Y":0.16681},
    "H":
        {"I":0.5776,
        "K":1.668401,
        "L":1.71757,
        "M":1.047437,
        "N":12.791207,
        "P":3.628382,
        "Q":18.060173,
        "R":10.201347,
        "S":2.362208,
        "T":1.516712,
        "V":0.387084,
        "W":0.407311,
        "Y":18.582271},
    "I":
        {"K":0.643605,
        "L":7.419497,
        "M":15.351776,
        "N":1.560076,
        "P":0.313221,
        "Q":0.248689,
        "R":0.76,
        "S":1.287195,
        "T":8.114273000000001,
        "V":30.292441,
        "W":0.427382,
        "Y":0.964144},
    "K":
        {"L":0.465418,
        "M":1.984497,
        "N":8.037099,
        "P":0.687399,
        "Q":9.426273,
        "R":20.745569,
        "S":1.507571,
        "T":3.068156,
        "V":0.394198,
        "W":0.283208,
        "Y":0.2793},
    "L":
        {"M":12.254646,
        "N":0.436212,
        "P":3.369566,
        "Q":2.252737,
        "R":1.182794,
        "S":1.8826,
        "T":0.865866,
        "V":5.596665,
        "W":1.685012,
        "Y":0.766034},
    "M":
        {"N":1.050805,
        "P":0.521765,
        "Q":1.451723,
        "R":1.36957,
        "S":0.90733,
        "T":6.719179,
        "V":9.654408999999999,
        "W":0.6397040000000001,
        "Y":0.603279},
    "N":
        {"P":0.387011,
        "Q":2.442836,
        "R":1.433276,
        "S":16.070798,
        "T":7.470881,
        "V":0.522749,
        "W":0.08801200000000001,
        "Y":2.22633},
    "P":
        {"Q":5.10954,
        "R":2.257455,
        "S":8.859674,
        "T":3.739588,
        "V":0.672198,
        "W":0.222302,
        "Y":0.361739},
    "Q":
        {"R":9.601862000000001,
        "S":1.743739,
        "T":1.664363,
        "V":0.571191,
        "W":0.5471549999999999,
        "Y":0.809408},
    "R":
        {"S":3.182254,
        "T":2.066158,
        "V":0.546484,
        "W":3.996952,
        "Y":0.7485810000000001},
    "S":
        {"T":15.18014,
        "V":1.298039,
        "W":0.987916,
        "Y":1.997292},
    "T":
        {"V":3.634797,
        "W":0.255953,
        "Y":0.638941},
    "V":
        {"W":0.761595,
        "Y":0.525762},
    "W":
        {"Y":2.376287}
};

//==================================================================================================================//


/**
 * @name models.protein.JC._GenerateRate
 * @private
 * @description Generate rate matrix for the JC69 (equal rates) model of protein evolution
 * @param {String} from - source amino acid
 * @param {String} to   - target amino acid
 * @param {String} namespace
 * @param {String} modelType
 */
function models.protein.JC._GenerateRate (from,to,namespace,modelType) {
    return models.protein.empirical._GenerateRate (models.protein.JC.empirical_Q, from,to,namespace,modelType);
}

 /**
 * @name models.protein.JC.ModelDescription
 * @description Create the baseline schema (dictionary) for the JC69 (equal rates) model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.JC.ModelDescription(type) {
    models.protein.JC.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.JC.ModelDescription.model_definition ["empirical-rates"] = models.protein.JC.empirical_Q;
    models.protein.JC.ModelDescription.model_definition ["frequency-estimator"] = "models.protein.JC.frequencies";
    models.protein.JC.ModelDescription.model_definition ["q_ij"] = "models.protein.JC._GenerateRate";
    return models.protein.JC.ModelDescription.model_definition;
}

/**
 * @name models.protein.JCF.ModelDescription
 * @description Create the baseline schema (dictionary) for the JC69+F model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
function models.protein.JCF.ModelDescription(type) {
    models.protein.JCF.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.JCF.ModelDescription.model_definition ["empirical-rates"] = models.protein.JC.empirical_Q;
    models.protein.JCF.ModelDescription.model_definition ["frequency-estimator"] = "frequencies.empirical.protein";
    models.protein.JCF.ModelDescription.model_definition ["parameters"] = {"global": {}, "local": {}, "empirical": 19};
    models.protein.JCF.ModelDescription.model_definition ["q_ij"] = "models.protein.JC._GenerateRate";
    return models.protein.JCF.ModelDescription.model_definition;
}


/**
 * @name models.protein.JC.frequencies
 * @param {Dictionary} Baseline JC69 model
 * @returns {Dictionary} Updated JC69 model with empirical frequencies
 * @description Define the empirical amino acid frequencies associated with the JC69 model of protein evolution
 */
function models.protein.JC.frequencies (model, namespace, datafilter) {
    model[terms.efv_estimate] =
      {{	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {	0.05}
        {   0.05}
    };
    model[terms.efv_estimate_name] = terms.freqs.predefined;
    (model["parameters"])["empirical"] = 0;
    return model;
}

/* Define a dictionary of equal amino-acid exchangeability rates for the JC69 model of protein evolution. */ 
models.protein.JC.empirical_Q = {
    "A":
        {"C":1.,
        "D":1.,
        "E":1.,
        "F":1.,
        "G":1.,
        "H":1.,
        "I":1.,
        "K":1.,
        "L":1.,
        "M":1.,
        "N":1.,
        "P":1.,
        "Q":1.,
        "R":1.,
        "S":11.,
        "T":11.,
        "V":1.,
        "W":1.,
        "Y":1.},
    "C":
        {"D":1.,
        "E":1.,
        "F":1.,
        "G":1.,
        "H":1.,
        "I":1.,
        "K":1.,
        "L":1.,
        "M":1.,
        "N":1.,
        "P":1.,
        "Q":1.,
        "R":1.,
        "S":1.,
        "T":1.,
        "V":1.,
        "W":1.,
        "Y":1.},
    "D":
        {"E":21.,
        "F":1.,
        "G":1.,
        "H":1.,
        "I":1.,
        "K":1.,
        "L":1.,
        "M":1.,
        "N":11.,
        "P":1.,
        "Q":1.,
        "R":1.,
        "S":1.,
        "T":1.,
        "V":1.,
        "W":1.,
        "Y":1.},
    "E":
        {"F":1.,
        "G":1.,
        "H":1.,
        "I":1.,
        "K":1.,
        "L":1.,
        "M":1.,
        "N":1.,
        "P":1.,
        "Q":11.,
        "R":1.,
        "S":1.,
        "T":1.,
        "V":1.,
        "W":1.,
        "Y":1.},
    "F":
        {"G":1.,
        "H":1.,
        "I":1.,
        "K":1.,
        "L":1.,
        "M":1.,
        "N":1.,
        "P":1.,
        "Q":1.,
        "R":1.,
        "S":1.,
        "T":1.,
        "V":1.,
        "W":1.,
        "Y":11.},
    "G":
        {"H":1.,
        "I":1.,
        "K":1.,
        "L":1.,
        "M":1.,
        "N":1.,
        "P":1.,
        "Q":1.,
        "R":1.,
        "S":1.,
        "T":1.,
        "V":1.,
        "W":1.,
        "Y":1.},
    "H":
        {"I":1.,
        "K":1.,
        "L":1.,
        "M":1.,
        "N":11.,
        "P":1.,
        "Q":11.,
        "R":11.,
        "S":1.,
        "T":1.,
        "V":1.,
        "W":1.,
        "Y":11.},
    "I":
        {"K":1.,
        "L":1.,
        "M":11.,
        "N":1.,
        "P":1.,
        "Q":1.,
        "R":1.,
        "S":1.,
        "T":1.,
        "V":31.,
        "W":1.,
        "Y":1.},
    "K":
        {"L":1.,
        "M":1.,
        "N":1.,
        "P":1.,
        "Q":1.,
        "R":21.,
        "S":1.,
        "T":1.,
        "V":1.,
        "W":1.,
        "Y":1.},
    "L":
        {"M":11.,
        "N":1.,
        "P":1.,
        "Q":1.,
        "R":1.,
        "S":1.,
        "T":1.,
        "V":1.,
        "W":1.,
        "Y":1.},
    "M":
        {"N":1.,
        "P":1.,
        "Q":1.,
        "R":1.,
        "S":1.,
        "T":1.,
        "V":1.,
        "W":1.,
        "Y":1.},
    "N":
        {"P":1.,
        "Q":1.,
        "R":1.,
        "S":11.,
        "T":1.,
        "V":1.,
        "W":1.,
        "Y":1.},
    "P":
        {"Q":1.,
        "R":1.,
        "S":1.,
        "T":1.,
        "V":1.,
        "W":1.,
        "Y":1.},
    "Q":
        {"R":1.,
        "S":1.,
        "T":1.,
        "V":1.,
        "W":1.,
        "Y":1.},
    "R":
        {"S":1.,
        "T":1.,
        "V":1.,
        "W":1.,
        "Y":1.},
    "S":
        {"T":11.,
        "V":1.,
        "W":1.,
        "Y":1.},
    "T":
        {"V":1.,
        "W":1.,
        "Y":1.},
    "V":
        {"W":1.,
        "Y":1.},
    "W":
        {"Y":1.}
};



/**
 * @name models.protein.Empirical._GenerateRate
 * @description Generates the r_ij component of q_ij := r_ij * time * freq_j
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