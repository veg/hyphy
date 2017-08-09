LoadFunctionLibrary("../protein.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");
LoadFunctionLibrary("../../all-terms.bf");

/** @module models.protein.empirical */

/**
 * @name models.protein.empirical.ModelDescription
 * @param {String} type
 * @returns {Dictionary} model description
 * @description Create the baseline schema (dictionary) for empirical protein model definitions
 */
function models.protein.empirical.ModelDescription(type) {
    return {
        terms.alphabet: models.protein.alphabet,
        terms.description: "General class of empirical substitution matrices for amino-acids",
        terms.model.canonical: 1, // is of the r_ij \times \pi_j form
        terms.model.reversible: 1,
        terms.model.efv_estimate_name: terms.frequencies.predefined,
        terms.parameters: {
            terms.global: {},
            terms.local: {},
            terms.model.empirical: 0
        },
        terms.model.type: type,
        terms.model.get_branch_length: "",
        terms.model.set_branch_length: "models.generic.SetBranchLength",
        terms.model.constrain_branch_length: "models.generic.constrain_branch_length",
        terms.model.frequency_estimator: "frequencies.empirical.protein",
        terms.model.q_ij: "",
        terms.model.time: "models.protein.generic.Time",
        terms.model.defineQ: "models.protein.empirical._DefineQ",
        terms.model.post_definition: "models.generic.post.definition"
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
    models.protein.WAG.ModelDescription.model_definition [terms.model.empirical_rates] = models.protein.WAG.empirical_Q;
    models.protein.WAG.ModelDescription.model_definition [terms.model.frequency_estimator] = "models.protein.WAG.frequencies";
    models.protein.WAG.ModelDescription.model_definition [terms.model.q_ij] = "models.protein.WAG._GenerateRate";
    return models.protein.WAG.ModelDescription.model_definition;
}

/**
 * @name models.protein.WAGF.ModelDescription
 * @description Create the baseline schema (dictionary) for the WAG+F model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
/*
function models.protein.WAGF.ModelDescription(type) {
    models.protein.WAGF.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.WAGF.ModelDescription.model_definition ["empirical-rates"] = models.protein.WAG.empirical_Q;
    models.protein.WAGF.ModelDescription.model_definition ["frequency-estimator"] = "frequencies.empirical.protein";
    models.protein.WAGF.ModelDescription.model_definition ["parameters"] = {"global": {}, "local": {}, "empirical": 19};
    models.protein.WAGF.ModelDescription.model_definition ["q_ij"] = "models.protein.WAG._GenerateRate";
    return models.protein.WAGF.ModelDescription.model_definition;
}*/

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

    model[terms.model.efv_estimate_name] = terms.frequencies.predefined;
    (model[terms.parameters])[terms.model.empirical] = 0;
    return model;
}


/* Define a dictionary of amino-acid exchangeability rates for the WAG model of protein evolution. Note that this dictionary has been **normalized** using the WAG default frequencies. */ 
models.protein.WAG.empirical_Q = {
 "A":{
   "C":0.02081175,
   "D":0.04424356,
   "E":0.09644885,
   "F":0.008490242,
   "G":0.1237844,
   "H":0.008127018999999999,
   "I":0.009834135000000001,
   "K":0.05899778,
   "L":0.03600239,
   "M":0.0182884,
   "N":0.02091646,
   "P":0.06909218,
   "Q":0.03502343,
   "R":0.02545459,
   "S":0.245933,
   "T":0.1358225,
   "V":0.1492591,
   "W":0.001708106,
   "Y":0.008912199000000001
  },
 "C":{
   "D":0.001813745,
   "E":0.001301055,
   "F":0.01605407,
   "G":0.02679532,
   "H":0.006383892,
   "I":0.008654049,
   "K":0.004819601,
   "L":0.03476936,
   "M":0.007992529,
   "N":0.0108821,
   "P":0.005254569,
   "Q":0.003809101,
   "R":0.02437562,
   "S":0.1027029,
   "T":0.03284827,
   "V":0.0745652,
   "W":0.01082647,
   "Y":0.02013312
  },
 "D":{
   "E":0.3762142,
   "F":0.001884863,
   "G":0.07562952000000001,
   "H":0.02386347,
   "I":0.002005993,
   "K":0.03123852,
   "L":0.007672926,
   "M":0.002123675,
   "N":0.2227414,
   "P":0.02036354,
   "Q":0.02377493,
   "R":0.00679797,
   "S":0.07819566,
   "T":0.02400406,
   "V":0.01133463,
   "W":0.001959249,
   "Y":0.01205807
  },
 "E":{
   "F":0.003272523,
   "G":0.04960369,
   "H":0.01461601,
   "I":0.006480045,
   "K":0.1682461,
   "L":0.01395734,
   "M":0.006450074,
   "N":0.0388587,
   "P":0.03277285,
   "Q":0.2108299,
   "R":0.02026676,
   "S":0.05143238,
   "T":0.0526847,
   "V":0.0438051,
   "W":0.00236373,
   "Y":0.007267292
  },
 "F":{
   "G":0.00436267,
   "H":0.01741975,
   "I":0.05389076,
   "K":0.005783216,
   "L":0.1913755,
   "M":0.02437025,
   "N":0.00394504,
   "P":0.007754001,
   "Q":0.003851615,
   "R":0.004740036,
   "S":0.03983115,
   "T":0.01100758,
   "V":0.04835584,
   "W":0.02309483,
   "Y":0.2389425
  },
 "G":{
   "H":0.006395123,
   "I":0.001548868,
   "K":0.02431859,
   "L":0.005546612,
   "M":0.003563543,
   "N":0.04617598,
   "P":0.01169843,
   "Q":0.0127224,
   "R":0.02698185,
   "S":0.09789926,
   "T":0.01446092,
   "V":0.01393229,
   "W":0.005087841,
   "Y":0.003835501
  },
 "H":{
   "I":0.007029141,
   "K":0.05796705,
   "L":0.04519012,
   "M":0.008272107000000001,
   "N":0.1623063,
   "P":0.03343772,
   "Q":0.1655236,
   "R":0.09862788,
   "S":0.05400277,
   "T":0.0303076,
   "V":0.008806540999999999,
   "W":0.003964322,
   "Y":0.1433978
  },
 "I":{
   "K":0.02108143,
   "L":0.2869017,
   "M":0.08714326,
   "N":0.02273747,
   "P":0.004799484,
   "Q":0.004391122,
   "R":0.008628942000000001,
   "S":0.02330635,
   "T":0.09337141,
   "V":0.5819514,
   "W":0.003208113,
   "Y":0.01555502
  },
 "K":{
   "L":0.02330296,
   "M":0.0191231,
   "N":0.1235674,
   "P":0.02674718,
   "Q":0.1501354,
   "R":0.246964,
   "S":0.07056185,
   "T":0.08881348999999999,
   "V":0.02272611,
   "W":0.002076079,
   "Y":0.004933538
  },
 "L":{
   "M":0.09935387,
   "N":0.005395922,
   "P":0.01997258,
   "Q":0.03351591,
   "R":0.02296714,
   "S":0.02515217,
   "T":0.02091482,
   "V":0.133956,
   "W":0.01004497,
   "Y":0.01475715
  },
 "M":{
   "N":0.008131996000000001,
   "P":0.008228767999999999,
   "Q":0.05956464,
   "R":0.03152741,
   "S":0.03603533,
   "T":0.0970828,
   "V":0.1531609,
   "W":0.007786238,
   "Y":0.01586107
  },
 "N":{
   "P":0.009369552999999999,
   "Q":0.05950219,
   "R":0.02932074,
   "S":0.28996,
   "T":0.1299922,
   "V":0.01460187,
   "W":0.001085813,
   "Y":0.04020457
  },
 "P":{
   "Q":0.03597839,
   "R":0.03135791,
   "S":0.1177049,
   "T":0.05093139,
   "V":0.02342947,
   "W":0.002104766,
   "Y":0.007998193000000001
  },
 "Q":{
   "R":0.140086,
   "S":0.07506641,
   "T":0.05493632,
   "V":0.0224171,
   "W":0.003257243,
   "Y":0.008430004
  },
 "R":{
   "S":0.08931696,
   "T":0.03550112,
   "V":0.01873906,
   "W":0.01757311,
   "Y":0.01412465
  },
 "S":{
   "T":0.2803409,
   "V":0.01731717,
   "W":0.007907568,
   "Y":0.0291351
  },
 "T":{
   "V":0.1032926,
   "W":0.001673848,
   "Y":0.01077852
  },
 "V":{
   "W":0.005516418,
   "Y":0.01165155
  },
 "W":{
   "Y":0.0920111
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
    models.protein.LG.ModelDescription.model_definition [terms.model.empirical_rates] = models.protein.LG.empirical_Q;
    models.protein.LG.ModelDescription.model_definition [terms.model.frequency_estimator] = "models.protein.LG.frequencies";
    models.protein.LG.ModelDescription.model_definition [terms.model.q_ij] = "models.protein.LG._GenerateRate";
    return models.protein.LG.ModelDescription.model_definition;
}

/**
 * @name models.protein.LGF.ModelDescription
 * @description Create the baseline schema (dictionary) for the LG+F model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
 /*
function models.protein.LGF.ModelDescription(type) {
    models.protein.LGF.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.LGF.ModelDescription.model_definition ["empirical-rates"] = models.protein.LG.empirical_Q;
    models.protein.LGF.ModelDescription.model_definition ["frequency-estimator"] = "frequencies.empirical.protein";
    models.protein.LGF.ModelDescription.model_definition ["parameters"] = {"global": {}, "local": {}, "empirical": 19};
    models.protein.LGF.ModelDescription.model_definition ["q_ij"] = "models.protein.LG._GenerateRate";
    return models.protein.LGF.ModelDescription.model_definition;
}
*/


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

    model[terms.model.efv_estimate_name] = terms.frequencies.predefined;
    (model[terms.parameters])[terms.model.empirical] = 0;
    return model;
}

/* Define a dictionary of amino-acid exchangeability rates for the LG model of protein evolution. Note that this dictionary has been **normalized** using the LG default frequencies. */ 
models.protein.LG.empirical_Q = {
 "A":{
   "C":0.03522956,
   "D":0.0229346,
   "E":0.08133689,
   "F":0.01174132,
   "G":0.1296008,
   "H":0.008776704999999999,
   "I":0.01018879,
   "K":0.03791848,
   "L":0.04285406,
   "M":0.02822381,
   "N":0.01271276,
   "P":0.05674114,
   "Q":0.04325807,
   "R":0.02601647,
   "S":0.3164948,
   "T":0.1247291,
   "V":0.1927457,
   "W":0.002385593,
   "Y":0.008181845
  },
 "C":{
   "D":0.003630821,
   "E":0.0002740351,
   "F":0.05115122,
   "G":0.03570948,
   "H":0.01566596,
   "I":0.0218034,
   "K":0.0009375764,
   "L":0.06438966,
   "M":0.02243974,
   "N":0.02428347,
   "P":0.003632027,
   "Q":0.003782507,
   "R":0.0327155,
   "S":0.1864267,
   "T":0.06666287,
   "V":0.1482198,
   "W":0.00884617,
   "Y":0.04355245
  },
 "D":{
   "E":0.41069,
   "F":0.0008060157,
   "G":0.05300146,
   "H":0.02267472,
   "I":0.0007269456,
   "K":0.01999816,
   "L":0.00163422,
   "M":0.0006414941,
   "N":0.2331202,
   "P":0.01900553,
   "Q":0.02334345,
   "R":0.007586211,
   "S":0.08303903999999999,
   "T":0.02482688,
   "V":0.002872194,
   "W":0.0003945694,
   "Y":0.005048546
  },
 "E":{
   "F":0.0008705765,
   "G":0.02188286,
   "H":0.01036699,
   "I":0.003010126,
   "K":0.1277224,
   "L":0.007552471,
   "M":0.004362376,
   "N":0.02487791,
   "P":0.02020781,
   "Q":0.1841385,
   "R":0.02227563,
   "S":0.04097288,
   "T":0.03524391,
   "V":0.01853676,
   "W":0.001027702,
   "Y":0.004485425
  },
 "F":{
   "G":0.00561965,
   "H":0.01668329,
   "I":0.0756681,
   "K":0.001690408,
   "L":0.2810447,
   "M":0.04516806,
   "N":0.004111401,
   "P":0.004551429,
   "Q":0.001599162,
   "R":0.003226682,
   "S":0.02422454,
   "T":0.009619268,
   "V":0.04952661,
   "W":0.03243575,
   "Y":0.2916085
  },
 "G":{
   "H":0.007618063,
   "I":0.0005919608,
   "K":0.02096479,
   "L":0.00479784,
   "M":0.003503711,
   "N":0.06602329999999999,
   "P":0.009489902,
   "Q":0.01195119,
   "R":0.02388046,
   "S":0.116496,
   "T":0.00756921,
   "V":0.005802412,
   "W":0.003544273,
   "Y":0.002043191
  },
 "H":{
   "I":0.007404237,
   "K":0.04927923,
   "L":0.03970833,
   "M":0.01111019,
   "N":0.207085,
   "P":0.02451727,
   "Q":0.2146863,
   "R":0.1485124,
   "S":0.06628340000000001,
   "T":0.03406144,
   "V":0.009003304,
   "W":0.007881539999999999,
   "Y":0.1983005
  },
 "I":{
   "K":0.01124222,
   "L":0.4493203,
   "M":0.1073075,
   "N":0.008794702,
   "P":0.003771706,
   "Q":0.003249348,
   "R":0.007772081,
   "S":0.004291965,
   "T":0.06026516,
   "V":0.04910478,
   "W":0.001473992,
   "Y":0.008688692
  },
 "K":{
   "L":0.01490483,
   "M":0.01648691,
   "N":0.09851189,
   "P":0.01880635,
   "Q":0.1442522,
   "R":0.3871668,
   "S":0.05012591,
   "T":0.06627711,
   "V":0.01401048,
   "W":0.0006587949,
   "Y":0.004929905
  },
 "L":{
   "M":0.1584993,
   "N":0.003142484,
   "P":0.01200011,
   "Q":0.02597806,
   "R":0.01847365,
   "S":0.0122045,
   "T":0.01766063,
   "V":0.1288122,
   "W":0.008179586000000001,
   "Y":0.01119695
  },
 "M":{
   "N":0.01703821,
   "P":0.004810887,
   "Q":0.07459796,
   "R":0.02962982,
   "S":0.0232297,
   "T":0.1177837,
   "V":0.1436375,
   "W":0.009190009000000001,
   "Y":0.01798497
  },
 "N":{
   "P":0.007795161,
   "Q":0.07563193999999999,
   "R":0.04601631,
   "S":0.268368,
   "T":0.116636,
   "V":0.006330976,
   "W":0.0005989957,
   "Y":0.02286955
  },
 "P":{
   "Q":0.02784403,
   "R":0.02035162,
   "S":0.08959077,
   "T":0.03331558,
   "V":0.02243022,
   "W":0.001255797,
   "Y":0.00334857
  },
 "Q":{
   "R":0.1718491,
   "S":0.08193787,
   "T":0.06297003,
   "V":0.01591156,
   "W":0.003117996,
   "Y":0.009615879000000001
  },
 "R":{
   "S":0.05745502,
   "T":0.03375392,
   "V":0.01292756,
   "W":0.007836037000000001,
   "Y":0.01174968
  },
 "S":{
   "T":0.3773224,
   "V":0.00744159,
   "W":0.003285156,
   "Y":0.01496724
  },
 "T":{
   "V":0.1655336,
   "W":0.00185899,
   "Y":0.009186346
  },
 "V":{
   "W":0.002501667,
   "Y":0.009316084000000001
  },
 "W":{
   "Y":0.1177739
  }
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
    models.protein.JTT.ModelDescription.model_definition [terms.model.empirical_rates] = models.protein.JTT.empirical_Q;
    models.protein.JTT.ModelDescription.model_definition [terms.model.frequency_estimator] = "models.protein.JTT.frequencies";
    models.protein.JTT.ModelDescription.model_definition [terms.model.q_ij] = "models.protein.JTT._GenerateRate";
    return models.protein.JTT.ModelDescription.model_definition;
}

/**
 * @name models.protein.JTTF.ModelDescription
 * @description Create the baseline schema (dictionary) for the JTT+F model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
 /*
function models.protein.JTTF.ModelDescription(type) {
    models.protein.JTTF.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.JTTF.ModelDescription.model_definition ["empirical-rates"] = models.protein.JTT.empirical_Q;
    models.protein.JTTF.ModelDescription.model_definition ["frequency-estimator"] = "frequencies.empirical.protein";
    models.protein.JTTF.ModelDescription.model_definition ["parameters"] = {"global": {}, "local": {}, "empirical": 19};
    models.protein.JTTF.ModelDescription.model_definition ["q_ij"] = "models.protein.JTT._GenerateRate";
    return models.protein.JTTF.ModelDescription.model_definition;
}*/


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

    model[terms.model.efv_estimate_name] = terms.frequencies.predefined;
    (model[terms.parameters])[terms.model.empirical] = 0;
    return model;
}

/* Define a dictionary of amino-acid exchangeability rates for the JTT model of protein evolution. Note that this dictionary has been **normalized** using the JTT default frequencies. */ 
models.protein.JTT.empirical_Q = {
 "A":{
   "C":0.01164983,
   "D":0.04242226,
   "E":0.06594220000000001,
   "F":0.005605013,
   "G":0.1300142,
   "H":0.005055569,
   "I":0.01901336,
   "K":0.02198075,
   "L":0.02824504,
   "M":0.01099041,
   "N":0.02373925,
   "P":0.09902242999999999,
   "Q":0.02285967,
   "R":0.02714587,
   "S":0.2651969,
   "T":0.2681624,
   "V":0.1940882,
   "W":0.00120894,
   "Y":0.004506008
  },
 "C":{
   "D":0.005415286,
   "E":0.003332529,
   "F":0.02749291,
   "G":0.04082289,
   "H":0.01666262,
   "I":0.007914733,
   "K":0.002915936,
   "L":0.01499622,
   "M":0.009581053000000001,
   "N":0.01333012,
   "P":0.006248431,
   "Q":0.003749032,
   "R":0.05207011,
   "S":0.1470474,
   "T":0.02749309,
   "V":0.04123968,
   "W":0.01582953,
   "Y":0.06831603999999999
  },
 "D":{
   "E":0.4801284,
   "F":0.001318116,
   "G":0.0950686,
   "H":0.0237263,
   "I":0.00609632,
   "K":0.01680615,
   "L":0.005602049,
   "M":0.004448682,
   "N":0.2361102,
   "P":0.006425849,
   "Q":0.0214193,
   "R":0.007908676,
   "S":0.04020279,
   "T":0.02487944,
   "V":0.02092512,
   "W":0.0008238323,
   "Y":0.014664
  },
 "E":{
   "F":0.001776388,
   "G":0.08335330000000001,
   "H":0.005602518,
   "I":0.005875793,
   "K":0.1030317,
   "L":0.008881953,
   "M":0.004099415,
   "N":0.02459647,
   "P":0.009701838000000001,
   "Q":0.1403343,
   "R":0.01626078,
   "S":0.02131682,
   "T":0.01940362,
   "V":0.03088188,
   "W":0.001639765,
   "Y":0.002049689
  },
 "F":{
   "G":0.003751538,
   "H":0.01042113,
   "I":0.04085083,
   "K":0.00145895,
   "L":0.2278042,
   "M":0.01021273,
   "N":0.003126321,
   "P":0.00750314,
   "Q":0.001875789,
   "R":0.003334736,
   "S":0.0644024,
   "T":0.008128382,
   "V":0.03939149,
   "W":0.007711647,
   "Y":0.1771572
  },
 "G":{
   "H":0.004635577,
   "I":0.002826581,
   "K":0.01605493,
   "L":0.006331483,
   "M":0.003052693,
   "N":0.03290136,
   "P":0.01051474,
   "Q":0.009497159999999999,
   "R":0.06941973,
   "S":0.1278738,
   "T":0.01854212,
   "V":0.03120506,
   "W":0.007801362,
   "Y":0.001695907
  },
 "H":{
   "I":0.00955641,
   "K":0.03124215,
   "L":0.04925195,
   "M":0.007718657,
   "N":0.1712807,
   "P":0.05770555,
   "Q":0.2333939,
   "R":0.1639271,
   "S":0.0507224,
   "T":0.02793385,
   "V":0.008086143,
   "W":0.001837774,
   "Y":0.1889208
  },
 "I":{
   "K":0.01205203,
   "L":0.2127567,
   "M":0.1131285,
   "N":0.0208902,
   "P":0.004981443,
   "Q":0.003213843,
   "R":0.01221257,
   "S":0.02763923,
   "T":0.1494435,
   "V":0.6328057,
   "W":0.001928334,
   "Y":0.009802181
  },
 "K":{
   "L":0.01334602,
   "M":0.01462393,
   "N":0.1076208,
   "P":0.01093234,
   "Q":0.1218169,
   "R":0.333364,
   "S":0.03237125,
   "T":0.05650736,
   "V":0.008234754,
   "W":0.001277824,
   "Y":0.002839562
  },
 "L":{
   "M":0.09030557,
   "N":0.005841096,
   "P":0.05358937,
   "Q":0.0291124,
   "R":0.01900652,
   "S":0.04042405,
   "T":0.01594697,
   "V":0.1169137,
   "W":0.007602722,
   "Y":0.007788057
  },
 "M":{
   "N":0.0140708,
   "P":0.008298109999999999,
   "Q":0.0187608,
   "R":0.02200785,
   "S":0.01948259,
   "T":0.1237496,
   "V":0.2016795,
   "W":0.002886323,
   "Y":0.006133368
  },
 "N":{
   "P":0.006154998,
   "Q":0.03156908,
   "R":0.02303155,
   "S":0.3450795,
   "T":0.1375939,
   "V":0.01092017,
   "W":0.000397107,
   "Y":0.02263447
  },
 "P":{
   "Q":0.06603124,
   "R":0.03627542,
   "S":0.1902389,
   "T":0.06887338,
   "V":0.01404214,
   "W":0.001003017,
   "Y":0.003677695
  },
 "Q":{
   "R":0.1542939,
   "S":0.03744234,
   "T":0.03065318,
   "V":0.01193211,
   "W":0.002468744,
   "Y":0.008229024999999999
  },
 "R":{
   "S":0.06833079,
   "T":0.03805319,
   "V":0.01141599,
   "W":0.01803412,
   "Y":0.007610617
  },
 "S":{
   "T":0.2795782,
   "V":0.02711589,
   "W":0.004457448,
   "Y":0.02030591
  },
 "T":{
   "V":0.0759305,
   "W":0.00115485,
   "Y":0.006495937
  },
 "V":{
   "W":0.003436295,
   "Y":0.005345272
  },
 "W":{
   "Y":0.02415905
  }
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
    models.protein.JC.ModelDescription.model_definition [terms.model.empirical_rates] = models.protein.JC.empirical_Q;
    models.protein.JC.ModelDescription.model_definition [terms.model.frequency_estimator] = "models.protein.JC.frequencies";
    models.protein.JC.ModelDescription.model_definition [terms.model.q_ij] = "models.protein.JC._GenerateRate";
    return models.protein.JC.ModelDescription.model_definition;
}

/**
 * @name models.protein.JCF.ModelDescription
 * @description Create the baseline schema (dictionary) for the JC69+F model of protein evolution
 * @returns {Dictionary} model description
 * @param {String} type
 */
 /*
function models.protein.JCF.ModelDescription(type) {
    models.protein.JCF.ModelDescription.model_definition = models.protein.empirical.ModelDescription(type);
    models.protein.JCF.ModelDescription.model_definition ["empirical-rates"] = models.protein.JC.empirical_Q;
    models.protein.JCF.ModelDescription.model_definition ["frequency-estimator"] = "frequencies.empirical.protein";
    models.protein.JCF.ModelDescription.model_definition ["parameters"] = {"global": {}, "local": {}, "empirical": 19};
    models.protein.JCF.ModelDescription.model_definition ["q_ij"] = "models.protein.JC._GenerateRate";
    return models.protein.JCF.ModelDescription.model_definition;
}*/


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

    model[terms.model.efv_estimate_name] = terms.frequencies.predefined;
    (model[terms.parameters])[terms.model.empirical] = 0;
    return model;
}

/* Define a dictionary of equal amino-acid exchangeability rates for the JC69 model of protein evolution.  Note that this dictionary has been **normalized** using the JC69 default frequencies. */ 
models.protein.JC.empirical_Q = {
    "A":
        {"C":0.05,
        "D":0.05,
        "E":0.05,
        "F":0.05,
        "G":0.05,
        "H":0.05,
        "I":0.05,
        "K":0.05,
        "L":0.05,
        "M":0.05,
        "N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "C":
        {"D":0.05,
        "E":0.05,
        "F":0.05,
        "G":0.05,
        "H":0.05,
        "I":0.05,
        "K":0.05,
        "L":0.05,
        "M":0.05,
        "N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "D":
        {"E":20.05,
        "F":0.05,
        "G":0.05,
        "H":0.05,
        "I":0.05,
        "K":0.05,
        "L":0.05,
        "M":0.05,
        "N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "E":
        {"F":0.05,
        "G":0.05,
        "H":0.05,
        "I":0.05,
        "K":0.05,
        "L":0.05,
        "M":0.05,
        "N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "F":
        {"G":0.05,
        "H":0.05,
        "I":0.05,
        "K":0.05,
        "L":0.05,
        "M":0.05,
        "N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "G":
        {"H":0.05,
        "I":0.05,
        "K":0.05,
        "L":0.05,
        "M":0.05,
        "N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "H":
        {"I":0.05,
        "K":0.05,
        "L":0.05,
        "M":0.05,
        "N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "I":
        {"K":0.05,
        "L":0.05,
        "M":0.05,
        "N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":30.05,
        "W":0.05,
        "Y":0.05},
    "K":
        {"L":0.05,
        "M":0.05,
        "N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":20.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "L":
        {"M":0.05,
        "N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "M":
        {"N":0.05,
        "P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "N":
        {"P":0.05,
        "Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "P":
        {"Q":0.05,
        "R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "Q":
        {"R":0.05,
        "S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "R":
        {"S":0.05,
        "T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "S":
        {"T":0.05,
        "V":0.05,
        "W":0.05,
        "Y":0.05},
    "T":
        {"V":0.05,
        "W":0.05,
        "Y":0.05},
    "V":
        {"W":0.05,
        "Y":0.05},
    "W":
        {"Y":0.05}
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
        models.protein.empirical._GenerateRate.p  [terms.model.rate_entry] = "" + (rateDict[fromChar])[toChar];
    } else {
        models.protein.empirical._GenerateRate.p  [terms.model.rate_entry] = "" + (rateDict[toChar])[fromChar];
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


/************************************ R_ij matrices ***************************************/


models.protein.empirical.WAG.empirical_R = {{0., 1.027039999999999953, 7.389980000000000437e-01, 1.58285000000000009, 2.104939999999999867e-01, 1.416719999999999979, 3.169540000000000135e-01, 1.933350000000000068e-01, 9.062649999999999872e-01, 3.979150000000000187e-01, 8.934959999999999569e-01, 5.098479999999999679e-01, 1.438549999999999995, 9.085980000000000167e-01, 5.515710000000000335e-01, 3.370789999999999953, 2.121109999999999829, 2.006009999999999849, 1.131329999999999975e-01, 2.407350000000000045e-01},
			{1.027039999999999953, 0., 3.029489999999999961e-02, 2.135199999999999945e-02, 3.980199999999999849e-01, 3.066740000000000022e-01, 2.489719999999999989e-01, 1.701350000000000084e-01, 7.403389999999999971e-02, 3.842869999999999897e-01, 3.904819999999999958e-01, 2.652559999999999918e-01, 1.094040000000000012e-01, 9.881790000000000018e-02, 5.281909999999999661e-01, 1.407659999999999911, 5.129839999999999955e-01, 1.002140000000000031, 7.170699999999999852e-01, 5.438330000000000108e-01},
			{7.389980000000000437e-01, 3.029489999999999961e-02, 0., 6.174159999999999648, 4.673039999999999833e-02, 8.655840000000000201e-01, 9.306759999999999478e-01, 3.943699999999999983e-02, 4.79854999999999976e-01, 8.480469999999999675e-02, 1.03753999999999999e-01, 5.429420000000000357, 4.239840000000000275e-01, 6.167829999999999702e-01, 1.473039999999999905e-01, 1.071760000000000046, 3.74865999999999977e-01, 1.523349999999999982e-01, 1.297669999999999935e-01, 3.257109999999999728e-01},
			{1.58285000000000009, 2.135199999999999945e-02, 6.174159999999999648, 0., 8.113389999999999491e-02, 5.677170000000000272e-01, 5.700250000000000039e-01, 1.273950000000000082e-01, 2.584429999999999783, 1.542630000000000112e-01, 3.151240000000000152e-01, 9.471979999999999844e-01, 6.823550000000000448e-01, 5.469470000000000276, 4.391570000000000196e-01, 7.049389999999999823e-01, 8.227649999999999686e-01, 5.887310000000000043e-01, 1.565570000000000017e-01, 1.963030000000000053e-01},
			{2.104939999999999867e-01, 3.980199999999999849e-01, 4.673039999999999833e-02, 8.113389999999999491e-02, 0., 4.993100000000000316e-02, 6.793709999999999471e-01, 1.059469999999999912, 8.88359999999999983e-02, 2.115169999999999995, 1.190630000000000077, 9.616210000000000035e-02, 1.614440000000000042e-01, 9.992080000000000406e-02, 1.027109999999999967e-01, 5.459310000000000551e-01, 1.719030000000000002e-01, 6.498920000000000252e-01, 1.529640000000000111, 6.454279999999999795},
			{1.416719999999999979, 3.066740000000000022e-01, 8.655840000000000201e-01, 5.677170000000000272e-01, 4.993100000000000316e-02, 0., 2.494099999999999928e-01, 3.045010000000000078e-02, 3.735580000000000012e-01, 6.130370000000000263e-02, 1.741000000000000048e-01, 1.125559999999999894, 2.435700000000000087e-01, 3.30052000000000012e-01, 5.846649999999999903e-01, 1.341820000000000013, 2.25833000000000006e-01, 1.872469999999999968e-01, 3.36982999999999977e-01, 1.036040000000000016e-01},
			{3.169540000000000135e-01, 2.489719999999999989e-01, 9.306759999999999478e-01, 5.700250000000000039e-01, 6.793709999999999471e-01, 2.494099999999999928e-01, 0., 1.381900000000000073e-01, 8.904320000000000013e-01, 4.994620000000000171e-01, 4.041409999999999725e-01, 3.956290000000000084, 6.961979999999999835e-01, 4.294109999999999872, 2.137150000000000105, 7.401689999999999658e-01, 4.733069999999999777e-01, 1.183580000000000049e-01, 2.625689999999999968e-01, 3.873439999999999994},
			{1.933350000000000068e-01, 1.701350000000000084e-01, 3.943699999999999983e-02, 1.273950000000000082e-01, 1.059469999999999912, 3.045010000000000078e-02, 1.381900000000000073e-01, 0., 3.238320000000000087e-01, 3.170970000000000066, 4.257460000000000022, 5.542359999999999509e-01, 9.992879999999999818e-02, 1.139170000000000044e-01, 1.869790000000000063e-01, 3.194400000000000017e-01, 1.458159999999999901, 7.821299999999999919, 2.124830000000000052e-01, 4.201699999999999879e-01},
			{9.062649999999999872e-01, 7.403389999999999971e-02, 4.79854999999999976e-01, 2.584429999999999783, 8.88359999999999983e-02, 3.735580000000000012e-01, 8.904320000000000013e-01, 3.238320000000000087e-01, 0., 2.575549999999999784e-01, 9.342759999999999954e-01, 3.012010000000000076, 5.568959999999999466e-01, 3.894899999999999807, 5.351420000000000066, 9.671300000000000452e-01, 1.386980000000000102, 3.054339999999999833e-01, 1.375049999999999883e-01, 1.332639999999999936e-01},
			{3.979150000000000187e-01, 3.842869999999999897e-01, 8.480469999999999675e-02, 1.542630000000000112e-01, 2.115169999999999995, 6.130370000000000263e-02, 4.994620000000000171e-01, 3.170970000000000066, 2.575549999999999784e-01, 0., 4.854020000000000223, 1.31528000000000006e-01, 4.158439999999999914e-01, 8.694889999999999564e-01, 4.976709999999999745e-01, 3.447390000000000176e-01, 3.266220000000000234e-01, 1.800340000000000051, 6.653090000000000392e-01, 3.986179999999999723e-01},
			{8.934959999999999569e-01, 3.904819999999999958e-01, 1.03753999999999999e-01, 3.151240000000000152e-01, 1.190630000000000077, 1.741000000000000048e-01, 4.041409999999999725e-01, 4.257460000000000022, 9.342759999999999954e-01, 4.854020000000000223, 0., 1.982210000000000083e-01, 1.713290000000000091e-01, 1.545260000000000078, 6.83162000000000047e-01, 4.939049999999999829e-01, 1.516119999999999912, 2.058450000000000113, 5.157059999999999977e-01, 4.284370000000000123e-01},
			{5.098479999999999679e-01, 2.652559999999999918e-01, 5.429420000000000357, 9.471979999999999844e-01, 9.616210000000000035e-02, 1.125559999999999894, 3.956290000000000084, 5.542359999999999509e-01, 3.012010000000000076, 1.31528000000000006e-01, 1.982210000000000083e-01, 0., 1.950810000000000044e-01, 1.543639999999999901, 6.353459999999999663e-01, 3.974229999999999929, 2.030060000000000198, 1.962460000000000038e-01, 7.191670000000000007e-02, 1.086000000000000076},
			{1.438549999999999995, 1.094040000000000012e-01, 4.239840000000000275e-01, 6.823550000000000448e-01, 1.614440000000000042e-01, 2.435700000000000087e-01, 6.961979999999999835e-01, 9.992879999999999818e-02, 5.568959999999999466e-01, 4.158439999999999914e-01, 1.713290000000000091e-01, 1.950810000000000044e-01, 0., 9.333719999999999795e-01, 6.794890000000000096e-01, 1.613280000000000047, 7.953839999999999799e-01, 3.148869999999999725e-01, 1.394050000000000011e-01, 2.160459999999999881e-01},
			{9.085980000000000167e-01, 9.881790000000000018e-02, 6.167829999999999702e-01, 5.469470000000000276, 9.992080000000000406e-02, 3.30052000000000012e-01, 4.294109999999999872, 1.139170000000000044e-01, 3.894899999999999807, 8.694889999999999564e-01, 1.545260000000000078, 1.543639999999999901, 9.333719999999999795e-01, 0., 3.035499999999999865, 1.028869999999999951, 8.579280000000000239e-01, 3.01281000000000021e-01, 2.157370000000000121e-01, 2.277099999999999957e-01},
			{5.515710000000000335e-01, 5.281909999999999661e-01, 1.473039999999999905e-01, 4.391570000000000196e-01, 1.027109999999999967e-01, 5.846649999999999903e-01, 2.137150000000000105, 1.869790000000000063e-01, 5.351420000000000066, 4.976709999999999745e-01, 6.83162000000000047e-01, 6.353459999999999663e-01, 6.794890000000000096e-01, 3.035499999999999865, 0., 1.224189999999999889, 5.544130000000000447e-01, 2.518489999999999895e-01, 1.163920000000000066, 3.81533000000000011e-01},
			{3.370789999999999953, 1.407659999999999911, 1.071760000000000046, 7.049389999999999823e-01, 5.459310000000000551e-01, 1.341820000000000013, 7.401689999999999658e-01, 3.194400000000000017e-01, 9.671300000000000452e-01, 3.447390000000000176e-01, 4.939049999999999829e-01, 3.974229999999999929, 1.613280000000000047, 1.028869999999999951, 1.224189999999999889, 0., 4.378020000000000245, 2.327390000000000014e-01, 5.237420000000000408e-01, 7.869930000000000536e-01},
			{2.121109999999999829, 5.129839999999999955e-01, 3.74865999999999977e-01, 8.227649999999999686e-01, 1.719030000000000002e-01, 2.25833000000000006e-01, 4.733069999999999777e-01, 1.458159999999999901, 1.386980000000000102, 3.266220000000000234e-01, 1.516119999999999912, 2.030060000000000198, 7.953839999999999799e-01, 8.579280000000000239e-01, 5.544130000000000447e-01, 4.378020000000000245, 0., 1.388230000000000075, 1.108640000000000042e-01, 2.911480000000000179e-01},
			{2.006009999999999849, 1.002140000000000031, 1.523349999999999982e-01, 5.887310000000000043e-01, 6.498920000000000252e-01, 1.872469999999999968e-01, 1.183580000000000049e-01, 7.821299999999999919, 3.054339999999999833e-01, 1.800340000000000051, 2.058450000000000113, 1.962460000000000038e-01, 3.148869999999999725e-01, 3.01281000000000021e-01, 2.518489999999999895e-01, 2.327390000000000014e-01, 1.388230000000000075, 0., 3.653689999999999993e-01, 3.147300000000000098e-01},
			{1.131329999999999975e-01, 7.170699999999999852e-01, 1.297669999999999935e-01, 1.565570000000000017e-01, 1.529640000000000111, 3.36982999999999977e-01, 2.625689999999999968e-01, 2.124830000000000052e-01, 1.375049999999999883e-01, 6.653090000000000392e-01, 5.157059999999999977e-01, 7.191670000000000007e-02, 1.394050000000000011e-01, 2.157370000000000121e-01, 1.163920000000000066, 5.237420000000000408e-01, 1.108640000000000042e-01, 3.653689999999999993e-01, 0., 2.48539000000000021},
			{2.407350000000000045e-01, 5.438330000000000108e-01, 3.257109999999999728e-01, 1.963030000000000053e-01, 6.454279999999999795, 1.036040000000000016e-01, 3.873439999999999994, 4.201699999999999879e-01, 1.332639999999999936e-01, 3.986179999999999723e-01, 4.284370000000000123e-01, 1.086000000000000076, 2.160459999999999881e-01, 2.277099999999999957e-01, 3.81533000000000011e-01, 7.869930000000000536e-01, 2.911480000000000179e-01, 3.147300000000000098e-01, 2.48539000000000021, 0.}};

models.protein.empirical.LG.empirical_R = {{0., 2.489084000000000074, 3.951439999999999952e-01, 1.038545000000000051, 2.537010000000000098e-01, 2.066040000000000099, 3.588580000000000103e-01, 1.49829999999999991e-01, 5.365180000000000504e-01, 3.95336999999999994e-01, 1.124034999999999895, 2.768180000000000085e-01, 1.177651000000000003, 9.698940000000000339e-01, 4.250929999999999986e-01, 4.727181999999999995, 2.139501000000000097, 2.547870000000000079, 1.807169999999999888e-01, 2.18958999999999987e-01},
			{2.489084000000000074, 0., 6.255600000000000049e-02, 3.49899999999999994e-03, 1.105250999999999983, 5.69265000000000021e-01, 6.405429999999999735e-01, 3.206269999999999953e-01, 1.326600000000000001e-02, 5.940069999999999517e-01, 8.936800000000000299e-01, 5.287680000000000158e-01, 7.538200000000000456e-02, 8.48079999999999945e-02, 5.345509999999999984e-01, 2.784478000000000009, 1.143480000000000052, 1.959290999999999894, 6.701279999999999459e-01, 1.165532000000000012},
			{3.951439999999999952e-01, 6.255600000000000049e-02, 0., 5.243870000000000253, 1.741600000000000092e-02, 8.449259999999999549e-01, 9.271139999999999937e-01, 1.068999999999999985e-02, 2.829590000000000161e-01, 1.507600000000000086e-02, 2.554800000000000126e-02, 5.076149000000000022, 3.944559999999999733e-01, 5.233860000000000179e-01, 1.239539999999999947e-01, 1.240275000000000016, 4.258600000000000163e-01, 3.796700000000000075e-02, 2.98899999999999999e-02, 1.351070000000000049e-01},
			{1.038545000000000051, 3.49899999999999994e-03, 5.243870000000000253, 0., 1.881100000000000133e-02, 3.488470000000000182e-01, 4.238810000000000078e-01, 4.426499999999999879e-02, 1.807177000000000033, 6.967299999999999882e-02, 1.737350000000000005e-01, 5.417119999999999713e-01, 4.194089999999999763e-01, 4.128591000000000122, 3.639700000000000157e-01, 6.119729999999999892e-01, 6.04544999999999999e-01, 2.450340000000000018e-01, 7.785200000000000453e-02, 1.200370000000000048e-01},
			{2.537010000000000098e-01, 1.105250999999999983, 1.741600000000000092e-02, 1.881100000000000133e-02, 0., 8.958599999999999897e-02, 6.821390000000000509e-01, 1.112727000000000022, 2.391799999999999829e-02, 2.592691999999999997, 1.798853000000000035, 8.952499999999999347e-02, 9.44640000000000063e-02, 3.585499999999999798e-02, 5.272199999999999803e-02, 3.618190000000000017e-01, 1.650010000000000088e-01, 6.546830000000000149e-01, 2.457120999999999889, 7.803901999999999894},
			{2.066040000000000099, 5.69265000000000021e-01, 8.449259999999999549e-01, 3.488470000000000182e-01, 8.958599999999999897e-02, 0., 3.114839999999999831e-01, 8.704999999999999197e-03, 2.966360000000000108e-01, 4.426100000000000173e-02, 1.395379999999999954e-01, 1.437645000000000062, 1.969609999999999972e-01, 2.679590000000000027e-01, 3.901919999999999833e-01, 1.739989999999999926, 1.29836000000000007e-01, 7.670100000000000529e-02, 2.684909999999999797e-01, 5.467899999999999844e-02},
			{3.588580000000000103e-01, 6.405429999999999735e-01, 9.271139999999999937e-01, 4.238810000000000078e-01, 6.821390000000000509e-01, 3.114839999999999831e-01, 0., 1.088820000000000066e-01, 6.972639999999999949e-01, 3.663170000000000037e-01, 4.424719999999999764e-01, 4.509237999999999857, 5.088510000000000533e-01, 4.813505000000000145, 2.426600999999999786, 9.900120000000000031e-01, 5.842619999999999481e-01, 1.190129999999999938e-01, 5.970539999999999736e-01, 5.306834000000000273},
			{1.49829999999999991e-01, 3.206269999999999953e-01, 1.068999999999999985e-02, 4.426499999999999879e-02, 1.112727000000000022, 8.704999999999999197e-03, 1.088820000000000066e-01, 0., 1.590689999999999882e-01, 4.145067000000000057, 4.273607000000000156, 1.915030000000000066e-01, 7.828100000000000336e-02, 7.285400000000000209e-02, 1.269909999999999928e-01, 6.410499999999999532e-02, 1.033738999999999963, 6.491069999999999895e-01, 1.116599999999999954e-01, 2.325230000000000075e-01},
			{5.365180000000000504e-01, 1.326600000000000001e-02, 2.829590000000000161e-01, 1.807177000000000033, 2.391799999999999829e-02, 2.966360000000000108e-01, 6.972639999999999949e-01, 1.590689999999999882e-01, 0., 1.375000000000000111e-01, 6.566039999999999655e-01, 2.145077999999999818, 3.903220000000000023e-01, 3.23429399999999978, 6.326067000000000107, 7.486829999999999874e-01, 1.136862999999999957, 1.852020000000000055e-01, 4.990599999999999897e-02, 1.319319999999999937e-01},
			{3.95336999999999994e-01, 5.940069999999999517e-01, 1.507600000000000086e-02, 6.967299999999999882e-02, 2.592691999999999997, 4.426100000000000173e-02, 3.663170000000000037e-01, 4.145067000000000057, 1.375000000000000111e-01, 0., 6.312357999999999691, 6.842700000000000171e-02, 2.490600000000000036e-01, 5.824570000000000025e-01, 3.018480000000000052e-01, 1.822870000000000046e-01, 3.029359999999999831e-01, 1.702744999999999953, 6.196319999999999606e-01, 2.996480000000000254e-01},
			{1.124034999999999895, 8.936800000000000299e-01, 2.554800000000000126e-02, 1.737350000000000005e-01, 1.798853000000000035, 1.395379999999999954e-01, 4.424719999999999764e-01, 4.273607000000000156, 6.566039999999999655e-01, 6.312357999999999691, 0., 3.710040000000000004e-01, 9.984899999999999332e-02, 1.672568999999999972, 4.8413299999999998e-01, 3.469599999999999906e-01, 2.020366000000000106, 1.898717999999999906, 6.961749999999999883e-01, 4.813060000000000116e-01},
			{2.768180000000000085e-01, 5.287680000000000158e-01, 5.076149000000000022, 5.417119999999999713e-01, 8.952499999999999347e-02, 1.437645000000000062, 4.509237999999999857, 1.915030000000000066e-01, 2.145077999999999818, 6.842700000000000171e-02, 3.710040000000000004e-01, 0., 1.617869999999999864e-01, 1.695751999999999926, 7.518780000000000463e-01, 4.00835800000000031, 2.000678999999999874, 8.36879999999999985e-02, 4.537599999999999967e-02, 6.120250000000000412e-01},
			{1.177651000000000003, 7.538200000000000456e-02, 3.944559999999999733e-01, 4.194089999999999763e-01, 9.44640000000000063e-02, 1.969609999999999972e-01, 5.088510000000000533e-01, 7.828100000000000336e-02, 3.903220000000000023e-01, 2.490600000000000036e-01, 9.984899999999999332e-02, 1.617869999999999864e-01, 0., 6.242940000000000156e-01, 3.32533000000000023e-01, 1.338132000000000099, 5.71467999999999976e-01, 2.965010000000000145e-01, 9.513099999999999334e-02, 8.961299999999999821e-02},
			{9.698940000000000339e-01, 8.48079999999999945e-02, 5.233860000000000179e-01, 4.128591000000000122, 3.585499999999999798e-02, 2.679590000000000027e-01, 4.813505000000000145, 7.285400000000000209e-02, 3.23429399999999978, 5.824570000000000025e-01, 1.672568999999999972, 1.695751999999999926, 6.242940000000000156e-01, 0., 2.807907999999999848, 1.223827999999999916, 1.080135999999999985, 2.103319999999999912e-01, 2.361989999999999923e-01, 2.573360000000000092e-01},
			{4.250929999999999986e-01, 5.345509999999999984e-01, 1.239539999999999947e-01, 3.639700000000000157e-01, 5.272199999999999803e-02, 3.901919999999999833e-01, 2.426600999999999786, 1.269909999999999928e-01, 6.326067000000000107, 3.018480000000000052e-01, 4.8413299999999998e-01, 7.518780000000000463e-01, 3.32533000000000023e-01, 2.807907999999999848, 0., 8.581509999999999971e-01, 5.789870000000000294e-01, 1.708870000000000111e-01, 5.936069999999999958e-01, 3.144399999999999973e-01},
			{4.727181999999999995, 2.784478000000000009, 1.240275000000000016, 6.119729999999999892e-01, 3.618190000000000017e-01, 1.739989999999999926, 9.900120000000000031e-01, 6.410499999999999532e-02, 7.486829999999999874e-01, 1.822870000000000046e-01, 3.469599999999999906e-01, 4.00835800000000031, 1.338132000000000099, 1.223827999999999916, 8.581509999999999971e-01, 0., 6.472279000000000337, 9.836899999999999811e-02, 2.488619999999999999e-01, 4.005469999999999864e-01},
			{2.139501000000000097, 1.143480000000000052, 4.258600000000000163e-01, 6.04544999999999999e-01, 1.650010000000000088e-01, 1.29836000000000007e-01, 5.842619999999999481e-01, 1.033738999999999963, 1.136862999999999957, 3.029359999999999831e-01, 2.020366000000000106, 2.000678999999999874, 5.71467999999999976e-01, 1.080135999999999985, 5.789870000000000294e-01, 6.472279000000000337, 0., 2.188158000000000047, 1.408250000000000057e-01, 2.45841000000000004e-01},
			{2.547870000000000079, 1.959290999999999894, 3.796700000000000075e-02, 2.450340000000000018e-01, 6.546830000000000149e-01, 7.670100000000000529e-02, 1.190129999999999938e-01, 6.491069999999999895e-01, 1.852020000000000055e-01, 1.702744999999999953, 1.898717999999999906, 8.36879999999999985e-02, 2.965010000000000145e-01, 2.103319999999999912e-01, 1.708870000000000111e-01, 9.836899999999999811e-02, 2.188158000000000047, 0., 1.895100000000000118e-01, 2.493130000000000068e-01},
			{1.807169999999999888e-01, 6.701279999999999459e-01, 2.98899999999999999e-02, 7.785200000000000453e-02, 2.457120999999999889, 2.684909999999999797e-01, 5.970539999999999736e-01, 1.116599999999999954e-01, 4.990599999999999897e-02, 6.196319999999999606e-01, 6.961749999999999883e-01, 4.537599999999999967e-02, 9.513099999999999334e-02, 2.361989999999999923e-01, 5.936069999999999958e-01, 2.488619999999999999e-01, 1.408250000000000057e-01, 1.895100000000000118e-01, 0., 3.151815000000000033},
			{2.18958999999999987e-01, 1.165532000000000012, 1.351070000000000049e-01, 1.200370000000000048e-01, 7.803901999999999894, 5.467899999999999844e-02, 5.306834000000000273, 2.325230000000000075e-01, 1.319319999999999937e-01, 2.996480000000000254e-01, 4.813060000000000116e-01, 6.120250000000000412e-01, 8.961299999999999821e-02, 2.573360000000000092e-01, 3.144399999999999973e-01, 4.005469999999999864e-01, 2.45841000000000004e-01, 2.493130000000000068e-01, 3.151815000000000033, 0.}};

models.protein.empirical.JTT.empirical_R = {{0., 5.744780000000000442e-01, 8.274449999999999861e-01, 1.06668099999999999, 1.382929999999999993e-01, 1.740159000000000011, 2.199699999999999989e-01, 3.616840000000000055e-01, 3.694370000000000154e-01, 3.10006999999999977e-01, 4.693950000000000067e-01, 5.579669999999999908e-01, 1.959599000000000091, 5.567250000000000254e-01, 5.316779999999999839e-01, 3.887094999999999967, 4.582564999999999777, 2.924160999999999788, 8.432900000000000118e-02, 1.394920000000000049e-01},
			{5.744780000000000442e-01, 0., 1.056249999999999967e-01, 5.390699999999999659e-02, 6.783350000000000213e-01, 5.463890000000000136e-01, 7.249980000000000313e-01, 1.505589999999999984e-01, 4.900899999999999701e-02, 1.645929999999999893e-01, 4.092020000000000102e-01, 3.133110000000000062e-01, 1.23652999999999999e-01, 9.130399999999999627e-02, 1.019843000000000055, 2.155330999999999886, 4.698229999999999906e-01, 6.213229999999999587e-01, 1.104181000000000079, 2.114851999999999954},
			{8.274449999999999861e-01, 1.056249999999999967e-01, 0., 7.766556999999999711, 3.25220000000000023e-02, 1.272434000000000065, 1.032342000000000093, 1.159680000000000016e-01, 2.824659999999999949e-01, 6.148599999999999899e-02, 1.900010000000000032e-01, 5.549529999999999852, 1.271639999999999993e-01, 5.216460000000000541e-01, 1.548990000000000089e-01, 5.89268000000000014e-01, 4.251590000000000091e-01, 3.152610000000000134e-01, 5.746600000000000319e-02, 4.539520000000000222e-01},
			{1.06668099999999999, 5.390699999999999659e-02, 7.766556999999999711, 0., 4.38289999999999999e-02, 1.115631999999999957, 2.437680000000000125e-01, 1.117729999999999974e-01, 1.731684000000000001, 9.748500000000000221e-02, 1.750839999999999896e-01, 5.781150000000000455e-01, 1.91993999999999998e-01, 3.417705999999999911, 3.184830000000000161e-01, 3.124489999999999768e-01, 3.315839999999999899e-01, 4.652709999999999901e-01, 1.143809999999999966e-01, 6.345199999999999452e-02},
			{1.382929999999999993e-01, 6.783350000000000213e-01, 3.25220000000000023e-02, 4.38289999999999999e-02, 0., 5.021199999999999969e-02, 4.534279999999999977e-01, 7.770899999999999475e-01, 2.452100000000000113e-02, 2.500293999999999794, 4.361809999999999854e-01, 7.348100000000000465e-02, 1.484830000000000039e-01, 4.568300000000000138e-02, 6.531399999999999706e-02, 9.439710000000000045e-01, 1.389039999999999997e-01, 5.9347799999999995e-01, 5.379220000000000113e-01, 5.484236000000000111},
			{1.740159000000000011, 5.463890000000000136e-01, 1.272434000000000065, 1.115631999999999957, 5.021199999999999969e-02, 0., 2.016959999999999864e-01, 5.376899999999999735e-02, 2.698400000000000243e-01, 6.949199999999999822e-02, 1.303789999999999949e-01, 7.733130000000000281e-01, 2.080809999999999882e-01, 2.312939999999999996e-01, 1.359652000000000083, 1.874295999999999962, 3.16861999999999977e-01, 4.701400000000000023e-01, 5.441799999999999971e-01, 5.249999999999999806e-02},
			{2.199699999999999989e-01, 7.249980000000000313e-01, 1.032342000000000093, 2.437680000000000125e-01, 4.534279999999999977e-01, 2.016959999999999864e-01, 0., 1.817880000000000051e-01, 5.250960000000000072e-01, 5.405710000000000237e-01, 3.296600000000000086e-01, 4.025777999999999857, 1.141961000000000004, 5.684079999999999799, 3.210671000000000053, 7.434579999999999522e-01, 4.773549999999999738e-01, 1.218270000000000047e-01, 1.281930000000000014e-01, 5.848399999999999821},
			{3.616840000000000055e-01, 1.505589999999999984e-01, 1.159680000000000016e-01, 1.117729999999999974e-01, 7.770899999999999475e-01, 5.376899999999999735e-02, 1.817880000000000051e-01, 0., 2.025619999999999921e-01, 2.335138999999999854, 4.831666000000000238, 4.910030000000000228e-01, 9.858000000000000096e-02, 7.827000000000000624e-02, 2.39194999999999991e-01, 4.051190000000000069e-01, 2.553805999999999798, 9.533943000000000723, 1.345099999999999907e-01, 3.034450000000000203e-01},
			{3.694370000000000154e-01, 4.900899999999999701e-02, 2.824659999999999949e-01, 1.731684000000000001, 2.452100000000000113e-02, 2.698400000000000243e-01, 5.250960000000000072e-01, 2.025619999999999921e-01, 0., 1.464810000000000001e-01, 6.245810000000000528e-01, 2.529516999999999793, 2.163450000000000095e-01, 2.966731999999999925, 6.529255000000000031, 4.744780000000000109e-01, 9.656409999999999716e-01, 1.240659999999999957e-01, 8.913400000000000489e-02, 8.790399999999999603e-02},
			{3.10006999999999977e-01, 1.645929999999999893e-01, 6.148599999999999899e-02, 9.748500000000000221e-02, 2.500293999999999794, 6.949199999999999822e-02, 5.405710000000000237e-01, 2.335138999999999854, 1.464810000000000001e-01, 0., 3.856905999999999946, 1.372889999999999944e-01, 1.060503999999999891, 7.090039999999999676e-01, 3.722610000000000086e-01, 5.925110000000000099e-01, 2.725139999999999785e-01, 1.761438999999999977, 5.303240000000000176e-01, 2.410940000000000027e-01},
			{4.693950000000000067e-01, 4.092020000000000102e-01, 1.900010000000000032e-01, 1.750839999999999896e-01, 4.361809999999999854e-01, 1.303789999999999949e-01, 3.296600000000000086e-01, 4.831666000000000238, 6.245810000000000528e-01, 3.856905999999999946, 0., 3.30720000000000014e-01, 1.642149999999999999e-01, 4.569010000000000016e-01, 4.310450000000000115e-01, 2.855639999999999845e-01, 2.114727999999999941, 3.03853300000000015, 2.01334000000000013e-01, 1.89870000000000011e-01},
			{5.579669999999999908e-01, 3.133110000000000062e-01, 5.549529999999999852, 5.781150000000000455e-01, 7.348100000000000465e-02, 7.733130000000000281e-01, 4.025777999999999857, 4.910030000000000228e-01, 2.529516999999999793, 1.372889999999999944e-01, 3.30720000000000014e-01, 0., 1.218039999999999956e-01, 7.688340000000000174e-01, 4.510950000000000237e-01, 5.057964000000000127, 2.351310999999999929, 1.645250000000000046e-01, 2.769999999999999893e-02, 7.006930000000000103e-01},
			{1.959599000000000091, 1.23652999999999999e-01, 1.271639999999999993e-01, 1.91993999999999998e-01, 1.484830000000000039e-01, 2.080809999999999882e-01, 1.141961000000000004, 9.858000000000000096e-02, 2.163450000000000095e-01, 1.060503999999999891, 1.642149999999999999e-01, 1.218039999999999956e-01, 0., 1.608125999999999944, 7.104890000000000372e-01, 2.788406000000000162, 1.176960999999999924, 2.115609999999999991e-01, 6.996499999999999941e-02, 1.138500000000000068e-01},
			{5.567250000000000254e-01, 9.130399999999999627e-02, 5.216460000000000541e-01, 3.417705999999999911, 4.568300000000000138e-02, 2.312939999999999996e-01, 5.684079999999999799, 7.827000000000000624e-02, 2.966731999999999925, 7.090039999999999676e-01, 4.569010000000000016e-01, 7.688340000000000174e-01, 1.608125999999999944, 0., 3.021994999999999987, 5.488070000000000448e-01, 5.238249999999999851e-01, 1.797709999999999864e-01, 1.72205999999999998e-01, 2.547449999999999992e-01},
			{5.316779999999999839e-01, 1.019843000000000055, 1.548990000000000089e-01, 3.184830000000000161e-01, 6.531399999999999706e-02, 1.359652000000000083, 3.210671000000000053, 2.39194999999999991e-01, 6.529255000000000031, 3.722610000000000086e-01, 4.310450000000000115e-01, 4.510950000000000237e-01, 7.104890000000000372e-01, 3.021994999999999987, 0., 1.00155100000000008, 6.502820000000000267e-01, 1.71995000000000009e-01, 1.257961000000000107, 2.356010000000000049e-01},
			{3.887094999999999967, 2.155330999999999886, 5.89268000000000014e-01, 3.124489999999999768e-01, 9.439710000000000045e-01, 1.874295999999999962, 7.434579999999999522e-01, 4.051190000000000069e-01, 4.744780000000000109e-01, 5.925110000000000099e-01, 2.855639999999999845e-01, 5.057964000000000127, 2.788406000000000162, 5.488070000000000448e-01, 1.00155100000000008, 0., 4.777646999999999977, 4.085320000000000062e-01, 3.109270000000000089e-01, 6.286079999999999446e-01},
			{4.582564999999999777, 4.698229999999999906e-01, 4.251590000000000091e-01, 3.315839999999999899e-01, 1.389039999999999997e-01, 3.16861999999999977e-01, 4.773549999999999738e-01, 2.553805999999999798, 9.656409999999999716e-01, 2.725139999999999785e-01, 2.114727999999999941, 2.351310999999999929, 1.176960999999999924, 5.238249999999999851e-01, 6.502820000000000267e-01, 4.777646999999999977, 0., 1.143979999999999997, 8.05560000000000026e-02, 2.01093999999999995e-01},
			{2.924160999999999788, 6.213229999999999587e-01, 3.152610000000000134e-01, 4.652709999999999901e-01, 5.9347799999999995e-01, 4.701400000000000023e-01, 1.218270000000000047e-01, 9.533943000000000723, 1.240659999999999957e-01, 1.761438999999999977, 3.03853300000000015, 1.645250000000000046e-01, 2.115609999999999991e-01, 1.797709999999999864e-01, 1.71995000000000009e-01, 4.085320000000000062e-01, 1.143979999999999997, 0., 2.396969999999999934e-01, 1.65473000000000009e-01},
			{8.432900000000000118e-02, 1.104181000000000079, 5.746600000000000319e-02, 1.143809999999999966e-01, 5.379220000000000113e-01, 5.441799999999999971e-01, 1.281930000000000014e-01, 1.345099999999999907e-01, 8.913400000000000489e-02, 5.303240000000000176e-01, 2.01334000000000013e-01, 2.769999999999999893e-02, 6.996499999999999941e-02, 1.72205999999999998e-01, 1.257961000000000107, 3.109270000000000089e-01, 8.05560000000000026e-02, 2.396969999999999934e-01, 0., 7.47889000000000026e-01},
			{1.394920000000000049e-01, 2.114851999999999954, 4.539520000000000222e-01, 6.345199999999999452e-02, 5.484236000000000111, 5.249999999999999806e-02, 5.848399999999999821, 3.034450000000000203e-01, 8.790399999999999603e-02, 2.410940000000000027e-01, 1.89870000000000011e-01, 7.006930000000000103e-01, 1.138500000000000068e-01, 2.547449999999999992e-01, 2.356010000000000049e-01, 6.286079999999999446e-01, 2.01093999999999995e-01, 1.65473000000000009e-01, 7.47889000000000026e-01, 0.}};

models.protein.empirical.JC.empirical_R = {{0., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.},
			{1., 0., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.},
			{1., 1., 0., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.},
			{1., 1., 1., 0., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.},
			{1., 1., 1., 1., 0., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.},
			{1., 1., 1., 1., 1., 0., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.},
			{1., 1., 1., 1., 1., 1., 0., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.},
			{1., 1., 1., 1., 1., 1., 1., 0., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.},
			{1., 1., 1., 1., 1., 1., 1., 1., 0., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.},
			{1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.},
			{1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 1., 1., 1., 1., 1., 1., 1., 1., 1.},
			{1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 1., 1., 1., 1., 1., 1., 1., 1.},
			{1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 1., 1., 1., 1., 1., 1., 1.},
			{1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 1., 1., 1., 1., 1., 1.},
			{1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 1., 1., 1., 1., 1.},
			{1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 1., 1., 1., 1.},
			{1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 1., 1., 1.},
			{1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 1., 1.},
			{1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 1.},
			{1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0.}};









