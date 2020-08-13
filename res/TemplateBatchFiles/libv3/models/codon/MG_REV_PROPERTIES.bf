RequireVersion ("2.5.7");

LoadFunctionLibrary("../codon.bf");
LoadFunctionLibrary("../DNA.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");
LoadFunctionLibrary("MG_REV.bf");
LoadFunctionLibrary("../protein.bf");

/** @module models.codon.MG_REV_PROP */

/**
 * @name models.codon.MG_REV_PROP.ModelDescription
 * @param {String} type
 * @param {String} code
 * @param {Dict}   properties
 *                 'name' -> dict {
                        'one letter amino-acid' -> property
                    }
 */



terms.model.residue_name_map = "residue_to_id_map";
terms.model.MG_REV_PROP.mean_prop = "mean_property_value";


terms.model.MG_REV_PROP.predefined = {
    "Atchley" : {
             "Factor I bipolar":{
               "A":-0.591,
               "C":-1.343,
               "D":1.05,
               "E":1.357,
               "F":-1.006,
               "G":-0.384,
               "H":0.336,
               "I":-1.239,
               "K":1.831,
               "L":-1.019,
               "M":-0.663,
               "N":0.945,
               "P":0.189,
               "Q":0.931,
               "R":1.538,
               "S":-0.228,
               "T":-0.032,
               "V":-1.337,
               "W":-0.595,
               "Y":0.26
              },
             "Factor II secondary structure":{
               "A":-1.302,
               "C":0.465,
               "D":0.302,
               "E":-1.453,
               "F":-0.59,
               "G":1.652,
               "H":-0.417,
               "I":-0.547,
               "K":-0.5610000000000001,
               "L":-0.987,
               "M":-1.524,
               "N":0.828,
               "P":2.081,
               "Q":-0.179,
               "R":-0.055,
               "S":1.399,
               "T":0.326,
               "V":-0.279,
               "W":0.008999999999999999,
               "Y":0.83
              },
             "Factor III volume":{
               "A":-0.733,
               "C":-0.862,
               "D":-3.656,
               "E":1.477,
               "F":1.891,
               "G":1.33,
               "H":-1.673,
               "I":2.131,
               "K":0.533,
               "L":-1.505,
               "M":2.219,
               "N":1.299,
               "P":-1.628,
               "Q":-3.005,
               "R":1.502,
               "S":-4.76,
               "T":2.213,
               "V":-0.544,
               "W":0.672,
               "Y":3.097
              },
             "Factor IV composition":{
               "A":1.57,
               "C":-1.02,
               "D":-0.259,
               "E":0.113,
               "F":-0.397,
               "G":1.045,
               "H":-1.474,
               "I":0.393,
               "K":-0.277,
               "L":1.266,
               "M":-1.005,
               "N":-0.169,
               "P":0.421,
               "Q":-0.503,
               "R":0.44,
               "S":0.67,
               "T":0.908,
               "V":1.242,
               "W":-2.128,
               "Y":-0.838
              },
             "Factor V charge":{
               "A":-0.146,
               "C":-0.255,
               "D":-3.242,
               "E":-0.837,
               "F":0.412,
               "G":2.064,
               "H":-0.078,
               "I":0.8159999999999999,
               "K":1.648,
               "L":-0.912,
               "M":1.212,
               "N":0.9330000000000001,
               "P":-1.392,
               "Q":-1.853,
               "R":2.897,
               "S":-2.647,
               "T":1.313,
               "V":-1.262,
               "W":-0.184,
               "Y":1.512
              }
     },
     "LCAP" : {
        "Chemical Composition" : {
            {-0.5868927454579467}
            {-0.5868460291700075}
            {-0.5868457923690357}
            {-0.586845591345447}
            {-0.5868454280017823}
            {1.442413127765356}
            {-0.02935869926260699}
            {0.4272890876194427}
            {-0.5868451756885411}
            {-0.3012025870591999}
            {0.2419262642781204}
            {0.6848994365857728}
            {1.312947120977898}
            {-0.1156419909884284}
            {1.384845017074097}
            {0.7270716263181436}
            {3.342731546150394}
            {-0.4010042423293907}
            {0.3415170202317384}
            {0.4699025925311998}
            },
        "Polarity" : {
            {-0.9956290326072353}
            {-1.107220007612595}
            {-0.9956292005349376}
            {-0.8100842857525447}
            {-0.7357356892564315}
            {0.4912036016809371}
            {0.04436400722996305}
            {0.2672975016276415}
            {0.08234475255665361}
            {-0.624379765568966}
            {0.9376066737521076}
            {0.9747140546900138}
            {1.383229010106464}
            {1.271792102330754}
            {1.904064473242525}
            {1.643233506555707}
            {-0.8846717344020271}
            {-0.9219033329616851}
            {0.9747140551029528}
            {0.4163765565466709}
            },

        "Volume" : {
            {1.257637282334762}
            {0.7687291690973896}
            {0.7687283894814801}
            {0.6290409911059951}
            {0.1401332510120709}
            {-1.070491916990394}
            {-1.058850814502372}
            {-0.3953354747308322}
            {-1.093771866791469}
            {1.350763381619076}
            {0.4195123000264476}
            {0.1634169621566906}
            {-0.5117403202399798}
            {0.9549794083529999}
            {-0.5583026469068318}
            {0.1168546775994032}
            {-0.535022037570798}
            {2.142318464323941}
            {1.071384899476407}
            {-1.745647766340383}
            },
        "Iso-electric Point" : {
            {-0.3057542421356372}
            {-0.02273844767848476}
            {-9.737035422433548e-05}
            {-0.1585867248775605}
            {-0.03405832632980638}
            {-0.192548376206011}
            {0.1583919144982974}
            {0.07914812431496633}
            {-0.01141750903277056}
            {-0.2038687480243035}
            {0.8885757410301605}
            {-0.2095283729063603}
            {-0.3453753360155064}
            {2.105545575730344}
            {-1.839706056109631}
            {-1.584990785596673}
            {-0.5378282523391487}
            {-0.07368049231084244}
            {2.682902758933077}
            {-0.02839872447370546}
            },
        "Hydropathy" : {
            {1.301968439295016}
            {1.636753002222122}
            {1.871099553721135}
            {1.000656143286975}
            {1.770671640107803}
            {0.09673165440139944}
            {-0.171098473395166}
            {0.1302104025093898}
            {0.9671789981061145}
            {-0.0706624429784429}
            {-0.7067573940958518}
            {-0.8071939952622403}
            {-0.807193639057381}
            {-0.9411095669881511}
            {-0.8071943355554064}
            {-0.8071934015198535}
            {1.201528247390539}
            {0.06325356334516005}
            {-1.14198047362342}
            {0.2306464200526206}
            }
     }
};


lfunction models.codon.MG_REV_PROP.ModelDescription(type, code, properties) {


    // piggyback on the standard MG_REV model for most of the code

    mg_base = models.codon.MG_REV.ModelDescription (type, code);
    mg_base[utility.getGlobalValue("terms.description")] = "The Muse-Gaut 94 codon-substitution model coupled with the general time reversible (GTR) model of nucleotide substitution, which allows incorporates amino-acid residue properties into the non-synonymous rates";
    mg_base[utility.getGlobalValue("terms.model.q_ij")] = "models.codon.MG_REV_PROP._GenerateRate";
    mg_base[utility.getGlobalValue("terms.model.residue_properties")] = models.codon.MG_REV_PROP._munge_properties(properties);
    mg_base[utility.getGlobalValue("terms.model.residue_name_map")] = parameters.ValidateIDs (utility.Keys (mg_base[utility.getGlobalValue("terms.model.residue_properties")]));
    mg_base[utility.getGlobalValue("terms.model.post_definition")] = "models.codon.MG_REV_PROP.post_definition";
    mg_base[utility.getGlobalValue("terms.model.set_branch_length")] = "models.codon.MG_REV_PROP.set_branch_length";

    return mg_base;
}

/**
 * @name models.codon.MG_REV_PROP._munge_properties
 * @param {String or AssociativeList} properties
 * @return {AssociativeList}
 * Convert a set of amino-acid properties into the expected format
 */


lfunction models.codon.MG_REV_PROP._munge_properties (properties) {

    if (Type (properties) == "String") {
        assert (^"terms.model.MG_REV_PROP.predefined" / properties, properties + " is not a valid predefined property set");
        return models.codon.MG_REV_PROP._munge_properties ((^"terms.model.MG_REV_PROP.predefined")[properties]);
    }

    assert (Type (properties) == "AssociativeList", "Properties definition must be an AssociativeArray or a String for predefined properties");

        /*
             each property must map to either a dictionary with '1 letter AA' -> value
             or to an array with 20 values
        */


    valid_aa = utility.MatrixToDict (^"models.protein.alphabet");

    for (key, value; in; properties) {
        if (Type (value) == "AssociativeList") {
            io.CheckAssertion ("utility.Array1D(`&value`)==20", "A dictionary of amino-acid properties must have 20 entries for " + key);
            for (_i = 0; _i < 20; _i += 1) {
                assert (value / (^"models.protein.alphabet")[_i],
                    "Property " + key + " is not defined for amino-acid " + (^"models.protein.alphabet")[_i]);
            }
        } else {
            if (Type (value) == "Matrix") {
                io.CheckAssertion ("utility.Array1D(`&value`)==20", "A matrix of amino-acid properties must have 20 entries for " + key);
                value.dict = {};
                for (_i = 0; _i < 20; _i += 1) {
                    value.dict [(^"models.protein.alphabet")[_i]] = value[_i];
                }
                properties[key] = value.dict;
            } else {
                assert (0, "Invalid entry for properties in class " + key);
            }
        }
    }

    assert (utility.Array1D(properties) > 0, "At least one valid property set must be defined");
    return properties;
}

lfunction models.codon.MG_REV_PROP._GenerateRate(fromChar, toChar, namespace, model_type, model) {
    return models.codon.MG_REV_PROP._GenerateRate_generic (fromChar, toChar, namespace, model_type,
        model[utility.getGlobalValue("terms.translation_table")],
        "alpha", utility.getGlobalValue("terms.parameters.synonymous_rate"),
        "beta",  utility.getGlobalValue("terms.parameters.nonsynonymous_rate"),
        "lambda",
        "omega",
        model[utility.getGlobalValue("terms.model.residue_properties")],
        model[utility.getGlobalValue("terms.model.residue_name_map")]
        );
}

/**
 * @name models.codon.MG_REV_PROP._GenerateRate
 * @param {Number} fromChar
 * @param {Number} toChar
 * @param {String} namespace
 * @param {String} model_type
 * @param {Matrix} _tt - translation table
 */


lfunction models.codon.MG_REV_PROP._GenerateRate_generic (fromChar, toChar, namespace, model_type, _tt, alpha, alpha_term, beta, beta_term, lambda, omega, properties, property_id_map) {

    _GenerateRate.p = {};
    _GenerateRate.diff = models.codon.diff.complete(fromChar, toChar);
    diff_count = utility.Array1D (_GenerateRate.diff);

    if (diff_count == 1) {

        _GenerateRate.p[model_type] = {};
        _GenerateRate.p[utility.getGlobalValue("terms.global")] = {};

        nuc_rate = "";

        for (i = 0; i < diff_count; i += 1) {
            if ((_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.from")] > (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.to")]) {
                nuc_p = "theta_" + (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.to")] + (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.from")];
            } else {
                nuc_p = "theta_" + (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.from")] +(_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.to")];
            }
            nuc_p = parameters.ApplyNameSpace(nuc_p, namespace);
            (_GenerateRate.p[utility.getGlobalValue("terms.global")])[terms.nucleotideRateReversible((_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.from")], (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.to")])] = nuc_p;

            nuc_rate = parameters.AppendMultiplicativeTerm (nuc_rate, nuc_p);
       }


        rate_entry = nuc_rate;

        if (_tt[fromChar] != _tt[toChar]) {
            if (model_type == utility.getGlobalValue("terms.global")) {
                aa_rate = {};
                prop_count = Abs (properties);
                prop_names = utility.Keys (properties);
                (_GenerateRate.p[model_type]) [^"terms.parameters.log_omega_ratio"] = parameters.ApplyNameSpace (omega, namespace);

                for (prop_name, prop_values; in; properties) {
                    prop_diff = Abs (prop_values[_tt[fromChar]] - prop_values[_tt[toChar]]);
                    term_rate = parameters.ApplyNameSpace(lambda + "_" + property_id_map[prop_name], namespace);
                    (_GenerateRate.p[model_type])[terms.propertyImportance (prop_name)] = term_rate;
                    aa_rate + ("`term_rate`*" + prop_diff);
                }

                //aa_rate = parameters.ApplyNameSpace(lambda, namespace);
                rate_entry += "*Min(10000,Exp(" + (_GenerateRate.p[model_type]) [^"terms.parameters.log_omega_ratio"] + "-(" + Join("+",aa_rate) + ")))";
             } else {
                aa_rate = {};
                prop_count = Abs (properties);
                prop_names = utility.Keys (properties);
                (_GenerateRate.p[model_type]) [beta_term] = beta;

                for (prop_name, prop_values; in; properties) {
                    prop_diff = Abs (prop_values[_tt[fromChar]] - prop_values[_tt[toChar]]);
                    term_rate = lambda + "_" + property_id_map[prop_name];
                    (_GenerateRate.p[model_type])[terms.propertyImportance (prop_name)] = term_rate;
                    aa_rate + ("`term_rate`*" + prop_diff);
                }
                rate_entry += "*Exp(-(" + Join("+",aa_rate) + "))";
                (_GenerateRate.p[model_type])[alpha_term] = alpha;
                rate_entry = "Min(10000," + rate_entry  + "*" + beta + ")";
             }
        } else {
            if (model_type == utility.getGlobalValue("terms.local")) {
                 (_GenerateRate.p[model_type])[alpha_term] = alpha;
                rate_entry += "*" + alpha;
           } else {
                _GenerateRate.p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate;
            }
        }

        _GenerateRate.p[utility.getGlobalValue("terms.model.rate_entry")] = rate_entry;
    }
    return _GenerateRate.p;
}

lfunction  models.codon.MG_REV_PROP.post_definition (model) {
    prop_range = {
        ^"terms.lower_bound": "-10",
        ^"terms.upper_bound": "10"
    };

    for (id; in ; model.GetParameters_RegExp (model, terms.propertyImportance ('') + "|" + ^"terms.parameters.log_omega_ratio")) {
        parameters.SetRange(id, prop_range);
        parameters.SetValue(id, 0.1);
    }

    models.generic.post.definition (model);
}

lfunction models.codon.MG_REV_PROP.set_branch_length(model, value, parameter) {

    if (utility.Has (model, ^"terms.model.MG_REV_PROP.mean_prop", "Number")) {
      properties = model.GetLocalParameters_RegExp(model, terms.propertyImportance (""));
      for (tag, id; in; properties) {
        parameters.SetValue (id, model[^"terms.model.MG_REV_PROP.mean_prop"]);
      }
      if (utility.Has (model, "fraction_same", "Number")) {
        fs = model["fraction_same"];
        models.codon.MG_REV.set_branch_length(model,value,parameter);
        weighted = model ["models.codon.MG_REV.set_branch_length"];
        if (Type (weighted) == "AssociativeList") {
          for (tag, id; in; properties) {
            parameters.SetValue (id, 0.);
          }
          models.codon.MG_REV.set_branch_length(model,value,parameter);
          unweighted = model ["models.codon.MG_REV.set_branch_length"];
          for (tag, value; in; weighted) {
            //console.log (fs * unweighted[tag] + (1-fs) * value);
            parameters.SetValue (tag, fs * unweighted[tag] + (1-fs) * value);
          }
          return 0;
        }
      }
    }
    return models.codon.MG_REV.set_branch_length(model,value,parameter);

}
