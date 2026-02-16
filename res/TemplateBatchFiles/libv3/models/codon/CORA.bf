RequireVersion ("2.5.90");

LoadFunctionLibrary("../codon.bf");
LoadFunctionLibrary("../DNA.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");
LoadFunctionLibrary("MG_REV.bf");
LoadFunctionLibrary("../protein.bf");

/** @module models.codon.CORA */

/**
 * @name models.codon.CORA.ModelDescription
 * @param {String} type
 * @param {String} code
 * @param {Dict}   classes
 *                 'one letter amino-acid' -> [property list]
 */



terms.models.residues = "aa_residues";
terms.models.conservative_type = "is_conservative";


terms.model.CORA.predefined = {
    "Weber and Whelan" : {
    "A": {{
        "Non-polar",
        "Small"
    }},
    "C": {{
        "Polar",
        "Small"
    }},
    "D": {{
        "Polar",
        "Small"
    }},
    "E": {{
        "Polar",
        "Large"
    }},
    "F": {{
        "Non-polar",
        "Large"
    }},
    "G": {{
        "Non-polar",
        "Small"
    }},
    "H": {{
        "Polar",
        "Large"
    }},
    "I": {{
        "Non-polar",
        "Large"
    }},
    "K": {{
        "Polar",
        "Large"
    }},
    "L": {{
        "Non-polar",
        "Large"
    }},
    "M": {{
        "Non-polar",
        "Large"
    }},
    "N": {{
        "Polar",
        "Small"
    }},
    "P": {{
        "Non-polar",
        "Small"
    }},
    "Q": {{
        "Polar",
        "Large"
    }},
    "R": {{
        "Polar",
        "Large"
    }},
    "S": {{
        "Polar",
        "Small"
    }},
    "T": {{
        "Polar",
        "Small"
    }},
    "V": {{
        "Non-polar",
        "Small"
    }},
    "W": {{
        "Polar",
        "Large"
    }},
    "Y": {{
        "Polar",
        "Large"
    }}
    }
};

//----------------------------------------------------------------------------------------------------------------


lfunction model.codon.CORA.prompt_and_define (type, code) {
    KeywordArgument ("attributes", "How to decide if a substitution is conservative or radical", "Weber and Whelan");
    
    attribute_set = io.SelectAnOption (
            {
                "Weber and Whelan":"Use conservative and radical substitution definitions based on polarity and size from Weber and Whelan [doi.org/10.1093/molbev/msz003]",
                "Custom":    "Load the set of amino-acid attributes from a file"
            }, 
            "How to decide if a substitution is conservative or radica");

    
    
    if (property_set == "Custom") {
        KeywordArgument ("attribute-file", "JSON file which defines amino-acid attributes");
        attribute_set = io.PromptUserForFilePathRead ("JSON file which defines amino-acid attributes");
        attribute_set = io.ParseJSON(attribute_set);
        console.log (">Loaded a set of `Abs(attribute_set)` properties");
     }
    
     return models.codon.CORA.ModelDescription(type, code, attribute_set);
}


lfunction models.codon.CORA.ModelDescription(type, code, properties) {


    // piggyback on the standard MG_REV model for most of the code

    mg_base = models.codon.MG_REV.ModelDescription (type, code);
    mg_base[utility.getGlobalValue("terms.description")] = "The Muse-Gaut 94 codon-substitution model coupled with the general time reversible (GTR) model of nucleotide substitution, which partitions non-synonymous substitutions in conservative and radical";
    mg_base[utility.getGlobalValue("terms.model.q_ij")] = "models.codon.CORA._GenerateRate";
    mg_base[utility.getGlobalValue("terms.model.residue_properties")] = models.codon.CORA._munge_properties(properties);
   
        
    return mg_base;
}

/**
 * @name models.codon.CORA._munge_properties
 * @param {String or AssociativeList} properties
 * @return {AssociativeList}
 * Convert a set of amino-acid properties into the expected format
 * This will also recenter the properties (mean 0, variance of 1)
 */


lfunction models.codon.CORA._munge_properties (properties) {

    if (Type (properties) == "String") {
        assert (^"terms.model.CORA.predefined" / properties, properties + " is not a valid predefined attribute set");
        return models.codon.CORA._munge_properties ((^"terms.model.CORA.predefined")[properties]);
    }

    assert (Type (properties) == "AssociativeList", "Properties definition must be an AssociativeArray or a String for predefined properties");

        /*
             each property must map to either a dictionary with '1 letter AA' -> value
             or to an array with 20 values
        */


    valid_aa = utility.MatrixToDict (^"models.protein.alphabet");
    attribute_count = -1;
    
    prop_set = {};

    for (key, value; in; properties) {
        if (Type (value) == "Matrix") {
            a_count = utility.Array1D(value);
            if (attribute_count < 0) {
                io.CheckAssertion ("`&a_count` >  0", "Must have >=1 attributes defined");
                attribute_count = a_count;
            } else {
               io.CheckAssertion ("`&a_count` == `&attribute_count`", "All amino-acids must have the same number of attributes");
            }
            
            prop_set[key] = value;
        } else {
            assert (0, "Invalid entry for properties for amino-acid " + key);
        }
    }
    
    
    
    io.CheckAssertion ("Abs(`&prop_set`)==20", "A dictionary of amino-acid attributes must have 20 entries");
        for (_i = 0; _i < 20; _i += 1) {
            assert (prop_set / (^"models.protein.alphabet")[_i],
                "Attribues are not defined for amino-acid " + (^"models.protein.alphabet")[_i]);
        }
        
    
    
    keys = Rows (prop_set);
    type = {20,20};
    
    key_map = {};
    
    for (i = 0; i < 20; i+=1) {
        key_map [keys[i]] = i;
        for (j = i+1; j < 20; j+=1) {
            conservative = 1;
            for (k = 0; k < attribute_count && conservative; k+=1) {
                conservative = (prop_set[keys[i]])[k] == (prop_set[keys[j]])[k];
            }
            type [i][j] =  conservative;
            type [j][i] =  conservative;
        }
    }    
    
    //assert (0);
    
    prop_set = {
            ^"terms.models.residues" : key_map,
            ^"terms.models.conservative_type" : type
     };
     
     
    return prop_set;
}

lfunction models.codon.CORA._GenerateRate(fromChar, toChar, namespace, model_type, model) {
    return models.codon.CORA._GenerateRate_generic (fromChar, toChar, namespace, model_type,
        model[utility.getGlobalValue("terms.translation_table")],
        "alpha", utility.getGlobalValue("terms.parameters.synonymous_rate"),
        "beta",  utility.getGlobalValue("terms.parameters.nonsynonymous_rate"),
        "omega", utility.getGlobalValue("terms.parameters.log_omega_ratio"),
        model[utility.getGlobalValue("terms.model.residue_properties")]
        );
}

/**
 * @name models.codon.CORA._GenerateRate
 * @param {Number} fromChar
 * @param {Number} toChar
 * @param {String} namespace
 * @param {String} model_type
 * @param {Matrix} _tt - translation table
 */


lfunction models.codon.CORA._GenerateRate_generic (fromChar, toChar, namespace, model_type, _tt, alpha, alpha_term, beta, beta_term, omega, omega_term, properties) {

 
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
            is_conservative = (properties[^"terms.models.conservative_type"])
                                [(properties[^"terms.models.residues"])[_tt[fromChar] ]]
                                [(properties[^"terms.models.residues"])[_tt[toChar] ]];
                                
            if (model_type == utility.getGlobalValue("terms.global")) {
            
                if (is_conservative) {
                    term_rate = parameters.ApplyNameSpace (omega + "_conservative", namespace);
                    (_GenerateRate.p[model_type]) [terms.omegaType (omega_term, "conservative substitutions") ] = term_rate;
                } else {
                    term_rate = parameters.ApplyNameSpace (omega + "_radical", namespace);
                    (_GenerateRate.p[model_type]) [terms.omegaType (omega_term, "radical substitutions") ] = term_rate;
                }
                


                rate_entry =  parameters.AppendMultiplicativeTerm (rate_entry, term_rate);
             } else {
                if (is_conservative) {
                    term_rate = beta + "_conservative";
                    (_GenerateRate.p[model_type]) [terms.omegaType (beta_term, "conservative substitutions") ] = term_rate;
                } else {
                    term_rate = beta + "_radical";
                    (_GenerateRate.p[model_type]) [terms.omegaType (beta_term, "radical substitutions") ] = term_rate;
                }
                

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
    
    //console.log (_GenerateRate.p);
    
    return _GenerateRate.p;
}

