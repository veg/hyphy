LoadFunctionLibrary("SelectionAnalyses/modules/io_functions.ibf");
LoadFunctionLibrary("MG_REV.bf");

/** @module models.codon.MG_REV */

/**
 * @name models.codon.MSS.ModelDescription
 * @param {String} type
 * @param {String} code
 * @param {Dict} codon_classes : codon => neutral class
  */
  
terms.model.MSS.syn_rate_within  = " within codon class ";
terms.model.MSS.syn_rate_between = " between codon classes ";
terms.model.MSS.between = "synonymous rate between codon classes";
terms.model.MSS.neutral = "neutral reference";
terms.model.MSS.codon_classes = "codon classes";
terms.model.MSS.normalize = "normalize rates";

//----------------------------------------------------------------------------------------------------------------


lfunction model.codon.MSS.prompt_and_define_freq (type, code, freq) {
    KeywordArgument ("mss-type", "How to partition synonymous codons into classes", "File");

    utility.ExecuteInGlobalNamespace ("fitter.frequency_type = '`freq`'");
    partitioning_option = io.SelectAnOption (
    {
        {"Full", "Each set of codons mapping to the same amino-acid class have a separate substitution rate (Valine == neutral)"}
        {"SynREV", "Each set of codons mapping to the same amino-acid class have a separate substitution rate (mean = 1)"}
        {"SynREV2", "Each pair of synonymous codons mapping to the same amino-acid class and separated by a transition have a separate substitution rate (no rate scaling))"}
        {"SynREV2g", "Each pair of synonymous codons mapping to the same amino-acid class and separated by a transition have a separate substitution rate (Valine == neutral). All between-class synonymous substitutions share a rate."}
        {"SynREVCodon", "Each codon pair that is exchangeable gets its own substitution rate (fully estimated, mean = 1)"}
        {"Random", "Random partition (specify how many classes; largest class = neutral)"}
        {"File", "Load a TSV partition from file (prompted for neutral class)"}
    },
    "Synonymous Codon Class Definitions");
    
    gc = models.codon.MapCode (code);
    
    if (partitioning_option == "Full") {
        bins    = {};
        mapping = {};
        tt = gc [^"terms.translation_table"];
        syn_codons = {};
        for (codon; in; gc [^"terms.sense_codons"]) {
            t = tt[codon];
            if (syn_codons / t == FALSE) {
                syn_codons[t] = {};
            }
            syn_codons[t] + codon;
        }
        
        nt = null;
        for (aa, codons; in; syn_codons) {
            if (Abs (codons) > 1) {
                for (c; in; codons) {
                    mapping[c] = aa+"_"+c;
                    if (aa == "V" && None == nt) {
                        nt = mapping[c];
                    }
                }
            }
        }
        
        return  models.codon.MSS.ModelDescription(type, code,
           {^"terms.model.MSS.codon_classes" : mapping, ^"terms.model.MSS.neutral" : "V"}
        );
    }
    
    if (partitioning_option == "SynREV" || partitioning_option == "SynREVFull" || partitioning_option == "SynREVCodon" ) {
        bins          = {};
        mapping       = {};
        mapping_codon = {};
        tt = gc [^"terms.translation_table"];
        for (codon; in; gc [^"terms.sense_codons"]) {
            mapping[codon] = tt[codon];
            mapping_codon [codon] = codon;
            bins[mapping[codon]] += 1;
        }
        if (partitioning_option == "SynREV") {
            return  models.codon.MSS.ModelDescription(type, code,
               {
                    ^"terms.model.MSS.codon_classes" : mapping, 
                    ^"terms.model.MSS.normalize" : TRUE
                    //^"terms.model.MSS.neutral" : "V"
               }
            );
        }
        if (partitioning_option == "SynREVFull") {
            return  models.codon.MSS.ModelDescription(type, code,
               {^"terms.model.MSS.codon_classes" : mapping}
            );
        }
        return  models.codon.MSS.ModelDescription(type, code,
           {^"terms.model.MSS.codon_classes" : mapping_codon, ^"terms.model.MSS.normalize" : TRUE}
        );
   }
   
   if (partitioning_option == "SynREV2" || partitioning_option == "SynREV2g") {
        bins    = {};
        aa_mapping = {};
        tt = gc [^"terms.translation_table"];
        
        for (codon; in; gc [^"terms.sense_codons"]) {
            AA = tt[codon];
            if (aa_mapping / AA == FALSE) {
                aa_mapping[AA] = {};
            }
            aa_mapping[AA] + codon;
        }
        
        neutral = "";
        mapping = {};
        
        for (aa, codons; in; aa_mapping) {
            N = Abs (codons);
            if (N >=2) {
                subsets = {};
                for (c; in; codons) {
                    matched = FALSE;
                    for (s; in; subsets) {
                        if (Abs (s) == 1) {
                            if ((s[0])[0] == c[0] && (s[0])[1] == c[1]) {
                                l3 = (s[0])[2];
                                if (l3 == 'A' && c[2] == 'G' || l3 == 'G' && c[2] == 'A' || l3 == 'C' && c[2] == 'T' || l3 == 'T' && c[2] == 'C') {
                                    s + c;
                                    matched = TRUE;
                                }
                            }
                        }   
                    }
                    if (!matched) {
                        subsets + {"0" : c__};
                    }
                }
                for (s; in; subsets) {
                    tag = aa + "_" + Join ("_", s);
                    for (c; in; s) {
                        mapping [c] = tag;
                    }
                    if (aa == "V" && Abs (neutral) == 0) {
                        neutral = tag;
                    }
                }
                
            } else {
                mapping[codons[0]] = aa;
            }
        }
        
        
        if (partitioning_option == "SynREV2g") {
            return  models.codon.MSS.ModelDescription(type, code,
               {^"terms.model.MSS.codon_classes" : mapping,  ^"terms.model.MSS.between" : "BETWEEN_CLASS"}
            );
        } else {
            return  models.codon.MSS.ModelDescription(type, code,
               {^"terms.model.MSS.codon_classes" : mapping}
            );
        }
        
   }
    
    if (partitioning_option == "Random") {
        KeywordArgument ("mss-classes", "How many codon rate classes");
        rc = io.PromptUser ("How many codon rate classes", 2,2,10,TRUE);
        bins    = {};
        mapping = {};
        for (codon; in; gc [^"terms.sense_codons"]) {
            mapping[codon] = "Class_" + Random (0,rc)$1;
            bins[mapping[codon]] += 1;
        }
        return  models.codon.MSS.ModelDescription(type, code,
           {^"terms.model.MSS.codon_classes" : mapping, ^"terms.model.MSS.neutral" : (Max(bins,1)["key"])}
        );
    }
    
    if (partitioning_option == "File") {
        KeywordArgument ("mss-file", "File defining the model partition");
        KeywordArgument ("mss-neutral", "Designation for the neutral substitution rate");
        return  models.codon.MSS.ModelDescription(type, code, models.codon.MSS.LoadClasses (null));
    }
    
    return {};
}

//----------------------------------------------------------------------------------------------------------------

lfunction model.codon.MSS.prompt_and_define (type, code) {
    return model.codon.MSS.prompt_and_define_freq (type, code, 'CF3x4');
}

//----------------------------------------------------------------------------------------------------------------
  
lfunction models.codon.MSS.ModelDescription(type, code, codon_classes) {


    m = Call ("models.codon.MG_REV.ModelDescription", type, code);
    
    if (^"fitter.frequency_type" == "F3x4") {
        m[utility.getGlobalValue("terms.model.frequency_estimator")] = "frequencies.empirical.F3x4";
    } else {
        if (^"fitter.frequency_type" == "F1x4") {
            m[utility.getGlobalValue("terms.model.frequency_estimator")] = "frequencies.empirical.F1x4";
        }
    }
    m[utility.getGlobalValue("terms.description")] = "The Muse-Gaut 94 codon-substitution model coupled with the general time reversible (GTR) model of nucleotide substitution, which allows multiple classes of synonymous substitution rates";
    m[utility.getGlobalValue("terms.model.q_ij")] = "models.codon.MSS._GenerateRate";
    m[utility.getGlobalValue("terms.model.MSS.codon_classes")] = codon_classes [^"terms.model.MSS.codon_classes"];
    m[utility.getGlobalValue("terms.model.MSS.neutral")] = codon_classes [^"terms.model.MSS.neutral"];
    if (codon_classes/utility.getGlobalValue("terms.model.MSS.between")) {
        m[^"terms.model.MSS.between"] = codon_classes [^"terms.model.MSS.between"];
    }
    
    
    if (codon_classes[utility.getGlobalValue("terms.model.MSS.normalize")]) {
        m[utility.getGlobalValue("terms.model.post_definition")] = "models.codon.MSS.post_definition";
    }
    
    return m;
}

//----------------------------------------------------------------------------------------------------------------

lfunction models.codon.MSS.post_definition (model) {
// TBD
    rates = model.GetParameters_RegExp (model,"^" + utility.getGlobalValue ("terms.parameters.synonymous_rate"));
    D = utility.Array1D (rates);
    w = 1 / D;
    
    vars = {1,D};
    weights = {1,D};
    terms = {1,D};
    
    
    i = 0;
    for (t,n; in; rates) {
        terms[i] = t;
        vars[i] = n;
        weights[i] = w;
        i += 1;
    } 
    ((model[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.global")]) * (parameters.ConstrainMeanOfSetWithTerms (vars, weights, terms, 1, model[^"terms.id"])[utility.getGlobalValue("terms.global")]);
    model = models.generic.post.definition (model);
}

//----------------------------------------------------------------------------------------------------------------

lfunction models.codon.MSS._GenerateRate (fromChar, toChar, namespace, model_type, model) {

    _GenerateRate.p = {};
    _GenerateRate.diff = models.codon.diff.complete(fromChar, toChar);
    diff_count = utility.Array1D (_GenerateRate.diff);

    omega_term = utility.getGlobalValue ("terms.parameters.omega_ratio");
    alpha_term = utility.getGlobalValue ("terms.parameters.synonymous_rate");
    beta_term  = utility.getGlobalValue ("terms.parameters.nonsynonymous_rate");
    nr = model[utility.getGlobalValue("terms.model.MSS.neutral")];
    omega      = "omega";
    alpha      = "alpha";
    beta       = "beta";
    between_rate = model[^"terms.model.MSS.between"];

    _tt = model[utility.getGlobalValue("terms.translation_table")];

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
                aa_rate = parameters.ApplyNameSpace(omega, namespace);
                (_GenerateRate.p[model_type])[omega_term] = aa_rate;
            } else {
                aa_rate = beta;
                (_GenerateRate.p[model_type])[beta_term] = aa_rate;
            }
            rate_entry += "*" + aa_rate;
        } else {

            class_from = (model[^"terms.model.MSS.codon_classes"])[fromChar];
            class_to   = (model[^"terms.model.MSS.codon_classes"])[toChar];
            
            
            if ((Abs (class_from) && Abs (class_to)) == FALSE) {
                class_from = (model[^"terms.model.MSS.codon_classes"])[fromChar+toChar];
                class_to = class_from;
            }
            
            assert (Abs (class_from) && Abs (class_to), "The neutral class for `fromChar` to `toChar` is not specified in the model definition");

            if (class_from == class_to) {
                if (class_from == nr) {
                    if (model_type == utility.getGlobalValue("terms.local")) {
                        codon_rate = alpha + "_" + class_from;
                        (_GenerateRate.p[model_type])[alpha_term + ^"terms.model.MSS.syn_rate_within" + class_from] = codon_rate;
                        rate_entry += "*" + codon_rate;
                    } else {
                        rate_entry = nuc_rate;
                    }
                } else {
                    if (model_type == utility.getGlobalValue("terms.local")) {
                        codon_rate = alpha + "_" + class_from;
                    } else {
                        codon_rate = parameters.ApplyNameSpace(alpha + "_" + class_from, namespace);
                    }
                    (_GenerateRate.p[model_type])[alpha_term + ^"terms.model.MSS.syn_rate_within" + class_from] = codon_rate;
                    rate_entry += "*" + codon_rate;
                }
            } else {
                if (class_from > class_to) {
                    codon_rate = class_to;
                    class_to = class_from;
                    class_from = codon_rate;
                }
                if (Abs (between_rate)) {
                    if (model_type == utility.getGlobalValue("terms.local")) {
                        codon_rate = between_rate;
                    } else {
                        codon_rate = parameters.ApplyNameSpace(between_rate, namespace);
                    }
                    (_GenerateRate.p[model_type])[^"terms.model.MSS.between"] = codon_rate;
            
                } else {
                    if (class_from + class_to == nr) {
                        //console.log ("NEUTRAL");
                        codon_rate  = 1;
                    } else {
                        if (model_type == utility.getGlobalValue("terms.local")) {
                            codon_rate = alpha + "_" + class_from + "_" + class_to;
                        } else {
                            codon_rate = parameters.ApplyNameSpace(alpha + "_" + class_from + "_" + class_to, namespace);
                        }
                        (_GenerateRate.p[model_type])[alpha_term + ^"terms.model.MSS.syn_rate_between" + class_from + " and "  + class_to] = codon_rate;
                    }
                }
                rate_entry += "*" + codon_rate;
            }
        }

        _GenerateRate.p[utility.getGlobalValue("terms.model.rate_entry")] = rate_entry;
     }

    return _GenerateRate.p;
}
 
 //----------------------------------------------------------------------------------------------------------------

lfunction models.codon.MSS.LoadClasses (file) {

    SetDialogPrompt ("A TSV file with three columns (AA, Codon, Class) which is used to partition synonymous substitutions into groups");
    classes = io.ReadDelimitedFile (file, "\t", TRUE);
    headers = utility.Array1D(classes[^'terms.io.header']);
    io.CheckAssertion("`&headers`>=3", "Expected a TSV file with at least 3 columns; 2nd column is the codon, 3rd is the class for this codon");
    codons_by_class = {};
    for (_record_; in; classes [^"terms.io.rows"]) {
        codons_by_class[_record_[1]] = _record_[2];
    }
    
    classes = utility.UniqueValues(codons_by_class);
    class_count = utility.Array1D(classes);
    io.CheckAssertion("`&class_count`>=2", "Expected at least 2 codon classes");

    choices = {class_count,2};
    for (i = 0; i < class_count; i += 1) {
        choices[i][0] = classes[i];
        choices[i][1] = "Codon class " + classes[i];
    }

    nr= io.SelectAnOption  (choices, "Select the codon class which will serve as the neutral rate reference (relative rate = 1)");
    
    return {^"terms.model.MSS.codon_classes" : codons_by_class, ^"terms.model.MSS.neutral" : nr};
}
