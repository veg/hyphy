RequireVersion ("2.5.71");

LoadFunctionLibrary("MG_REV_PROPERTIES.bf");
LoadFunctionLibrary("BS_REL.bf");

//----------------------------------------------------------------------------------------------------------------


lfunction model.codon.MG_REV_PROPERTIES_BSREL.prompt_and_define (type, code) {
    KeywordArgument ("property-set", "How to partition synonymous codons into classes", "Atchley");
    
    property_set = io.SelectAnOption (
            {
                "Atchley":"Use the five properties derived from a factor analysis of 500 amino-acid properties [Table 2 in PNAS (2005) 102(18) 6395-6400 doi: 10.1073/pnas.0408677102]",
                "LCAP":"Use the five properties defined in the Conant and Stadler LCAP model [Mol Biol Evol (2009) 26 (5): 1155-1161. doi: 10.1093/molbev/msp031]",
                "3PROP" : "Use three primary properties: Hydrophobicity (Kyle-Doolitle), Volume (cubic angstroms), Isoelectric Point (pI)",
                "Random-2" : "Two random properties (for null hypothesis testing)",
                "Random-3" : "Three random properties (for null hypothesis testing)",
                "Random-4" : "Four random properties (for null hypothesis testing)",
                "Random-5" : "Five random properties (for null hypothesis testing)",
                "Custom":"Load the set of properties from a file"
            }, 
            "The set of properties to use in the model");

    
    
    if (property_set == "Custom") {
        KeywordArgument ("property-file", "JSON file which defines amino-acid properties");
        property_set = io.PromptUserForFilePathRead ("JSON file which defines amino-acid properties");
        property_set = io.ParseJSON(property_set);
        console.log (">Loaded a set of `Abs(property_set)` properties");
     }
    
     KeywordArgument ("components", "How many model components?", "2");    
     components = io.PromptUser ("How many model components?", 2, 1, 10, TRUE);
     return models.codon.MG_REV_PROPERTIES_BSREL.ModelDescription(type, code, property_set, components);
}


lfunction models.codon.MG_REV_PROPERTIES_BSREL.ModelDescription(type, code, properties, components) {

    io.CheckAssertion ("`&components` >= 1 && `&components` <= 10", "must have between 1 and 10 components in call to models.codon.MG_REV_PROPERTIES_BSREL.ModelDescription");
    mg_base = models.codon.MG_REV_PROPERTIES.ModelDescription (type, code, properties);
    mg_base[  utility.getGlobalValue("terms.model.reversible") ] = TRUE;
    mg_base[  utility.getGlobalValue("terms.model.canonical") ] = "EXPLICIT_FORM_MATRIX_EXPONENTIAL";
    mg_base[  utility.getGlobalValue("terms.model.q_ij") ] = "";
    mg_base[  utility.getGlobalValue("terms.model.components") ] = components;
    mg_base[  utility.getGlobalValue("terms.model.defineQ") ] = "models.codon.MG_REV_PROPERTIES_BSREL._DefineQ";
    mg_base [utility.getGlobalValue("terms.description")] = "The Muse-Gaut 94 codon-substitution model coupled with the general time reversible (GTR) model of nucleotide substitution, which  incorporates amino-acid residue properties into the non-synonymous rates and allows branch-site mixture support";
    mg_base [ utility.getGlobalValue("terms.model.post_definition") ] = "models.codon.MG_REV_PROPERTIES_BSREL.post_definition";
    return mg_base;
}

lfunction  models.codon.MG_REV_PROPERTIES_BSREL.post_definition (model) {
    prop_range = {
        ^"terms.lower_bound": "-10",
        ^"terms.upper_bound": "10"
    };
    
    for (id; in ; model.GetParameters_RegExp (model, terms.propertyImportance ('', '') + "|" + ^"terms.parameters.log_omega_ratio")) {
        parameters.SetRange(id, prop_range);
        parameters.SetValue(id, 0.1);
    }

    return models.codon.BS_REL.post_definition (model);
}

lfunction models.codon.MG_REV_PROPERTIES_BSREL._DefineQ(bs_rel, namespace) {

    rate_matrices = {};

    bs_rel [utility.getGlobalValue("terms.model.q_ij")] = &rate_generator;
    bs_rel [utility.getGlobalValue("terms.mixture.mixture_components")] = {};

    _aux = parameters.GenerateSequentialNames (namespace + ".bsrel_mixture_aux", bs_rel[utility.getGlobalValue("terms.model.components")] - 1, "_");
    _wts = parameters.helper.stick_breaking (_aux, None);
    mixture = {};
    
    

    for (component = 1; component <= bs_rel[utility.getGlobalValue("terms.model.components")]; component += 1) {
       key = "component_" + component;
       ExecuteCommands ("
        function rate_generator (fromChar, toChar, namespace, model_type, model) {
               return models.codon.MG_REV_PROPERTIES._GenerateRate_generic (fromChar, toChar, namespace, model_type, model[utility.getGlobalValue('terms.translation_table')],
                ^'models.codon.BS_REL.rate_term', utility.getGlobalValue('terms.parameters.synonymous_rate'),
                'beta_`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.nonsynonymous_rate'), component),
                'lambda`component`', '' + component,
                'omega`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.log_omega_ratio'), component),
                model[utility.getGlobalValue('terms.model.residue_properties')],
                model[utility.getGlobalValue('terms.model.residue_name_map')])
            }"
       );

       if ( component < bs_rel[utility.getGlobalValue("terms.model.components")]) {
            model.generic.AddGlobal ( bs_rel, _aux[component-1], terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), component ));
            parameters.DeclareGlobalWithRanges (_aux[component-1], 0.5, 0, 1);
       }
        models.codon.generic.DefineQMatrix(bs_rel, namespace);
       rate_matrices [key] = bs_rel[utility.getGlobalValue("terms.model.rate_matrix")];
       (bs_rel [^'terms.mixture.mixture_components'])[key] = _wts [component-1];
    }


    bs_rel[utility.getGlobalValue("terms.model.rate_matrix")] = rate_matrices;
    parameters.SetConstraint(((bs_rel[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.global")])[terms.nucleotideRate("A", "G")], "1", "");


    return bs_rel;
}



