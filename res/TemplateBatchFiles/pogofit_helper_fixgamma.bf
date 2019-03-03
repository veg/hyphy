
/********************************************************************************************************************/
/********************************************* MODEL DEFINITIONS ****************************************************/
/********************************************************************************************************************/

/**
 * @name models.protein.Baseline.ModelDescription.withGamma
 * @description Define baseline (standard matrix) model w/ +F and *no* four-category gamma rate variation
 */
function pogofit.Baseline.ModelDescription(type){
    def = Call( models.protein.empirical.plusF_generators[pogofit.baseline_model], type);
    return def;
}
/**
 * @name models.protein.Baseline.ModelDescription.withGamma
 * @description Define baseline (standard matrix) model w/ +F and *yes* four-category gamma rate variation
 */
function pogofit.Baseline.ModelDescription.withGamma(type){
    def = pogofit.Baseline.ModelDescription(type);
	def [utility.getGlobalValue("terms.model.rate_variation")] = rate_variation.types.Gamma.factory ({utility.getGlobalValue("terms.rate_variation.bins") : 4});
    return def;
}


/**
 * @name pogofit.REV.ModelDescription.freqs
 * @description Define REV model frequencies as empirical
 */
function pogofit.REV.ModelDescription.freqs (model, namespace, datafilter) {
    model[terms.efv_estimate] = pogofit.shared_EFV;
    model[terms.model.efv_estimate_name] = terms.frequencies.predefined;
    (model[terms.parameters])[terms.model.empirical] = 0;
    return model;
}

/**
 * @name pogofit.REV.ModelDescription
 * @description Define a REV model with constant site rates
 */
function pogofit.REV.ModelDescription (type) {
    def = models.protein.REV.ModelDescription(type);
    if (Type (pogofit.shared_EFV) == "Matrix") {
        def [terms.model.frequency_estimator] = "pogofit.REV.ModelDescription.freqs";
    }
    return def;
}

/**
 * @name pogofit.REV.ModelDescription.withGamma
 * @description Define a REV model with Gamma rate variation
 */
function pogofit.REV.ModelDescription.withGamma (type) {
    def = pogofit.REV.ModelDescription(type);
    def [utility.getGlobalValue("terms.model.rate_variation")] = rate_variation.types.Gamma.factory ({utility.getGlobalValue("terms.rate_variation.bins"): 4});
    return def;
}
/********************************************************************************************************************/
/********************************************************************************************************************/
/********************************************************************************************************************/



/********************************************************************************************************************/
/********************************************** FITTING FUNCTIONS ***************************************************/
/********************************************************************************************************************/

function pogofit.fitGTR_fixalpha (current_results) {

    filter_info    = {};
    trees = {};
    initial_values = {terms.global : {}, terms.branch_length : {}};
    index_to_file_name   = {};
    
    
    // SLKP: create constrain branch length constraint list
    //bl_constraints = {};
    
    for (file_index = 0; file_index < pogofit.file_list_count; file_index += 1) {
        file_path = pogofit.file_list [file_index];
        dataset_name = "pogofit.msa.file_" + file_index;
        data_info = alignments.ReadNucleotideDataSet (dataset_name, file_path);
        data_info = alignments.EnsureMapping (dataset_name, data_info);

        partition_specification = { "0" : {terms.data.name : "all", terms.data.filter_string : "", term.data.tree : ((current_results[file_index])[terms.fit.trees])[0]}};


        this_bl = (current_results[terms.branch_length])[file_index];
        this_tree = (current_results[terms.fit.trees])[file_index];

        filter_info [file_index] = (alignments.DefineFiltersForPartitions (partition_specification,
                                                                            dataset_name ,
                                                                            dataset_name,
                                                                            data_info))[0];
        trees [file_index] = {terms.trees.newick :  this_tree};
        (initial_values[terms.branch_length])[file_index] = this_bl;
       // bl_constraints [file_index] = terms.model.branch_length_constrain;
        
        
    }
    filter_names = utility.Map (filter_info, "_value_", "_value_[terms.data.name]");

    utility.ToggleEnvVariable ("AUTO_PARALLELIZE_OPTIMIZE", 1);
    utility.ToggleEnvVariable ("OPTIMIZATION_METHOD", 0);
    utility.ToggleEnvVariable ("VERBOSITY_LEVEL", 1);
    
    // set intial values to user chosen matrix (same as baseline)
    for (l1 = 0; l1 < 20; l1 += 1) {
        for (l2 = l1 + 1; l2 < 20; l2 += 1) {
            rate_term = terms.aminoacidRate (models.protein.alphabet[l1],models.protein.alphabet[l2]);
            (initial_values[terms.global]) [rate_term] = {terms.fit.MLE : (pogofit.initial_rates[models.protein.alphabet[l1]])[models.protein.alphabet[l2]]}; 
        }
    }
    alpha_term = "Gamma distribution shape parameter";
    alpha = ((current_results[terms.global])[alpha_term])[terms.fit.MLE];
    (initial_values[terms.global]) [alpha_term] = {terms.fit.MLE : alpha , terms.fix : TRUE};

    
    pogofit.rev_mle = estimators.FitSingleModel_Ext (
                                        filter_names,
                                        trees,
                                        pogofit.rev_model,
                                        initial_values,
                                        {
                                          // terms.run_options.model_type: terms.global,
                                           //terms.run_options.proportional_branch_length_scaler: bl_constraints,
                                            terms.run_options.retain_lf_object : TRUE
                                        }
                                   );                         
    /*   
    // Uncomment these lines if you'd like to save the NEXUS LF.                        
    lf_id = pogofit.rev.mle[terms.likelihood_function];
    Export(pogofit.finalphase_LF, ^lf_id);
    fprintf(pogofit.final_likelihood_function, pogofit.finalphase_LF);
    */
    pogofit.rev_mle - terms.likelihood_function;
    
    // Save the rev.mle into the analysis_results, and cache it.
    (^"pogofit.analysis_results")[pogofit.final_phase] = pogofit.rev_mle;

    console.log (""); // clear past the optimization progress line
    utility.SetEnvVariable ("VERBOSITY_LEVEL", 0);
    utility.ToggleEnvVariable ("AUTO_PARALLELIZE_OPTIMIZE", None);
    utility.ToggleEnvVariable ("OPTIMIZATION_METHOD", None);


    // Trees as dictionary for compatibility with rest of the output.
    pogofit.rev_mle[terms.fit.trees] = utility.SwapKeysAndValues(utility.MatrixToDict(pogofit.rev_mle[terms.fit.trees]));
    
    return pogofit.rev_mle;

}



function pogofit.fitBaselineTogether () {


    filter_info    = {};
    trees = {};
    index_to_file_name   = {};
    
    for (file_index = 0; file_index < pogofit.file_list_count; file_index += 1) {
        file_path = pogofit.file_list [file_index];
        dataset_name = "pogofit.msa.file_" + file_index;
        
        file_info = alignments.ReadNucleotideDataSet (dataset_name, file_path);
        file_info = alignments.EnsureMapping(dataset_name, file_info);

        utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", "/dev/null");
        ExecuteCommands ('partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (file_info[terms.data.partitions], name_mapping)',
                         {"0" : "Y"});
        utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", None);

        tree = utility.Map (partitions_and_trees, "_value_", '_value_[terms.data.tree]');

        tree_with_lengths = (tree["0"])[terms.trees.newick_with_lengths];
        
        /* SAVE TO JSON */
        (pogofit.input[pogofit.options.dataset_information])[file_index] = { 
                                        terms.json.sequences: file_info[terms.data.sequences],
                                        terms.json.sites: file_info[terms.data.sites],
                                        terms.json.trees: tree_with_lengths,
                                        terms.json.file: file_path};
  

        partition_specification = { "0" : {
                                            terms.data.name : "all", 
                                            terms.data.filter_string : "", 
                                            terms.data.tree : tree_with_lengths}};
             
        filter_info [file_index] = (alignments.DefineFiltersForPartitions (partition_specification,
                                                                            dataset_name ,
                                                                            dataset_name,
                                                                            file_info))[0];                                                        
        trees [file_index] = {terms.trees.newick : tree_with_lengths};        
    }

    utility.SetEnvVariable ("VERBOSITY_LEVEL", 1);
    utility.ToggleEnvVariable ("AUTO_PARALLELIZE_OPTIMIZE", 1);
    utility.ToggleEnvVariable ("OPTIMIZATION_METHOD", 0);

    pogofit.baseline_mle = estimators.FitSingleModel_Ext (
                                            utility.Map (filter_info, "_value_", "_value_[terms.data.name]"),
                                            trees,
                                            pogofit.baseline_model_desc,
                                            None,
                                            None
                                           );
    
    // Save the rev.mle into the analysis_results
    (^"pogofit.analysis_results")[pogofit.baseline_phase] = pogofit.baseline_mle;

    console.log (""); // clear past the optimization progress line
    utility.SetEnvVariable ("VERBOSITY_LEVEL", 0);
    utility.ToggleEnvVariable ("AUTO_PARALLELIZE_OPTIMIZE", None);
    utility.ToggleEnvVariable ("OPTIMIZATION_METHOD", None);


    // Trees as dictionary for compatibility with rest of the output.
    pogofit.baseline_mle[terms.fit.trees] = utility.SwapKeysAndValues(utility.MatrixToDict(pogofit.baseline_mle[terms.fit.trees]));
    
    return pogofit.baseline_mle;

}




/********************************************************************************************************************/
/********************************************************************************************************************/
/********************************************************************************************************************/



/********************************************************************************************************************/
/********************************************** UTILITY FUNCTIONS ***************************************************/
/********************************************************************************************************************/

lfunction pogofit.startTimer(timers, key) {
    timers[key] = {
        terms.timers.timer: Time(1),
    };

}
lfunction pogofit.stopTimer(timers, key) {
    (timers[key])[terms.timers.timer] = Time(1) - (timers[key])[terms.timers.timer];
}




// From the fitted model results and create rate dictionary which can be used for .fitted_model
function pogofit.extract_rates() {

    rij = {};
    for (l1 = 0; l1 < 20 - 1; l1 += 1) {
        rij[models.protein.alphabet[l1]] = {};
        for (l2 = l1 + 1; l2 < 20; l2 += 1) {
    
            rate_search = terms.aminoacidRate(models.protein.alphabet[l1], models.protein.alphabet[l2]);       
            (rij[models.protein.alphabet[l1]])[models.protein.alphabet[l2]] = ((pogofit.gtr_fit[terms.global])[rate_search])[terms.fit.MLE];
        }
    }
    return rij
}


// From the fitted model results and create rate dictionary which can be used for .fitted_model
function pogofit.extract_rates_imputation() {


    pogofit.thekeys =  utility.Keys(pogofit.site_counts);
    summation = 0.0;
    for (x = 0; x < utility.Array1D(pogofit.thekeys); x+=1)
    {
        this_key = pogofit.thekeys[x];
        summation += ((1. + pogofit.site_counts[this_key]) * pogofit.tree_lengths[this_key] );            
    }
    
    rij = {};
    for (l1 = 0; l1 < 20 - 1; l1 += 1) {
        rij[models.protein.alphabet[l1]] = {};
        for (l2 = l1 + 1; l2 < 20; l2 += 1) {
            rate_search = terms.aminoacidRate(models.protein.alphabet[l1], models.protein.alphabet[l2]);
            this_rate = ((pogofit.gtr_fit[terms.global])[rate_search])[terms.fit.MLE];
            if (this_rate == 0.0) 
            {
                efv_sum = pogofit.final_efv[l1] + pogofit.final_efv[l2];
                this_rate = 1./ ((pogofit.final_efv[l1] + pogofit.final_efv[l2]) * summation);
            }
            (rij[models.protein.alphabet[l1]])[models.protein.alphabet[l2]] = this_rate;
        }
    }
    return rij
}




function pogofit.write_model_to_file() {

    pogofit.external_order   = {{"A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"}};
    pogofit.hyphy_order_dict = utility.MatrixToDict(models.protein.alphabet); 

    if (pogofit.output_format == pogofit.output_hyphy)
    {
        pogofit.save_hyphy_model();
    }

    if (pogofit.output_format == pogofit.output_paml)
    {
        pogofit.save_paml_model();
    }
    if (pogofit.output_format == pogofit.output_raxml)
    {
        pogofit.save_raxml_model();
    }
    if (pogofit.output_format == pogofit.output_all)
    {
        pogofit.save_raxml_model();
        pogofit.save_paml_model();
        pogofit.save_hyphy_model();
    }
}


function pogofit.save_hyphy_model(){
    fprintf(pogofit.output_model_prefix + pogofit.hyphy_model_ext, CLEAR_FILE, "Rij = " + pogofit.final_rij + ";");
    fprintf(pogofit.output_model_prefix + pogofit.hyphy_model_ext, "\n\n\n");
    fprintf(pogofit.output_model_prefix + pogofit.hyphy_model_ext, "EFV = " + pogofit.final_efv + ";");
}

function pogofit.save_raxml_model(){

    pogofit.raxml_output = "";
    pogofit.raxml_file = pogofit.output_model_prefix + pogofit.raxml_model_ext;
    
    for (i = 0; i < 20; i +=1)
    {
        for (j = 0; j < 20; j += 1)
        {
        
            if (i == j)
            {
                pogofit.raxml_output += "0.0\n";
            }
            else
            {
                pogofit.aa1 = pogofit.external_order[i];
                pogofit.aa2 = pogofit.external_order[j];
        
                if (pogofit.aa1 < pogofit.aa2)
                {
                    pogofit.rate = (pogofit.final_rij[pogofit.aa1])[pogofit.aa2];
                }
                else
                {
                    pogofit.rate = (pogofit.final_rij[pogofit.aa2])[pogofit.aa1];
                }
                pogofit.raxml_output += pogofit.rate;
                pogofit.raxml_output += "\n"; 
            }
        }    
    }
    for (i = 0; i < 20; i += 1)
    {
        pogofit.raxml_output += pogofit.final_efv[ pogofit.hyphy_order_dict[pogofit.external_order[i]] ];
        // strip hack
        if (i <= 18){
            pogofit.raxml_output += "\n";
        }
    }

    fprintf(pogofit.raxml_file, CLEAR_FILE, pogofit.raxml_output);
}




function pogofit.save_paml_model(){

    pogofit.paml_output = "";
    pogofit.paml_file = pogofit.output_model_prefix + pogofit.paml_model_ext;

    for (i = 1; i < 20; i +=1)
    {
        pogofit.row = "";
        for (j = 0; j < i; j += 1)
        {
        
            pogofit.aa1 = pogofit.external_order[i];
            pogofit.aa2 = pogofit.external_order[j];
        
            if (pogofit.aa1 < pogofit.aa2){
                pogofit.rate = (pogofit.final_rij[pogofit.aa1])[pogofit.aa2];
            }
            else{
                pogofit.rate = (pogofit.final_rij[pogofit.aa2])[pogofit.aa1];
            }
            
            pogofit.row += pogofit.rate;
            // strip hack
            if (j != (i-1)){
                pogofit.row += " ";
            }
        }
    
        pogofit.paml_output += pogofit.row + "\n";
    }

    pogofit.paml_output += "\n";
    for (i = 0; i < 20; i += 1)
    {
        pogofit.paml_output += pogofit.final_efv[ pogofit.hyphy_order_dict[pogofit.external_order[i]] ];
        // strip hack
        if (i <= 18){
            pogofit.paml_output += " ";
        }
    }
    fprintf(pogofit.paml_file, CLEAR_FILE, pogofit.paml_output);
}
/********************************************************************************************************************/
/********************************************************************************************************************/
/********************************************************************************************************************/


