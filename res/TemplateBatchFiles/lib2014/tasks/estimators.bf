LoadFunctionLibrary ("../models/model_functions.bf");
LoadFunctionLibrary ("../models/terms.bf");
LoadFunctionLibrary ("../models/DNA/GTR.bf");

function estimators.copyGlobals2 (key, value) {
    (estimators.extractMLEs.results ["global"])[key] = {"ID" : value, "MLE" : Eval (value)};
}

function estimators.copyGlobals (key, value) {
    ((value["parameters"])["global"])["estimators.copyGlobals2"][""];
}

function estimators.setGlobals2 (key, value) {
    if (Type (initial_values[key]) == "AssociativeList") {
        ExecuteCommands ("`value` = " + (initial_values [key])["MLE"]);
    }
}

function estimators.setGlobals (key, value) {
    ((value["parameters"])["global"])["estimators.setGlobals2"][""];
}

function estimators.extractBranchLength (tree, node, model) {
    if (Abs(model["get_branch_length"])) {
        // this is where we handle non-standard cases
    } 
    return Eval ("BranchLength (`tree`, \"`node`\")");
}

function estimators.applyBranchLength (tree, node, model, length) {
    
}

function estimators.extractMLEs (likelihood_function_id, model_descriptions) {
    ExecuteCommands ("GetString (estimators.extractMLEs.lfInfo, `likelihood_function_id`,-1)");  
    estimators.extractMLEs.results = {};
    estimators.extractMLEs.partitions = Rows (estimators.extractMLEs.lfInfo["Trees"]);
   
    // copy global variables first 
    
    estimators.extractMLEs.results ["global"] = {};
    model_descriptions ["estimators.copyGlobals"][""];
    estimators.extractMLEs.results ["branch lengths"] = {};
    
    for (estimators.extractMLEs.i = 0; estimators.extractMLEs.i < estimators.extractMLEs.partitions; estimators.extractMLEs.i  += 1) {
        _tree_name = (estimators.extractMLEs.lfInfo["Trees"])[estimators.extractMLEs.i];
       
        ExecuteCommands ("GetInformation (estimators.extractMLEs.map, `_tree_name`);");  
        estimators.extractMLEs.branch_names = Rows (estimators.extractMLEs.map);
        (estimators.extractMLEs.results ["branch lengths"])[estimators.extractMLEs.i] = {};
        for (estimators.extractMLEs.b = 0; estimators.extractMLEs.b < Abs(estimators.extractMLEs.map); estimators.extractMLEs.b += 1) {
            _branch_name = estimators.extractMLEs.branch_names[estimators.extractMLEs.b];
            ((estimators.extractMLEs.results ["branch lengths"])[estimators.extractMLEs.i])[_branch_name] = {"MLE" :
                estimators.extractBranchLength (_tree_name, _branch_name, model_descriptions[estimators.extractMLEs.map[_branch_name]])};
        }   
    }   
    
    return estimators.extractMLEs.results;
}

function estimators.applyExistingEstimates (likelihood_function_id, model_descriptions, initial_values) {

    ExecuteCommands ("GetString (estimators.extractMLEs.lfInfo, `likelihood_function_id`,-1)");  
    estimators.extractMLEs.results = {};
    estimators.extractMLEs.partitions = Rows (estimators.extractMLEs.lfInfo["Trees"]);
   
    // copy global variables first 
    
    estimators.extractMLEs.results ["global"] = {};
    model_descriptions ["estimators.setGlobals"][""];
        
    for (estimators.extractMLEs.i = 0; estimators.extractMLEs.i < estimators.extractMLEs.partitions; estimators.extractMLEs.i  += 1) {
        if (Type (initial_values [estimators.extractMLEs.i]) == "AssociativeList") {
            _tree_name = (estimators.extractMLEs.lfInfo["Trees"])[estimators.extractMLEs.i];
       
            ExecuteCommands ("GetInformation (estimators.extractMLEs.map, `_tree_name`);");  
            estimators.extractMLEs.branch_names = Rows (estimators.extractMLEs.map);

            for (estimators.extractMLEs.b = 0; estimators.extractMLEs.b < Abs(estimators.extractMLEs.map); estimators.extractMLEs.b += 1) {
                _branch_name = estimators.extractMLEs.branch_names[estimators.extractMLEs.b];
                _existing_estimate = (initial_values [estimators.extractMLEs.i])[_branch_name];
                if (Type (_existing_estimate) == "AssociativeList") {
                    estimators.applyBranchLength (_tree_name, _branch_name, model_descriptions[estimators.extractMLEs.map[_branch_name]], _existing_estimate["MLE"]);
                }
            }   
            
        }
    }   
}

function estimators.fitGTR  (data_filter, tree, initial_values) {
	// create a nucleotide filter first
	
	DataSetFilter estimators.fitGTR.nuc_data = CreateFilter (^data_filter, 1);
	estimators.fitGTR.model =  model.define_model ("models.DNA.GTR.modelDescription", "estimators.fitGTR.gtr", {"0" : "terms.global"}, "estimators.fitGTR.nuc_data", None); 

    model.applyModelToTree ("estimators.fitGTR.tree", tree, {"default" : estimators.fitGTR.model} , None);
    
    LikelihoodFunction estimators.fitGTR.likelihoodFunction = (estimators.fitGTR.nuc_data, estimators.fitGTR.tree);
    
    if (Type (initial_values) == "AssociativeList") {
        utility.toggleEnvVariable ("USE_LAST_RESULTS", 1);
        estimators.applyExistingEstimates ("estimators.fitGTR.likelihoodFunction", {"estimators.fitGTR.gtr" : estimators.fitGTR.model}, initial_values)
    }
    
    Optimize (estimators.fitGTR.mles, estimators.fitGTR.likelihoodFunction);
    if (Type (initial_values) == "AssociativeList") {
        utility.toggleEnvVariable ("USE_LAST_RESULTS", None);
    }
        
    estimators.fitGTR.results = estimators.extractMLEs ("estimators.fitGTR.likelihoodFunction", {"estimators.fitGTR.gtr" : estimators.fitGTR.model});
    
    estimators.fitGTR.results["LogL"]            = estimators.fitGTR.mles[1][0];
    estimators.fitGTR.results["parameters"] = estimators.fitGTR.mles[1][1] + 3;
    
    //io.spoolLF ("estimators.fitGTR.likelihoodFunction", "/Volumes/home-raid/Desktop/test", None);
    
    DeleteObject (estimators.fitGTR.likelihoodFunction);
    
    return estimators.fitGTR.results;
}

function estimators.fitMGREV  (codon_data, tree, option, initial_values) {
	// create a nucleotide filter first
		
    estimators.fitMGREV.filter = codon_data["dataset"];
		
	DataSetFilter estimators.fitMGREV.codon_data = CreateFilter (^estimators.fitMGREV.filter, 3, "", "", codon_data["stop"]);
	estimators.fitMGREV.model =  model.define_model ("models.codon.MG_REV.modelDescription", "estimators.fitMGREV.mg", {"0" : "terms.global", "1" : codon_data["code"]}, "estimators.fitMGREV.codon_data", None); 

    fprintf (stdout, estimators.fitMGREV.model, "\n");

    return None;

    model.applyModelToTree ("estimators.fitGTR.tree", tree, {"default" : estimators.fitGTR.model} , None);
    
    LikelihoodFunction estimators.fitGTR.likelihoodFunction = (estimators.fitGTR.nuc_data, estimators.fitGTR.tree);
    
    if (Type (initial_values) == "AssociativeList") {
        utility.toggleEnvVariable ("USE_LAST_RESULTS", 1);
        estimators.applyExistingEstimates ("estimators.fitGTR.likelihoodFunction", {"estimators.fitGTR.gtr" : estimators.fitGTR.model}, initial_values)
    }
    
    Optimize (estimators.fitGTR.mles, estimators.fitGTR.likelihoodFunction);
    
    if (Type (initial_values) == "AssociativeList") {
        utility.toggleEnvVariable ("USE_LAST_RESULTS", None);
    }
        
    estimators.fitGTR.results = estimators.extractMLEs ("estimators.fitGTR.likelihoodFunction", {"estimators.fitGTR.gtr" : estimators.fitGTR.model});
    
    estimators.fitGTR.results["LogL"]            = estimators.fitGTR.mles[1][0];
    estimators.fitGTR.results["parameters"] = estimators.fitGTR.mles[1][1] + 3;
    
    //io.spoolLF ("estimators.fitGTR.likelihoodFunction", "/Volumes/home-raid/Desktop/test", None);
    
    DeleteObject (estimators.fitGTR.likelihoodFunction);*/
    
    return estimators.fitMGREV.results;
}
