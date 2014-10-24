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
    __init_value = (initial_values["global"])[key];
    if (Type (__init_value) == "AssociativeList") {
        if (__init_value["fix-me"]) {
            ExecuteCommands ("`value` := " + __init_value["MLE"]);       
        } else {
            ExecuteCommands ("`value` = " + __init_value["MLE"]);
        }
    }
}

function estimators.setGlobals (key, value) {
    ((value["parameters"])["global"])["estimators.setGlobals2"][""];
}

function estimators.extractBranchInformation.copy_local (key, value) {
    estimators.extractBranchLength.result [key] = {"ID": value, "MLE": Eval (estimators.extractBranchLength.parameter_tag + "." + value)};
}

function estimators.extractBranchInformation (tree, node, model) {
    estimators.extractBranchLength.result = {}; 
    
    if (Abs(model["get-branch-length"])) {
        estimators.extractBranchLength.result ["MLE"] = utility.callFunction (model["get-branch-length"], {"0" : "model", "1" : "tree",  "2": "node"});
    } else {
        estimators.extractBranchLength.result ["MLE"] = Eval ("BranchLength (`tree`, \"`node`\")");
    }
    
    estimators.extractBranchLength.parameter_tag = tree + "." + node;
    (model.parameters.local (model))["estimators.extractBranchInformation.copy_local"][""];
    
    return estimators.extractBranchLength.result;
}

function estimators.applyBranchLength (tree, node, model, length) {
    utility.callFunction (model["set-branch-length"], {"0" : "model", "1" : length, "2": parameters.quote(tree + "." + node)});
}

function estimators.fixSubsetOfEstimates.helper (key, value) {
    value["fix-me"] = 1;
}

function estimators.fixSubsetOfEstimates.helper_condition (key) {
    return Type (variables[key]) != Unknown;
}

function estimators.fixSubsetOfEstimates (estimates, variables) {
    (estimates["global"])["estimators.fixSubsetOfEstimates.helper"]["estimators.fixSubsetOfEstimates.helper_condition"];
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
            ((estimators.extractMLEs.results ["branch lengths"])[estimators.extractMLEs.i])[_branch_name] = 
                estimators.extractBranchInformation (_tree_name, _branch_name, model_descriptions[estimators.extractMLEs.map[_branch_name]]);
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
        if (Type ((initial_values["branch lengths"])[estimators.extractMLEs.i]) == "AssociativeList") {
        
            _tree_name = (estimators.extractMLEs.lfInfo["Trees"])[estimators.extractMLEs.i];
       
            ExecuteCommands ("GetInformation (estimators.extractMLEs.map, `_tree_name`);");  
            estimators.extractMLEs.branch_names = Rows (estimators.extractMLEs.map);

            for (estimators.extractMLEs.b = 0; estimators.extractMLEs.b < Abs(estimators.extractMLEs.map); estimators.extractMLEs.b += 1) {
                _branch_name = estimators.extractMLEs.branch_names[estimators.extractMLEs.b];
                _existing_estimate = ((initial_values["branch lengths"]) [estimators.extractMLEs.i])[_branch_name];
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
	estimators.fitGTR.model =  model.generic.define_model ("models.DNA.GTR.modelDescription", "estimators.fitGTR.gtr", {"0" : "terms.global"}, "estimators.fitGTR.nuc_data", None); 

    model.applyModelToTree ("estimators.fitGTR.tree", tree, {"default" : estimators.fitGTR.model} , None);
    
    LikelihoodFunction estimators.fitGTR.likelihoodFunction = (estimators.fitGTR.nuc_data, estimators.fitGTR.tree);
    
    if (Type (initial_values) == "AssociativeList") {
        utility.toggleEnvVariable ("USE_LAST_RESULTS", 1);
        estimators.applyExistingEstimates ("estimators.fitGTR.likelihoodFunction", {"estimators.fitGTR.gtr" : estimators.fitGTR.model}, initial_values);
    }
    
    Optimize (estimators.fitGTR.mles, estimators.fitGTR.likelihoodFunction);
    if (Type (initial_values) == "AssociativeList") {
        utility.toggleEnvVariable ("USE_LAST_RESULTS", None);
    }
        
    estimators.fitGTR.results = estimators.extractMLEs ("estimators.fitGTR.likelihoodFunction", {"estimators.fitGTR.gtr" : estimators.fitGTR.model});
    
    estimators.fitGTR.results["LogL"]            = estimators.fitGTR.mles[1][0];
    estimators.fitGTR.results["parameters"] = estimators.fitGTR.mles[1][1] + 3;
    
    DeleteObject (estimators.fitGTR.likelihoodFunction);
    
    return estimators.fitGTR.results;
}

function estimators.fitMGREV.set_partition_omega (key, value) {
    Eval  ("estimators.fitMGREV.tree.`key`.`estimators.fitMGREV.alpha` = 0.1");
    ExecuteCommands ("estimators.fitMGREV.tree.`key`.`estimators.fitMGREV.beta`:=estimators.fitMGREV.tree.`key`.`estimators.fitMGREV.alpha`*" + estimators.fitMGREV.partitioned_omega.parameters [value]);
}

function estimators.fitMGREV  (codon_data, tree, option, initial_values) {
		
    estimators.fitMGREV.filter = codon_data["dataset"];
    estimators.fitMGREV.partitioned_omega = 0;
		
	DataSetFilter estimators.fitMGREV.codon_data = CreateFilter (^estimators.fitMGREV.filter, 3, "", "", codon_data["stop"]);
	
	estimators.fitMGREV.model =  model.generic.define_model ("models.codon.MG_REV.modelDescription", 
	                                                         "estimators.fitMGREV.mg", 
	                                                         {"0" : parameters.quote (option["model-type"]), "1" : codon_data["code"]}, 
	                                                         "estimators.fitMGREV.codon_data", 
	                                                         None);
	                                                          
    model.applyModelToTree ("estimators.fitMGREV.tree", tree, {"default" : estimators.fitMGREV.model} , None);

    if (option["model-type"] == terms.local && option["partitioned-omega"]) {
        estimators.fitMGREV.model_list = tree["model_list"];
        if (Columns(estimators.fitMGREV.model_list) > 1) {
            estimators.fitMGREV.partitioned_omega = 1;  
        }
    }

    LikelihoodFunction estimators.fitMGREV.likelihoodFunction = (estimators.fitMGREV.codon_data, estimators.fitMGREV.tree);
    
    if (estimators.fitMGREV.partitioned_omega) {
        estimators.fitMGREV.partitioned_omega.parameters = (estimators.fitMGREV.model["parameters"])["global"];
        
        for (estimators.fitMGREV.i = 0; estimators.fitMGREV.i < Columns(estimators.fitMGREV.model_list); estimators.fitMGREV.i += 1) {
            estimators.fitMGREV.partitioned_omega.parameters [estimators.fitMGREV.model_list[estimators.fitMGREV.i]] = "estimators.fitMGREV.mg.omega_" + estimators.fitMGREV.model_list[estimators.fitMGREV.i];
        }
        
        parameters.declareGlobal (estimators.fitMGREV.partitioned_omega.parameters, None);

        estimators.fitMGREV.lp    = model.parameters.local (estimators.fitMGREV.model);
        estimators.fitMGREV.beta  = estimators.fitMGREV.lp[terms.nonsynonymous_rate];
        estimators.fitMGREV.alpha = estimators.fitMGREV.lp[terms.synonymous_rate];
     
        (tree["model_map"])["estimators.fitMGREV.set_partition_omega"][""];
        
        (estimators.fitMGREV.model["parameters"])["global"] = estimators.fitMGREV.partitioned_omega.parameters;
        
     }

    if (Type (initial_values) == "AssociativeList") {
        utility.toggleEnvVariable ("USE_LAST_RESULTS", 1);
        estimators.applyExistingEstimates ("estimators.fitMGREV.likelihoodFunction", {"estimators.fitMGREV.mg" : estimators.fitMGREV.model}, initial_values);
    }
    
    // io.spoolLF ("estimators.fitMGREV.likelihoodFunction", "/Volumes/home-raid/Desktop/test", None);
    Optimize (estimators.fitMGREV.mles, estimators.fitMGREV.likelihoodFunction);
    
    if (Type (initial_values) == "AssociativeList") {
        utility.toggleEnvVariable ("USE_LAST_RESULTS", None);
    }
    
    estimators.fitMGREV.results = estimators.extractMLEs ("estimators.fitMGREV.likelihoodFunction", {"estimators.fitMGREV.mg" : estimators.fitMGREV.model});
    estimators.fitMGREV.results["LogL"]            = estimators.fitMGREV.mles[1][0];
    estimators.fitMGREV.results["parameters"]      = estimators.fitMGREV.mles[1][1] + 9;
    
    //io.spoolLF ("estimators.fitMGREV.likelihoodFunction", "/Volumes/home-raid/Desktop/test", None);
    
    DeleteObject (estimators.fitMGREV.likelihoodFunction);
    
    return estimators.fitMGREV.results;
}
