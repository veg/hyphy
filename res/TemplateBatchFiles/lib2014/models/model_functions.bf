LoadFunctionLibrary ("DNA.bf");
LoadFunctionLibrary ("parameters.bf");
LoadFunctionLibrary ("frequencies.bf");
LoadFunctionLibrary ("../UtilityFunctions.bf");
//------------------------------------------------------------------------------ 


function model.applyModelToTree (id, tree, model_list, rules) {


	if (Type (rules) == "AssociativeList") {
	    // this has the form 
	    // model id : list of branches to apply the model (as a string matrix)
	    // OR 
	    // DEFAULT : model id
	    
	    if (Abs (rules["DEFAULT"])) {
            ExecuteCommands ("UseModel (" + rules["DEFAULT"] + ");
                              Tree `id` = " + tree["string"] + ";
                              ");	    
	    }
	    
	    model.applyModelToTree.ids = Rows (rules);
	    for (model.applyModelToTree.k = 0; model.applyModelToTree.k < Abs (rules); model.applyModelToTree.k += 1) {
	        model.applyModelToTree.name = model.applyModelToTree.ids[model.applyModelToTree.k];
	        if ( model.applyModelToTree.name != "DEFAULT") {
                model.applyModelToTree.list = rules[model.applyModelToTree.name];
                for (model.applyModelToTree.b = 0; model.applyModelToTree.b < Columns (model.applyModelToTree.list); model.applyModelToTree.b += 1) {
                    ExecuteCommands ("SetParameter (`id`." + model.applyModelToTree.list[model.applyModelToTree.b] + ",MODEL," + model.applyModelToTree.name + ")");
                }
            }
	    }
	    
	} else {
		model.applyModelToTree.modelID = model_list[model_list ["INDEXORDER"][0]];
		ExecuteCommands ("UseModel (" + model.applyModelToTree.modelID["id"] + ");
						  Tree `id` = " + tree["string"] + ";
						  ");
	}
}
//------------------------------------------------------------------------------ 

function model.define_model (model_spec, id, type, data_filter, estimator_type) {
	
	model.define_model.model = utility.callFunction (model_spec, None);
	models.generic.attachFilter (model.define_model.model, data_filter);
	
	model.define_model.model = utility.callFunction(model.define_model.model ["defineQ"], {"0" :   parameters.quote (type),
																						   "1" :    parameters.quote (id)});
	model.define_model.model ["matrix id"] = "`id`_" + terms.rate_matrix;
	model.define_model.model ["efv id"] = "`id`_" + terms.efv_matrix;
	model.define_model.model ["id"] = id;
	
	
	if (estimator_type != None) {
		model.define_model.model ["frequency_estimator"] = estimator_type;
	} 
		
	utility.callFunction (model.define_model.model ["frequency_estimator"], {"0": "model.define_model.model", 
													    "1":  parameters.quote(id),
													    "2":   parameters.quote(data_filter)}); // this sets the EFV field
													  
													  
	parameters.stringMatrixToFormulas (model.define_model.model ["matrix id"],model.define_model.model[terms.rate_matrix]);
	ExecuteCommands (model.define_model.model ["efv id"]   + " = " + model.define_model.model[terms.efv_estimate]);
		
	ExecuteCommands ("Model `id` = (" + model.define_model.model ["matrix id"] + "," + model.define_model.model ["efv id"] + "," + model.define_model.model ["canonical"] + ")");
	
	
	return model.define_model.model;
}

//------------------------------------------------------------------------------ 

function models.generic.attachFilter (model, filter) {
	GetDataInfo (_givenAlphabet, *filter, "CHARACTERS");
	__alphabet = model ["alphabet"];
	assert (Columns (__alphabet) == Columns (_givenAlphabet) && model.matchAlphabets (_givenAlphabet, __alphabet), "The declared model alphabet '" + __alphabet + "' does not match the `filter` filter: '" + _givenAlphabet + "'");
	
	model ["alphabet"] = _givenAlphabet;
	model ["data"] = filter;
	return model;
}


//------------------------------------------------------------------------------ 

lfunction model.matchAlphabets (a1, a2) {
	_validStates = {};
	for (_k = 0; _k < Columns (a1); _k += 1) {
		_validStates [a1[_k]] = _k;
	}
	for (_k = 0; _k < Columns (a2); _k += 1) {
		if (_validStates [a2[_k]] != _k) {
			return 0;
		}
	}
	return 1;

}