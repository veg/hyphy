LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("terms.bf");

function parameters.applyNameSpace (id, namespace) {
	if (Type (namespace) == "String") {
		if (Abs (namespace) > 0) {
			return namespace + "." + id;
		}
	}
	return id;
}

function parameters.declareGlobal (id, cache) {
	if (Abs (id)) {
		if (Type (cache) == "AssociativeList") {
			if (Abs (cache[id]) == 0) {
				return;
			}
		}
		ExecuteCommands ("global `id` = 1;");
	}
}

function parameters.quote (arg) {
	return "\"" + arg + "\"";
}

function parameters.addMultiplicativeTerm (matrix, term) {
	
	if (Abs (term) > 0) {
		__N = Rows (matrix);
	
		for (__r = 0; __r < __N; __r+=1) {
			for (__c = 0; __c < __N; __c+=1) {
				if (__r != __c) {
					if (Abs (matrix[__r][__c])) {
						matrix[__r][__c] += "*" + term;
					} else {
						matrix[__r][__c] =  term;
					}
				}
			}
		}
	}
	
	return matrix;
}

function parameters.stringMatrixToFormulas (id, matrix) {
	__N = Rows (matrix);
	
	ExecuteCommands ("`id` = {__N,__N}");
	
	for (__r = 0; __r < __N; __r+=1) {
		for (__c = 0; __c < __N; __c+=1) {
		
			if (__r != __c && Abs (matrix[__r][__c])) {
				ExecuteCommands ("`id`[__r][__c] := " + matrix[__r][__c]);
			}
		}
	}
	
}

function parameters.setRange (id, ranges) {
	if (Abs (id)) {
		if (Type (ranges) == "AssociativeList") {
			if (Abs (ranges[terms.lower_bound])) {
				ExecuteCommands ("`id` :> " + ranges[terms.lower_bound]);
			} 
			if (Abs (ranges[terms.upper_bound])) {
				ExecuteCommands ("`id` :< " + ranges[terms.upper_bound]);
			} 
		}
	}
}

function parameters.fixValue (id, value) {
	if (Abs (id)) {
		ExecuteCommands ("`id` := " + value);
	}
}