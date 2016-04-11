LoadFunctionLibrary ("parameters.bf");

rate_variation.types = {
    "Gamma+I": "rate_variation.types.gamma_i"
};

function rate_variation.add (model, specs) {

    rate_variation.add.spec = rate_variation.types [specs["type"]];
    assert (Type (rate_variation.add.spec) == "String", 
            "Unsupported rate variation type '`type`' in call to rate_variation. Use one of " + Rows (rate_variation.types));
            
    model  [terms.rate_variation.bins] = specs["bins"];  
    model ["original q call-back"]    = model["defineQ"];     
    model ["defineQ"] =  rate_variation.add.spec;
    return model;
}

function rate_variation.multiply_in (matrix, variable) {
    
}



function  rate_variation.types.gamma_i (model, namespace) {
    
    __global_cache = {};
	
	
	
	rate_variation.types.gamma_i.alpha = parameters.applyNameSpace ("alpha", namespace);
	parameters.declareGlobal (rate_variation.types.gamma_i.alpha, __global_cache);
    
    rate_variation.types.gamma_i.range = {};
    rate_variation.types.gamma_i.range [terms.lower_bound] = 0.001;
    rate_variation.types.gamma_i.range [terms.upper_bound] = 100;
    
	parameters.setRange (rate_variation.types.gamma_i.alpha, 
                         rate_variation.types.gamma_i.range);


	rate_variation.types.gamma_i.inv_p = parameters.applyNameSpace ("inv_p", namespace);
	parameters.declareGlobal (rate_variation.types.gamma_i.inv_p, __global_cache);
	
	parameters.set_value (rate_variation.types.gamma_i.inv_p, 0.1);
	
	parameters.setRange (rate_variation.types.gamma_i.inv_p, terms.range01);

	rate_variation.types.gamma_i.beta = parameters.applyNameSpace ("beta", namespace);
    parameters.setConstraint (rate_variation.types.gamma_i.beta, 
                                "(1-`rate_variation.types.gamma_i.inv_p`)*`rate_variation.types.gamma_i.alpha`", 
                                "global");

    rate_variation.types.gamma_i.range [terms.upper_bound] = 200;
	parameters.setRange (rate_variation.types.gamma_i.beta, 
                         rate_variation.types.gamma_i.range);
	
	((model["parameters"])["global"])[terms.rate_variation.gamma_alpha] = rate_variation.types.gamma_i.alpha;
	((model["parameters"])["global"])[terms.rate_variation.gamma_beta] = rate_variation.types.gamma_i.beta;
	((model["parameters"])["global"])[terms.rate_variation.gamma_p_inv] = rate_variation.types.gamma_i.inv_p;
	
	rate_variation.types.gamma_i.bins = 0 + model [terms.rate_variation.bins];
	
	
	rate_variation.types.gamma_i.main = parameters.applyNameSpace ("gamma_i", namespace);
	rate_variation.types.gamma_i.cats = parameters.applyNameSpace ("gamma_i.freqs", namespace);
	
	ExecuteCommands (rate_variation.types.gamma_i.cats + " = {1 + rate_variation.types.gamma_i.bins, 1}");
    ExecuteCommands (rate_variation.types.gamma_i.cats + "[0] := `rate_variation.types.gamma_i.inv_p`");
    	
	for (rate_variation.types.gamma_i.i=1; rate_variation.types.gamma_i.i<=rate_variation.types.gamma_i.bins;rate_variation.types.gamma_i.i += 1){
        ExecuteCommands (rate_variation.types.gamma_i.cats + 
            "[rate_variation.types.gamma_i.i] := (1-`rate_variation.types.gamma_i.inv_p`)/" + 
            rate_variation.types.gamma_i.bins);
            
    }
    	
	
	ExecuteCommands ("category `rate_variation.types.gamma_i.main` = (" +
	                (1+rate_variation.types.gamma_i.bins) + "," + 
	                (rate_variation.types.gamma_i.cats) + 
	                ",MEAN, 
					(1-`rate_variation.types.gamma_i.inv_p`)*GammaDist(_x_,`rate_variation.types.gamma_i.alpha`,`rate_variation.types.gamma_i.beta`)*(_x_>0), 
					(1-`rate_variation.types.gamma_i.inv_p`)*CGammaDist(_x_,`rate_variation.types.gamma_i.alpha`,`rate_variation.types.gamma_i.beta`)*(_x_>0)+`rate_variation.types.gamma_i.inv_p`, 
					0, 
			  	    1e25,
			  	    (1-`rate_variation.types.gamma_i.inv_p`)*CGammaDist(_x_,`rate_variation.types.gamma_i.alpha`+1,`rate_variation.types.gamma_i.beta`)*(`rate_variation.types.gamma_i.alpha`/`rate_variation.types.gamma_i.beta`)*(_x_>0)
			  	 )");
			  	 
     ExecuteCommands ("GetInformation (info, `rate_variation.types.gamma_i.main`)");
			  	 
     rate_variation.types.gamma_i.q = utility.callFunction (model["original q call-back"],   
                                                           {"0" : model,
                                                            "1" : parameters.quote (namespace)});
       
     rate_variation.types.gamma_i.q [terms.rate_matrix] = parameters.addMultiplicativeTerm (rate_variation.types.gamma_i.q[terms.rate_matrix], rate_variation.types.gamma_i.main, 0);
     return rate_variation.types.gamma_i.q;                                               
}