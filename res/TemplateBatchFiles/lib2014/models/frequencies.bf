LoadFunctionLibrary ("terms.bf");
LoadFunctionLibrary ("parameters.bf");
LoadFunctionLibrary ("../UtilityFunctions.bf");
LoadFunctionLibrary ("model_functions.bf");

function frequencies.equal (model, namespace, datafilter) {
	__N = Abs (model["alphabet"]);
	model[terms.efv_estimate]      = {__N,1}["1/__N"];
	model[terms.efv_estimate_name] = terms.freqs.equal;
	return model;
}

function frequencies.empirical.nucleotide (model, namespace, datafilter) {
	model = frequencies._aux.empirical.singlechar (model, namespace, datafilter);
	model[terms.efv_estimate_name] = terms.freqs.4x1;
	return model;
}

function frequencies.empirical.corrected.CF3x4 (model, namespace, datafilter) {

    __dimension = model.dimension (model); 
    __alphabet = model ["alphabet"];
    
    assert ( Type (model[terms.rate_matrix]) == "Matrix" && Rows (model[terms.rate_matrix]) == __dimension && Columns (model[terms.rate_matrix]) == __dimension,
            "`terms.rate_matrix` must be defined prior to calling frequencies.empirical.corrected.CF3x4");


    GetDataInfo (_givenAlphabet, *datafilter, "CHARACTERS");
    
 	utility.toggleEnvVariable ("COUNT_GAPS_IN_FREQUENCIES", 0);
	HarvestFrequencies (__f, *datafilter, 3,1,1);
	utility.toggleEnvVariable ("COUNT_GAPS_IN_FREQUENCIES", None);
    
    __estimates = frequencies._aux.CF3x4 (__f, model["bases"], __alphabet, model["stop"]);
 	model[terms.efv_estimate]      = __estimates ["codons"];
    __estimates = __estimates["bases"];
    
     for (_rowChar = 0; _rowChar < __dimension; _rowChar +=1 ){
		for (_colChar = _rowChar + 1; _colChar < __dimension; _colChar += 1) {
			
			__diff = models.codon.diff (__alphabet[_rowChar], __alphabet[_colChar]);
            if (None != __diff) {
                (model[terms.rate_matrix])[_rowChar][_colChar] += "*" + (__estimates[__diff["to"]])[__diff["position"]];
                (model[terms.rate_matrix])[_colChar][_rowChar] += "*" + (__estimates[__diff["from"]])[__diff["position"]];
            }
		}	
	}
	
	model[terms.efv_estimate_name] = terms.freqs.CF3x4;	
	return model;
}

function frequencies.mle (model, namespace, datafilter) {
	assert (0, "frequencies.mle is TBD");
}

//--- AUX FUNCTIONS FROM THIS POINT ON ---//

function frequencies._aux.empirical.singlechar (model, namespace, datafilter) {
	utility.toggleEnvVariable ("COUNT_GAPS_IN_FREQUENCIES", 0);
	HarvestFrequencies (__f, *datafilter, 1,1,1);
	utility.toggleEnvVariable ("COUNT_GAPS_IN_FREQUENCIES", None);
	model[terms.efv_estimate] = __f;
	return model;
}


function frequencies._aux.CF3x4 (observed_3x4,base_alphabet,sense_codons, stop_codons) {

   frequencies._aux.CF3x4.p = {};
   
   frequencies._aux.CF3x4.args = {};
   
   for (frequencies._aux.CF3x4.k = 0; frequencies._aux.CF3x4.k < 3; frequencies._aux.CF3x4.k += 1) {
        frequencies._aux.CF3x4.p [frequencies._aux.CF3x4.k] = parameters.generate_sequential_names ("frequencies._aux.CF3x4.p" + frequencies._aux.CF3x4.k,3,None);    
        parameters.declareGlobal (frequencies._aux.CF3x4.p [frequencies._aux.CF3x4.k], None);
        parameters.setRange (frequencies._aux.CF3x4.p [frequencies._aux.CF3x4.k], terms.range01);
        frequencies._aux.CF3x4.args + (Join (",",frequencies._aux.CF3x4.p [frequencies._aux.CF3x4.k]));
   }   
    
    frequencies._aux.CF3x4.args = Join (",", frequencies._aux.CF3x4.args);
    
    frequencies._aux.CF3x4.n = {}; 
    
    for (frequencies._aux.CF3x4.k = 0; frequencies._aux.CF3x4.k < 3; frequencies._aux.CF3x4.k += 1) {
        frequencies._aux.CF3x4.n[frequencies._aux.CF3x4.k] = parameters.generate_attributed_names ("frequencies._aux.CF3x4.n" +frequencies._aux.CF3x4.k,base_alphabet,None);
        parameters.setConstraint (frequencies._aux.CF3x4.n[frequencies._aux.CF3x4.k], 
                                  parameters.helper.stick_breaking (frequencies._aux.CF3x4.p[frequencies._aux.CF3x4.k], observed_3x4[-1][frequencies._aux.CF3x4.k]), 
                                  "global");
    }

   
	frequencies._aux.stop_count = Columns (stop_codons);
    
    frequencies._aux.CF3x4.stop_correction = {};
    
    

    
    for (frequencies._aux.CF3x4.i = 0; frequencies._aux.CF3x4.i < Columns (base_alphabet); frequencies._aux.CF3x4.i += 1) {
        frequencies._aux.CF3x4.stop_correction [base_alphabet[frequencies._aux.CF3x4.i]] = {{"","",""}};
    }
    
    frequencies._aux.CF3x4.denominator = "1";
    
    for (frequencies._aux.CF3x4.i = 0; frequencies._aux.CF3x4.i < frequencies._aux.stop_count; frequencies._aux.CF3x4.i += 1) {
        frequencies._aux.CF3x4.sc = stop_codons[frequencies._aux.CF3x4.i];
        
        (frequencies._aux.CF3x4.stop_correction[frequencies._aux.CF3x4.sc[0]])[0] += 
                "-frequencies._aux.CF3x4.n1_" +  frequencies._aux.CF3x4.sc[1] +
                "*frequencies._aux.CF3x4.n2_" +  frequencies._aux.CF3x4.sc[2];

        (frequencies._aux.CF3x4.stop_correction[frequencies._aux.CF3x4.sc[1]])[1] += 
                "-frequencies._aux.CF3x4.n0_" +  frequencies._aux.CF3x4.sc[0] +
                "*frequencies._aux.CF3x4.n2_" +  frequencies._aux.CF3x4.sc[2];

        (frequencies._aux.CF3x4.stop_correction[frequencies._aux.CF3x4.sc[2]])[2] += 
                "-frequencies._aux.CF3x4.n0_" +  frequencies._aux.CF3x4.sc[0] +
                "*frequencies._aux.CF3x4.n1_" +  frequencies._aux.CF3x4.sc[1];
                
        frequencies._aux.CF3x4.denominator += "-frequencies._aux.CF3x4.n0_" + frequencies._aux.CF3x4.sc[0] +    
                                              "*frequencies._aux.CF3x4.n1_" + frequencies._aux.CF3x4.sc[1] +
                                              "*frequencies._aux.CF3x4.n2_" + frequencies._aux.CF3x4.sc[2];
    }
	
	
	parameters.setConstraint ("frequencies._aux.CF3x4.denominator", frequencies._aux.CF3x4.denominator, 1);
		    
	frequencies._aux.N = {Columns (base_alphabet),3};
	frequencies._aux.res = {};
	frequencies._aux.codons = {Columns (sense_codons), 1};
	
	for (frequencies._aux.CF3x4.i = 0; frequencies._aux.CF3x4.i < Columns (sense_codons); frequencies._aux.CF3x4.i += 1) {
	    frequencies._aux.CF3x4.sc = {3,1};
		for (frequencies._aux.CF3x4.pos = 0; frequencies._aux.CF3x4.pos < 3; frequencies._aux.CF3x4.pos += 1) {
		    frequencies._aux.CF3x4.sc [frequencies._aux.CF3x4.pos] = "frequencies._aux.CF3x4.n"
                             +frequencies._aux.CF3x4.pos+"_"+(sense_codons[frequencies._aux.CF3x4.i])[frequencies._aux.CF3x4.pos];
        }	    
        ExecuteCommands ("frequencies._aux.codons[frequencies._aux.CF3x4.i] := " + Join ("*", frequencies._aux.CF3x4.sc) + "/frequencies._aux.CF3x4.denominator");
	}
	
	for (frequencies._aux.CF3x4.i = 0; frequencies._aux.CF3x4.i < Columns (base_alphabet); frequencies._aux.CF3x4.i += 1) {
	    frequencies._aux.CF3x4.n = base_alphabet[frequencies._aux.CF3x4.i];
	    frequencies._aux.res [frequencies._aux.CF3x4.n] = {3,1};
		for (frequencies._aux.CF3x4.pos = 0; frequencies._aux.CF3x4.pos < 3; frequencies._aux.CF3x4.pos += 1) {
		    
		    frequencies._aux.CF3x4.sc = (frequencies._aux.CF3x4.stop_correction[frequencies._aux.CF3x4.n])[frequencies._aux.CF3x4.pos];
		    
			if (Abs (frequencies._aux.CF3x4.sc)) {
                frequencies._aux.CF3x4.sc = "*(1" + frequencies._aux.CF3x4.sc + ")";		
            }
            
            ExecuteCommands( "frequencies._aux.N[" + frequencies._aux.CF3x4.i + "][" + frequencies._aux.CF3x4.pos + "] := frequencies._aux.CF3x4.n"
                             +frequencies._aux.CF3x4.pos+"_"+frequencies._aux.CF3x4.n
                             +frequencies._aux.CF3x4.sc+"/frequencies._aux.CF3x4.denominator");
                             
            ExecuteCommands ("(frequencies._aux.res[frequencies._aux.CF3x4.n])[frequencies._aux.CF3x4.pos] := frequencies._aux.CF3x4.n"
                             +frequencies._aux.CF3x4.pos+"_"+frequencies._aux.CF3x4.n);
		}
	}
	
	
	
    ExecuteCommands ("Optimize (frequencies._aux.CF3x4.p, frequencies._aux._CF3x4_minimizer( "  +
                     frequencies._aux.CF3x4.args + "))");

	return {"codons" : Eval("frequencies._aux.codons"), "bases" : frequencies._aux.res};
}

function frequencies._aux._CF3x4_minimizer (p11,p12,p13,p21,p22,p23,p31,p32,p33) {
	frequencies._aux._CF3x4_minimizer.error = frequencies._aux.N-observed_3x4;
	return  - (+ frequencies._aux._CF3x4_minimizer.error$frequencies._aux._CF3x4_minimizer.error);
}
