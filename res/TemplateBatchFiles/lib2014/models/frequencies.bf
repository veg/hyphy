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

    GetDataInfo (_givenAlphabet, *datafilter, "CHARACTERS");
    
 	utility.toggleEnvVariable ("COUNT_GAPS_IN_FREQUENCIES", 0);
	HarvestFrequencies (__f, *datafilter, 3,1,1);
	utility.toggleEnvVariable ("COUNT_GAPS_IN_FREQUENCIES", None);
    
    frequencies._aux.CF3x4 (__f, model["bases"], model["stop"]);
	
	model[terms.efv_estimate] = __f;
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


function frequencies._aux.CF3x4 (observed_3x4,base_alphabet,stop_codons) {

   frequencies._aux.CF3x4.p = {};
   for (frequencies._aux.CF3x4.k = 0; frequencies._aux.CF3x4.k < 3; frequencies._aux.CF3x4.k += 1) {
        frequencies._aux.CF3x4.p [frequencies._aux.CF3x4.k] = parameters.generate_sequential_names ("frequencies._aux.CF3x4.p" + frequencies._aux.CF3x4.k,3,None);    
        parameters.declareGlobal (frequencies._aux.CF3x4.p [frequencies._aux.CF3x4.k], None);
        parameters.setRange (frequencies._aux.CF3x4.p [frequencies._aux.CF3x4.k], terms.range01);
   }   
    
    frequencies._aux.CF3x4.n = {}; 
    
    for (frequencies._aux.CF3x4.k = 0; frequencies._aux.CF3x4.k < 3; frequencies._aux.CF3x4.k += 1) {
        frequencies._aux.CF3x4.n[frequencies._aux.CF3x4.k] = parameters.generate_attributed_names ("frequencies._aux.CF3x4.n" +frequencies._aux.CF3x4.k,base_alphabet,None);
        parameters.setConstraint (frequencies._aux.CF3x4.n[frequencies._aux.CF3x4.k], 
                                  parameters.helper.stick_breaking (frequencies._aux.CF3x4.p[frequencies._aux.CF3x4.k], observed_3x4[-1][frequencies._aux.CF3x4.k]), 
                                  "global");
    }

   
	frequencies._aux.stop_count = Columns (stop_codons);
    fprintf (stdout, frequencies._aux.stop_count, "\n");
	
	return 0;
	
	charMap   		= {"A":0,"C":1,"G":2, "T": 3};
	revMap			= {{"A","C","G","T"}};
	stopComposition = {stopCount,3};
	
	SDef =""; SDef * 128; SDef * "1";
	
	for (i = 0; i < stopCount; i = i+1)
	{
		SDef * "-";
		for (j = 0; j < 3; j = j+1)
		{
			stopComposition[i][j] = charMap[stopCodons[4*i+j]];
			if (j)
			{
				SDef * "*";
			}
			SDef * ("n"+(j+1)+stopCodons[4*i+j]);
		}
	}
	
	SDef * 0;
	
	ExecuteCommands ("global S:=`SDef`");
	
	SDef = {4,3};
	for (i = 0; i<4; i=i+1)
	{
		for (j = 0; j<3; j=j+1)
		{
			SDef[i][j] = "";
		}
	}
	
	for (k = 0; k < stopCount; k = k+1)
	{
		SDef[stopComposition[k][0]][0] = SDef[stopComposition[k][0]][0] + 
										 "-" +
										 "n2" + revMap[stopComposition[k][1]] +
										 "*n3" + revMap[stopComposition[k][2]];
		SDef[stopComposition[k][1]][1] = SDef[stopComposition[k][1]][1] + 
										 "-" +
										 "n1" + revMap[stopComposition[k][0]] +
										 "*n3" + revMap[stopComposition[k][2]];
		SDef[stopComposition[k][2]][2] = SDef[stopComposition[k][2]][2] + 
										 "-" +
										 "n1" + revMap[stopComposition[k][0]] +
										 "*n2" + revMap[stopComposition[k][1]];
	}

	frequencies._aux.N = {4,3};
	
	for (i = 0; i<4; i=i+1) {
		for (j = 0; j<3; j=j+1)
		{
			if (Abs (SDef[i][j]))
			{
				ExecuteCommands ("N[i][j] := n"+(j+1)+revMap[i] + "*(1" + SDef[i][j] + ")/S;");
			}
			else
			{
				ExecuteCommands ("N[i][j] := n"+(j+1)+revMap[i] + "/S;");			
			}
		}
	}
	

	Optimize (res, frequencies._aux._CF3x4_minimizer(p11,p12,p13,p21,p22,p23,p31,p32,p33));
		
	return {{n1A__,n2A__,n3A__}{n1C__,n2C__,n3C__}{n1G__,n2G__,n3G__}{n1T__,n2T__,n3T__}};

	return None;
}

function frequencies._aux._CF3x4_minimizer (p11,p12,p13,p21,p22,p23,p31,p32,p33) {
	error = N-observed_3x4;
	return  - (+ error$error);
}
