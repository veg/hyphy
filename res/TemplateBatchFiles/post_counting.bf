ChoiceList  (response,"Subsitution Bias",1,NO_SKIP,
			 "None","All nucleotide substitutions have the same probability.",
			 "Custom","Nucleotide substitutions are weighted by their relative rates.",
			 "Load","Load a saved matrix of substitution biases");			 
			 
if (response<0)
{
	return;
}

useCustomCountingBias = 0;

function promptRate (pString)
{
	rateValue = -1;
	while (rateValue < 0)
	{
		fprintf (stdout, "\nEnter the rate for ", pString, " subsitutions (>=0):");
		fscanf (stdin,"String",rateValue);
	}
	return rateValue;
}

if (response)
{
	if (response == 1)
	{
		ACRATE = promptRate ("A<->C");
		AGRATE = promptRate ("A<->G");
		ATRATE = promptRate ("A<->T");
		CGRATE = promptRate ("C<->G");
		CTRATE = promptRate ("C<->T");
		GTRATE = promptRate ("G<->T");
		
		_EFV_MATRIX0_ = {{0,ACRATE__,AGRATE__,ATRATE__}
						{ACRATE__,0,CGRATE__,CTRATE__}
						{AGRATE__,CGRATE__,0,GTRATE__}
						{ATRATE__,CTRATE__,GTRATE__,0}};
	}
	else
	{
		SetDialogPrompt ("Load saved bias matrix:");
		fscanf (PROMPT_FOR_FILE,"NMatrix", _EFV_MATRIX0_);
	}

	_EFV_MATRIX1_ = _EFV_MATRIX0_;				
	_EFV_MATRIX2_ = _EFV_MATRIX0_;
	
	useCustomCountingBias = 1;
}

#include "SGEmulator.bf";

useCustomCountingBias = 0;
