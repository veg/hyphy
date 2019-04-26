AnalysisCount   = 10;
console.log('Im here');

AnalysisChoices = {AnalysisCount,2};

AnalysisChoices [0][0] = "CountSubstitutions.bf";
AnalysisChoices [0][1] = "Select a genetic code; count the number of synonymous and nonsynonymous substitutions by codon position and by transition/transversion.";

AnalysisChoices [1][0] = "NeutralExpectation.bf";
AnalysisChoices [1][1] = "Estimate the proportion of synononyous and nonsynonymous substitutions expected under neutral evolution for varying base compositions and transition/transversion biases.";

AnalysisChoices [2][0] = "EffectOfTopology.bf";
AnalysisChoices [2][1] = "Investigate the effect of incorrect tree topology on the global estimates of the dN/dS ratio.";

AnalysisChoices [3][0] = "WhatsInTheMean.bf";
AnalysisChoices [3][1] = "A demonstration of how vastly different rate distributions with the same mean can not be distinguished by a simple mean comparison.";

AnalysisChoices [4][0] = "dSdN.bf";
AnalysisChoices [4][1] = "Compute and display dS and dN trees based on a codon model fit.";

AnalysisChoices [5][0] = "NucleotideBiases.bf";
AnalysisChoices [5][1] = "Investigate the effect of nucleotide bias corrections on the global estimates of the dN/dS ratio.";

AnalysisChoices [6][0] = "LRT.bf";
AnalysisChoices [6][1] = "Test the hypothesis of neutrality using a likelihood ratio test and parametric bootstrap to estimate the distribution of the test statistic.";

AnalysisChoices [7][0] = "ErrorEstimates.bf";
AnalysisChoices [7][1] = "Fit a local MG94 model to an alignment and obtain error estimates for the omega ratio on every branch using three different methods.";

AnalysisChoices [8][0] = "LocalvsGlobal.bf";
AnalysisChoices [8][1] = "Use a LRT to test for variation in omega between tree branches.";

AnalysisChoices [9][0] = "BranchAPriori.bf";
AnalysisChoices [9][1] = "Decide if a pre-selected branch (or branches) evolve non-neutrally. One of the approaches is to use a uniform background, and the second one is to use a general background.";

ChoiceList (runWhat, "Example Analysis", 1, SKIP_NONE, AnalysisChoices);

if (runWhat>=0)
{
	console.log('now im in the if');
	fprintf (stdout, 'AnalysisChoices[runWhat][0]: ', AnalysisChoices[runWhat][0], '\n');
	ExecuteAFile ("./"+AnalysisChoices [runWhat][0]);
}
