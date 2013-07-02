/* 
	Reproduce Knudsen and Hein's grammar for step-loop RNA,
	termed 'G6' in Dowell and Eddy (2004) BMC Bioinformatics
		S -> LS | L
		L -> >F< | .
		F -> >F< | LS
	where > < correspond to a base-pair (base and its complement)
	and . corresponds to an unpaired base
*/


/* load NLP function macros */
#include "SCFG.ibf";

/* 
	Anatomy of a global variable declaration:
	------------------------------------------------------------------
		global [variable name] = [initial value];
		[variable name] :> [is constrained to be greater than value];
		[variable name] :< [is constrained to be less than value];
*/

/* probability of chain extension */
global extend = 1.8;		extend :> 0.;	extend :< 100.;
global stop = 0.01;			stop :> 0.;		stop :< 100.;
global dS := extend + stop;

/* probability of steming */
global grow_L = 0.01; 	grow_L :> 0.;	grow_L :< 100.;
global stop_L = 1.99;		stop_L :> 0.;	stop_L :< 100.;
global dL := grow_L + stop_L;

global grow_F = 1.32; 	grow_F :> 0.;	grow_F :< 100.;
global stop_F = 0.44;		stop_F :> 0.;	stop_F :< 100.;
global dF := grow_F + stop_F;



/*
	How to declare a production rule:
	-------------------------------------------------------
		A branching rule assumes the form X->YZ.  The non-terminal
		symbols X, Y, and Z can take any string token in HBL so 
		long as your notation is consistent.  
		
		e.g., letting X="foo", Y="bar", and Z="Baz",
		add_branching_rule("foo","bar","Baz","");
		
		Empty quotes indicate that the rule is deterministic.
		Otherwise the quotes must contain a string representation
		of a formula acting on declared global variables.
		
		Any operator implemented in HBL can be used:
			arithmetic		+ - * /
			polynomial		^ Sqrt
			trigonometric	Sin Cos Tan
			transcendental	Exp Log
*/


/* S production rules */
add_branching_rule("S","L","S","extend/dS");
add_branching_rule("S",">","F<","stop/dS*grow_L/dL");	/* enforce CNF */
add_terminal_rule("S",".","stop/dS*stop_L/dL");


/* L production rules */
add_branching_rule("L",">","F<","grow_L/dL");
add_terminal_rule("L",".","stop_L/dL");


/* F production rules */
add_branching_rule("F",">","F<","grow_F/dF");
add_branching_rule("F","L","S","stop_F/dF");
add_branching_rule("F<","F","<","");


add_terminal_rule(">",">","");
add_terminal_rule("<","<","");


/* output_graphviz(0); */

make_SCFG ("G6");
attach_corpus("G6");

//train_grammar("G6", "cg");
//report_grammar("G6");

/*parse_corpus("G6");*/
/*
sim = simulate_corpus("G6", 10);
fprintf (stdout, sim, "\n");
*/
