ncases = 500;

fprintf (stdout, "\nTesting BayesianGraphicalModel discrete graph functionality\n\n");


ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"bayesgraph.ibf");

fprintf (stdout, "Loaded bayesgraph include file\n");


/* import Bayesian network structure and parameters 
	from XMLBIF (XML Bayesian Interchange Format)
	as an associative list
*/

fprintf (stdout, "Import ALARM network from XMLBIF file...");

import_xmlbif ("alarm.xml", "Alarm");

if ( (Rows(Alarm))[0] == "Hypovolemia") {
	fprintf (stdout, "[PASSED]\n");
} else {
	fprintf (stdout, "[FAILED]\n");
}


// adjacency matrix of network
adjmat = list2adjmat(Alarm);



/* this object contains all the info we need to simulate data */
fprintf (stdout, "Simulate ", ncases, " cases from network object...");
sim = simulate_data (Alarm, ncases);
if (Rows(sim) == ncases) {
	fprintf (stdout, "[PASSED]\n");
} else {
	fprintf (stdout, "[FAILED]\n");
}


/* keys of associative list are variable (node) names */
names = Rows(Alarm);

/* a Bayesian Graphical Model object in HyPhy is constructed with 
	a single (associative list) argument
*/
nodes={};
for (i = 0; i < Abs(Alarm); i=i+1)
{
	/* add_discrete_node (	node name, 
							max. # parents, 
							prior sample size, 
							# levels) 
	*/
	nodes[Abs(nodes)] = add_discrete_node (names[i], 2, 0, (Alarm[names[i]])["Levels"]);
}


num_nodes = Abs(nodes);

/* construct BGM */
fprintf (stdout, "Create BGM object...\n");

BayesianGraphicalModel alarm_bgm = (nodes);

GetString (bgm_names_list, BayesianGraphicalModel, -1);	// returns names of all BGMs

// THIS TEST IS BROKEN
/*
lLength = Rows(bgm_names_list) * Columns(bgm_names_list);
for (_i = 0; _i < lLength; _i += 1) {
	if (bgm_names_list[_i] == "alarm_bgm") {
		fprintf (stdout, "[PASSED]\n");
		break;
	}
}
if (_i == lLength) {
	fprintf (stdout, "[FAILED]\n");
}
*/

/*
	Assign data set to BGM.
	attach_data ( 	BGM identifier,
					data matrix,
					Gibbs imputation #steps,
					"		"		burnin,
					"		"		#samples)
 */
 
fprintf (stdout, "Attaching data and caching node scores...");
attach_data ("alarm_bgm", sim, 0, 0, 0);
cache = get_node_score_cache("alarm_bgm");
if (Abs(cache) == 111) {
	fprintf (stdout, "[PASSED]\n");
} else {
	fprintf (stdout, "[FAILED]\n");
}


/* graph structural MCMC */
fprintf (stdout, "RUNNING GRAPH-MCMC\n");

result0 = graph_MCMC ("alarm_bgm", 100000, 100000, 100, 0);


temp = check_edgelist (result0, adjmat, 0.8);
fprintf (stdout, "\tTrue positives = ", temp[0], "\n");
fprintf (stdout, "\tFalse negatives = ", temp[1], "\n");
fprintf (stdout, "\tFalse positives = ", temp[2], "\n");
fprintf (stdout, "\tTrue negatives = ", temp[3], "\n");

sens = temp[0]/(temp[0]+temp[1]);
spec = temp[3]/(temp[2]+temp[3]);

fprintf (stdout, "\tSensitivity (TP/TP+FN) = ", sens, "\n");
fprintf (stdout, "\tSpecificity (TN/TN+FP) = ", spec, "\n");

fprintf (stdout, "Specificity > 75% and specificity > 90% for cutoff = 0.8 ... ");
if (sens > 0.75 && spec > 0.9) {
	fprintf (stdout, "[PASSED]\n");
} else {
	fprintf (stdout, "[FAILED]\n");
}



display_MCMC_chain (result0);
write_edgelist("TestBGM.graphMCMC.edges", result0, num_nodes, 1);
mcmc_graph_to_dotfile("TestBGM.graphMCMC.dot", 0.6, result0, nodes);



/* node order permutation MCMC */

fprintf (stdout, "RUNNING ORDER-MCMC\n");

result1 = order_MCMC ("alarm_bgm", 10000, 10000, 100);

temp = check_edgelist (result1, adjmat, 0.8);
fprintf (stdout, "\tTrue positives = ", temp[0], "\n");
fprintf (stdout, "\tFalse negatives = ", temp[1], "\n");
fprintf (stdout, "\tFalse positives = ", temp[2], "\n");
fprintf (stdout, "\tTrue negatives = ", temp[3], "\n");

sens = temp[0]/(temp[0]+temp[1]);
spec = temp[3]/(temp[2]+temp[3]);

fprintf (stdout, "\tSensitivity (TP/TP+FN) = ", sens, "\n");
fprintf (stdout, "\tSpecificity (TN/TN+FP) = ", spec, "\n");

fprintf (stdout, "Specificity > 75% and specificity > 90% for cutoff = 0.8 ... ");
if (sens > 0.75 && spec > 0.9) {
	fprintf (stdout, "[PASSED]\n");
} else {
	fprintf (stdout, "[FAILED]\n");
}


write_edgelist("TestBGM.orderMCMC.edges", result1, num_nodes, 1);
mcmc_graph_to_dotfile("TestBGM.orderMCMC.dot", 0.6, result1, nodes);





