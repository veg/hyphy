ExecuteAFile 			("Kernel_support.ibf");
read_kernel_matrix 		(0);

SetDialogPrompt 		("Write component variance explained to:");
DEFAULT_FILE_SAVE_NAME = "KPCA_components.csv";
fprintf 		 		(PROMPT_FOR_FILE,CLEAR_FILE,KEEP_OPEN,"Component,% Explained Variance");
varianceCompPath		= LAST_FILE_PATH;

SetDialogPrompt 		  ("Write rotated datapoints to:");
DEFAULT_FILE_SAVE_NAME  = "KPCA_projections.csv";
fprintf 		 		  (PROMPT_FOR_FILE,CLEAR_FILE,KEEP_OPEN,"Point");
projectionPath			= LAST_FILE_PATH;

kernelPCAOut 			= perform_kernel_PCA (kernel_matrix);

normEV					= kernelPCAOut["E-Vectors"];
components				= Rows (normEV);
eval					= kernelPCAOut["E-Values"];
rotated 				= kernelPCAOut["Rotated"];

fprintf (stdout, "Found ", components, " significant components\n");

columnHeaders = {1,components};

for (k=0; k<components; k=k+1)
{
	fprintf (varianceCompPath, "\n", k+1, ",", 100*eval[k]);
	fprintf (projectionPath, ",Component_", k+1);
	columnHeaders[k] = "Component "+(k+1) + "[" + (100*eval[k]) + "%]";
}

fprintf (varianceCompPath, CLOSE_FILE);

for (k=0; k<points; k=k+1)
{
	fprintf (projectionPath, "\n", names[k+1]);
	for (k2=0; k2<components; k2=k2+1)
	{
		fprintf (projectionPath, ",", rotated[k][k2]);
	}
}

fprintf (projectionPath,CLOSE_FILE);

OpenWindow (CHARTWINDOW,{{"PCA projections"}
			{"columnHeaders"}
			{"rotated"}
			{"Scatterplot"}
			{columnHeaders[0]}
			{columnHeaders[1]}
			{"1st component"}
			{""}
			{"2nd component"}
			{"0"}
			{""}
			{"-1;-1"}
			{"10;1.309;0.785398"}
			{"Times:12:0;Times:10:0;Times:12:2"}
			{"0;0;16777215;16777215;0;0;6579300;11842740;13158600;14474460;0;3947580;16777215;16711680;6845928;16771158;2984993;9199669;7018159;1460610;16748822;11184810;14173291"}
			{"16,0,0"}
			},
			"671;681;70;70");
