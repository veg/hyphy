//------------------------------------------------------------------------------------------
function runAnalysisOnManyFiles (pathList, analysisName, verboseFlag, callbackInput, callbackOutput)
{
	assert (MPI_NODE_COUNT > 1, "This analysis requires an MPI environment");
	_MPI_NODE_STATUS = {MPI_NODE_COUNT-1,1}["-1"]; 

	for (_fileLine = 0; _fileLine < Columns (pathList); _fileLine += 1)
	{
		_filePath = pathList[_fileLine];
		if ((!_filePath) == 0)
		{
			if (verboseFlag)
			{
				fprintf (stdout, "[MPITools: Warning] Filepath ", pathList[_fileLine], " couldn't be open for reading and will be skipped\n");
			}
			continue;
		}
		_analysisOptions = Eval("`callbackInput` (_fileLine, pathList[_fileLine])");
		if (Abs(_analysisOptions) == 0) {
                        if (verboseFlag)
                        {
                                fprintf (stdout, "[MPITools: Warning] Filepath ", pathList[_fileLine], " was marked as cached and will be skipped\n");
                        }
			continue;
		}
		_SendAnMPIJob   (_fileLine, analysisName, _analysisOptions, callbackOutput,verboseFlag);
	}
	
	_fileLine = +(_MPI_NODE_STATUS["_MATRIX_ELEMENT_VALUE_>=0"]);
	while (_fileLine > 0)
	{
		_ReceiveAnMPIJob (callbackOutput, verboseFlag);
		_fileLine = _fileLine - 1;
	}

	return 0;
}

//------------------------------------------------------------------------------------------
function _SendAnMPIJob (jobNumber, fileID, options, callbackOutput, verboseFlag)
{
	for (_mpiNode = 0; _mpiNode < MPI_NODE_COUNT-1; _mpiNode += 1)
	{
		if (_MPI_NODE_STATUS[_mpiNode] < 0)
		{
			break;
		}
	}
	if (_mpiNode == MPI_NODE_COUNT-1)
	{
		_mpiNode = _ReceiveAnMPIJob (callbackOutput, verboseFlag);
	}
	
	_MPI_NODE_STATUS[_mpiNode] = jobNumber + 1;
	if (verboseFlag)
	{
		fprintf (stdout, "[MPITools] Sending filepath ", pathList[jobNumber], " (ID ", jobNumber+1, ") to node ", _mpiNode+1, "\n");
	}
	MPISend (_mpiNode+1,"GLOBAL_FPRINTF_REDIRECT = \"/dev/null\"; LoadFunctionLibrary (\"`fileID`\"," + options + ")");
	return 0;
}

//------------------------------------------------------------------------------------------
function _ReceiveAnMPIJob (callbackOutput, verboseFlag)
{
	MPIReceive (-1,fromNode,result);
	fromNode += (-1);
	
	doneID   = _MPI_NODE_STATUS[fromNode]-1;
	_MPI_NODE_STATUS[fromNode] = -1;

	if (verboseFlag)
	{
		fprintf (stdout, "[MPITools] Received filepath ", pathList[doneID], " (ID ", doneID+1, ") from node ", fromNode+1, "\n");
	}
	
	returnAVL = Eval(result);
	
	if (Abs(callbackOutput))
	{
		ExecuteCommands ("`callbackOutput`(doneID,returnAVL,pathList[doneID])");
	}
	
	return fromNode;
}
