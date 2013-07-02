/*

    This is a file which does NN interfacing with HBL model structures.

    September 2003, Sergei L Kosakovsky Pond

*/

#include "HYNetInterface.h"
#include "baseobj.h"
#include "parser.h"

_String   ModelTrainNNFlag          ("TRAIN_MODEL_NN"),
          ModelLoadNNFlag             ("LOAD_MODEL_NN"),
          ModelNNTrainingSteps        ("MODEL_NN_TRAINING_STEPS"),
          ModelNNVerificationSample   ("MODEL_NN_VERIFICATION_STEPS"),
          ModelNNPrecision            ("MODEL_NN_PRECISION"),
          ModelNNLoops                ("MODEL_MAX_LOOPS"),
          ModelNNHiddenNodes      ("MODEL_NN_HIDDEN_NODES"),
          ModelNNFile             ("MODEL_NN_FILE_PATH"),
          ModelNNLearningRate     ("MODEL_NN_LEARNING_RATE"),
          ModelNNPersistenceRate  ("MODEL_NN_PERSISTENCE_RATE");



/*______________________________________________________________________*/

void         TrainModelNN   (_String* model, _String* matrix)
{
    _String         errMsg;

    long            modelIdx = modelNames.Find(model);

    _Parameter      verbI;

    checkParameter (VerbosityLevelString, verbI, 0.0);

    char            buffer [128];

    if (modelIdx < 0) {
        errMsg = *model & " did not refer to an existring model";
    } else {
        _Variable*    boundsMatrix              =   FetchVar  (LocateVarByName (*matrix));

        if (boundsMatrix && (boundsMatrix->ObjectClass() == MATRIX)) {
            _Matrix  *    bmatrix               =   (_Matrix*)      boundsMatrix->GetValue ();

            if (bmatrix->IsAStringMatrix() && (bmatrix->GetVDim () == 3)) {
                _Variable*    modelMatrix           =   LocateVar       (modelMatrixIndices.lData[modelIdx]);
                _SimpleList   modelVariableList;
                {
                    _AVLList mvla (&modelVariableList);
                    modelMatrix->ScanForVariables       (mvla, true);
                    mvla.ReorderList();
                }

                if (bmatrix->GetHDim () == modelVariableList.lLength) {
                    // now map model variables to bounds matrix
                    _SimpleList          variableMap;
                    _String             *myName;

                    for (long k = 0; k < modelVariableList.lLength; k++) {
                        myName = ((_FString*)bmatrix->GetFormula(k,0)->Compute())->theString;
                        long      vID    = LocateVarByName (*myName);

                        if (vID < 0) {
                            break;
                        }

                        vID = variableNames.GetXtra (vID);
                        vID = modelVariableList.Find(vID);

                        if (vID < 0) {
                            break;
                        }

                        variableMap << vID;
                    }

                    if (variableMap.lLength == modelVariableList.lLength) {
                        _Matrix     vBounds (variableMap.lLength,2, false, true);

                        long        k2 = 0;

                        for (; k2 < variableMap.lLength; k2++) {
                            _Parameter lb = ((_FString*)bmatrix->GetFormula(k2,1)->Compute())->theString->toNum(),
                                       ub = ((_FString*)bmatrix->GetFormula(k2,2)->Compute())->theString->toNum();

                            if ( ub>lb || k2) {
                                vBounds.Store (k2,0,lb);
                                vBounds.Store (k2,1,ub);
                                if (ub<=lb && vBounds (k2-1,0) <= vBounds (k2-1,1) && (!CheckEqual(vBounds (k2-1,0),0.0) || !CheckEqual(vBounds (k2-1,1),1.0))) {
                                    break;
                                }
                            }

                        }
                        if (k2 == modelVariableList.lLength) {
                            // set up the sampling now
                            _String             fName       = ProcessLiteralArgument (&ModelNNFile,nil);
                            FILE*               nnFile      = doFileOpen (fName.getStr(), "w");
                            if (nnFile) {
                                _Matrix*            modelMatrix = (_Matrix*) LocateVar(modelMatrixIndices.lData[modelIdx])->GetValue();

                                _Parameter          mainSteps,
                                                    checkSteps,
                                                    errorTerm,
                                                    loopMax,
                                                    hiddenNodes,
                                                    absError,
                                                    nn1,
                                                    nn2;

                                long                fullDimension = modelMatrix->GetHDim() * modelMatrix->GetVDim();


                                checkParameter      (ModelNNTrainingSteps,      mainSteps,      10000.0);
                                checkParameter      (ModelNNVerificationSample, checkSteps,     500.0);
                                checkParameter      (ModelNNPrecision,          errorTerm,      0.01);
                                checkParameter      (ModelNNTrainingSteps,      loopMax,        10);
                                checkParameter      (ModelNNHiddenNodes,        hiddenNodes,    5);
                                checkParameter      (ModelNNLearningRate,       nn1,            .3);
                                checkParameter      (ModelNNPersistenceRate,    nn2,            .1);

                                Net**               matrixNet = new Net* [fullDimension] ;

                                for (long i = 0; i < fullDimension; i++) {
                                    checkPointer (matrixNet [i] = new Net (variableMap.lLength,(long)hiddenNodes,1,errorTerm,nn1,nn2,100,200,true));
                                    //matrixNet[i]->verbose = true;
                                }

                                checkPointer        (matrixNet);

                                _List               tIn,
                                                    tOut;

                                FILE*               varSamples = doFileOpen ("variableSamples.out", "w");

                                fprintf (varSamples, "%s" ,LocateVar(modelVariableList.lData[0])->GetName()->getStr());
                                for (long vc = 1; vc < modelVariableList.lLength; vc++) {
                                    fprintf (varSamples, ",%s" ,LocateVar(modelVariableList.lData[variableMap.lData[vc]])->GetName()->getStr());
                                }

                                fprintf (varSamples, "\n");

                                for (long itCount = 0; itCount < loopMax; itCount ++) {
                                    if (verbI > 5) {
                                        snprintf (buffer, sizeof(buffer), "\nNeural Network Pass %ld. Building a training set...\n", itCount);
                                        BufferToConsole (buffer);
                                    }

                                    while   (tIn.countitems() < mainSteps) {
                                        NNMatrixSampler     (0, vBounds, modelVariableList, variableMap, modelMatrix, tIn, tOut);
                                    }

                                    _Matrix             inData (mainSteps, variableMap.lLength, false, true);
                                    _Parameter          *md = inData.theData;

                                    for (long matrixC = 0; matrixC < mainSteps; matrixC++) {
                                        _Parameter  *   ed = ((_Matrix*)tIn (matrixC))->theData;
                                        fprintf (varSamples, "\n%g",*ed);
                                        *md = *ed;
                                        ed++;
                                        md++;
                                        for (long entryC = 1; entryC < variableMap.lLength; entryC++, ed++, md++) {
                                            *md = *ed;
                                            fprintf (varSamples, ",%g", *md);
                                        }
                                    }

                                    tIn.Clear();

                                    if (verbI > 5) {
                                        BufferToConsole ( "Done Building Training Set. Training...\n");
                                    }

                                    long lastDone = 0;

                                    for (long cellCount = 0; cellCount < fullDimension; cellCount++) {
                                        Net* thisCell = matrixNet[cellCount];
                                        _Matrix outVector (mainSteps, 1, false, true);

                                        for (long oc = 0; oc < mainSteps; oc++) {
                                            outVector.theData[oc] = ((_Matrix*)tOut(oc))->theData[cellCount];
                                        }

                                        thisCell->studyAll (inData.theData, outVector.theData, mainSteps);

                                        long    nowDone = (cellCount+1)*100./fullDimension;
                                        if (nowDone > lastDone) {
                                            snprintf (buffer, sizeof(buffer),"%ld%% done\n", lastDone = nowDone);
                                            BufferToConsole (buffer);
                                        }
                                    }
                                    tOut.Clear();

                                    if (verbI > 5) {
                                        BufferToConsole ( "Done Training. Resampling...\n");
                                    }

                                    _PMathObj  tObj = _Constant(0).Time();
                                    _Parameter time1 = tObj->Value(),
                                               time2;

                                    while   (tIn.countitems() < checkSteps) {
                                        NNMatrixSampler     (0, vBounds, modelVariableList, variableMap, modelMatrix, tIn, tOut);
                                    }

                                    absError = 0.0;

                                    DeleteObject (tObj);
                                    tObj = _Constant(0).Time();
                                    time2 = tObj->Value();

                                    if (verbI > 5) {
                                        snprintf (buffer, sizeof(buffer),"Done Resampling in %g seconds. Computing Error...\n", time2-time1);
                                        BufferToConsole (buffer);
                                    }

                                    _Parameter  maxValT,
                                                maxValE;

                                    for (long verCount = 0; verCount < checkSteps; verCount++) {
                                        _Parameter*  inData     = ((_Matrix*)tIn(verCount))->theData,
                                                     *  outData  = ((_Matrix*)tOut(verCount))->theData;

                                        for (long cellCount = 0; cellCount < fullDimension; cellCount++) {
                                            Net         *thisCell = matrixNet[cellCount];

                                            _Parameter estVal     = thisCell->eval(inData)[0],
                                                       trueVal      = outData[cellCount],
                                                       localError;

                                            localError = estVal-trueVal;

                                            if (localError < 0) {
                                                localError = -localError;
                                            }

                                            if (absError < localError) {
                                                maxValT  = trueVal;
                                                maxValE  = estVal;
                                                absError = localError;
                                            }
                                        }
                                    }

                                    DeleteObject (tObj);
                                    tObj = _Constant(0).Time();
                                    time1 = tObj->Value();
                                    DeleteObject (tObj);


                                    if (verbI > 5) {
                                        snprintf (buffer, sizeof(buffer), "Done Error Checking in %g seconds. Got max abs error %g on the pair %g %g\n", time1-time2, absError, maxValT, maxValE);
                                        BufferToConsole (buffer);
                                    }
                                    if (absError <= errorTerm) {
                                        break;
                                    }
                                }

                                if (absError > errorTerm) {
                                    ReportWarning (_String("Couldn't achive desired precision in TrainModelNN. Achieved error of ") & absError);
                                }
                                fclose  (varSamples);
                                fprintf (nnFile,"{{\n\"%s\"", LocateVar(modelVariableList.lData[0])->GetName()->getStr());
                                _Matrix newBounds (modelVariableList.lLength, 2, false, true);
                                if (vBounds(0,0)>vBounds(0,1)) {
                                    newBounds.Store (variableMap.lData[0],0,0.);
                                    newBounds.Store (variableMap.lData[0],1,1.);
                                } else {
                                    newBounds.Store (variableMap.lData[0],0,vBounds(0,0));
                                    newBounds.Store (variableMap.lData[0],1,vBounds(0,1));
                                }
                                for (long varCounter = 1; varCounter < modelVariableList.lLength; varCounter ++) {
                                    fprintf (nnFile,",\n\"%s\"", LocateVar(modelVariableList.lData[varCounter])->GetName()->getStr());
                                    if (vBounds(varCounter,0)>vBounds(varCounter,1)) {
                                        newBounds.Store (variableMap.lData[varCounter],0,0.);
                                        newBounds.Store (variableMap.lData[varCounter],1,1.);
                                    } else {
                                        newBounds.Store (variableMap.lData[varCounter],0,vBounds(varCounter,0));
                                        newBounds.Store (variableMap.lData[varCounter],1,vBounds(varCounter,1));
                                    }
                                }

                                fprintf (nnFile,"\n}}\n");
                                newBounds.toFileStr (nnFile);


                                for (long i2 = 0; i2 < fullDimension; i2++) {
                                    matrixNet[i2]->save(nnFile);
                                    delete matrixNet [i2];
                                }

                                fclose (nnFile);
                                delete              matrixNet;
                            } else {
                                errMsg = _String ("Failed to open ") & fName & " for writing";
                            }
                        } else {
                            errMsg = _String ("Invalid variable bounds in row ") & (k2+1) & " of the bounds matrix";
                        }
                    } else {
                        errMsg = *myName & " was not one of the model parameters";
                    }

                } else {
                    errMsg = *matrix & " must be a have the same number of rows as the number of model parameters";
                }

            } else {
                errMsg = *matrix & " must be a string matrix with 3 columns";
            }
        } else {
            errMsg = *matrix & " was not the identifier of a valid matrix variable";
        }



    }

    if (errMsg.sLength) {
        errMsg = errMsg & _String(" in call to TrainModelNN.");
        WarnError (errMsg);
    }
}

/*______________________________________________________________________*/

void        NNMatrixSampler (long index, _Matrix& bnds, _SimpleList& varList, _SimpleList& reidx, _Matrix* modelMatrix, _List& trainIn, _List& trainOut)
{
    _Parameter lb = bnds (index,0),
               ub = bnds (index,1);

    if (ub<=lb)
        // freq parameter
    {
        ub = 1.;
        lb = 0.;
        for (long k = index-1; k>=0; k--) {
            ub -= LocateVar (varList.lData[reidx.lData[k]])->Value();
            if (bnds (k,1) > bnds (k,0)) {
                break;
            }
        }

        if ((index == varList.lLength - 1)||(bnds (index+1,1) > bnds (index+1,0))) {
            lb = ub;
        }
    }

    _Constant cv (lb + genrand_real2 () * (ub-lb));
    LocateVar (varList.lData[reidx.lData[index]])->SetValue (&cv);

    if (index == varList.lLength - 1) {
        _Matrix*  inMx = new _Matrix (index+1,1,false,true);

        checkPointer (inMx);

        for (long m = 0; m <= index; m++) {
            inMx->Store (m,0,LocateVar (varList.lData[reidx.lData[m]])->Value());
        }

        trainIn << inMx;

        DeleteObject (inMx);

        _Matrix  *em = ((_Matrix*)modelMatrix->ComputeNumeric())->Exponentiate();
        trainOut <<  em;
        DeleteObject (em);
    } else {
        NNMatrixSampler (index+1, bnds, varList, reidx, modelMatrix, trainIn, trainOut);
    }
}
