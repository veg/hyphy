/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2009
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon              (apoon@cfenet.ubc.ca)

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#include "baseobj.h"
#include "list.h"
#include "matrix.h"
#include "datasetfilter.h"
#include "datasetfilternumeric.h"

_DataSetFilterNumeric::_DataSetFilterNumeric (_Matrix* freqs, _List& values, _DataSet *ds, long cc)
{
    unitLength      = 1;
    categoryCount   = cc;

    SetData         (ds);

    _SimpleList     baseFreqs;

    freqs->ConvertToSimpleList (baseFreqs);
    dimension =                ((_Matrix*)values(0))->GetVDim();

    theNodeMap.Populate         (ds->GetNames().lLength,0,1);
    theOriginalOrder.Populate   (((_Matrix*)values(0))->GetHDim()/categoryCount,0,1);

    //theMap.Populate           (theFrequencies.lLength,0,1);
    //theOriginalOrder.Populate     (theFrequencies.lLength,0,1);
    //duplicateMap.Populate     (theFrequencies.lLength,0,1);


    /*CreateMatrix (&probabilityVectors, theNodeMap.lLength, shifter,false,true, false);

    _Parameter   *storeHere = probabilityVectors.theData;
    for (long spec = 0; spec < theNodeMap.lLength; spec++)
    {
        _Matrix * specMatrix = (_Matrix*)values(spec);
        for (long site = 0; site < theFrequencies.lLength; site++)
            for (long state = 0; state < dimension; state++,storeHere++)
                //probabilityVectors.theData [shifter*spec + site*dimension+state]
                *storeHere= specMatrix->theData[site*dimension+state];
    }*/

    _List        siteScores;
    _AVLListXL   siteIndices(&siteScores);

    duplicateMap.RequestSpace (baseFreqs.lLength+1);

    //bool       startD = false;

    char buffer[255];

    for (long site =0; site <baseFreqs.lLength; site++) {
        _Parameter      testV = 0.0;

        for (long k=0; k<theNodeMap.lLength; k++) // sweep down the columns
            for (long state = 0; state < dimension; state++) {
                testV += ((_Matrix*)(((_Matrix**)values.lData)[k]))->theData[site*dimension+state];
            }

        snprintf (buffer, sizeof(buffer), "%20.18g", testV);
        _String     testS (buffer);
        long        f = siteIndices.Find (&testS);

        _SimpleList * sameScore = nil;

        if (f>=0) {
            sameScore = (_SimpleList*)siteIndices.GetXtra (f);
            for (long k = 0; k<sameScore->lLength; k++) {
                bool fit = true;
                f        = sameScore->lData[k];


                for (long spec=0; spec<theNodeMap.lLength && fit; spec++) { // sweep down the columns
                    _Matrix * specMatrix =(_Matrix*)(((_Matrix**)values.lData)[spec]);
                    for (long state = 0; state < dimension; state++)
                        if (specMatrix->theData[site*dimension+state]!=specMatrix->theData[theMap.lData[f]*dimension+state]) {
                            fit = false;
                            break;
                        }
                }

                if (fit) {
                    theFrequencies[f]+=baseFreqs[site];
                    duplicateMap<<f;
                    f = 0;
                    break;
                } else {
                    f = -1;
                }
            }
        }
        if (f==-1) {
            if (!sameScore) {
                checkPointer        (sameScore = new _SimpleList);
                if (siteIndices.Insert  (testS.makeDynamic(),(long)sameScore,false) < 0) {
                    _String eh ("WTF?");
                    StringToConsole(eh);
                }
            }

            (*sameScore) << theFrequencies.lLength;

            duplicateMap<<theFrequencies.lLength;
            theFrequencies<<baseFreqs[site];
            theMap<<site;
        }
    }

    siteIndices.Clear(true);
    shifter         = theFrequencies.lLength*dimension;
    categoryShifter = shifter*theNodeMap.lLength;

    CreateMatrix (&probabilityVectors, theNodeMap.lLength, shifter*categoryCount,false,true, false);
    _Parameter   *storeHere    = probabilityVectors.theData;

    long      refShifter = 0;
    for (long cc = 0; cc < categoryCount; cc++, refShifter += theOriginalOrder.lLength * dimension) {
        for (long spec = 0; spec < theNodeMap.lLength; spec++) {
            _Matrix * specMatrix = (_Matrix*)values(spec);
            for (long site = 0; site < theFrequencies.lLength; site++)
                for (long state = 0; state < dimension; state++,storeHere++) {
                    *storeHere = specMatrix->theData[refShifter + theMap.lData[site]*dimension+state];
                }
        }
    }
}

BaseRef _DataSetFilterNumeric::makeDynamic (void)
{
    _DataSetFilterNumeric * r = new _DataSetFilterNumeric();
    checkPointer            (r);
    r->CopyFilter           (this);
    r->probabilityVectors.Duplicate(&probabilityVectors);
    return r;
}

_Parameter * _DataSetFilterNumeric::getProbabilityVector (long spec, long site, long categoryID)
{
    return probabilityVectors.theData + categoryID * categoryShifter + spec * shifter + site * dimension;
}

bool     _DataSetFilterNumeric::CompareTwoSites (unsigned long, unsigned long, unsigned long)
{
    return false;
}


