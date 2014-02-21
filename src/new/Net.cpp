// copyright Oliver Serang, 2003
// this document may not be duplicated, not redistributed without
// written consent of the author

#include "SerangNet.h"
#include "legacy_parser.h"
//#include <iostream.h>

Net::Net(int innum, int hiddennum, int outnum, _Parameter eps, _Parameter c, _Parameter m,int d, int t, bool v)
{
    init(innum, hiddennum, outnum, eps);
    coef=c;
    mom=m;
    timeout=t;
    verbose=v;
    density=d;
}

Net::Net(int innum, int hiddennum, int outnum, _Parameter eps)
{
    init(innum, hiddennum, outnum, eps);
    coef=.3;
    mom=.1;
    timeout=600;
    verbose=false;
    density=100;
}

void Net::init(int innum, int hiddennum, int outnum, _Parameter eps)
{
//  cout << innum << " " << hiddennum << " " << outnum << endl;
    inNum=innum;
    hiddenNum=hiddennum;
    outNum=outnum;
    in=new Node[innum+1];
    hidden=new Node[hiddennum+1];
    out= new _Parameter[outnum];
    temp= new _Parameter[outnum];
    int k, i;
    for (k=0; k<=innum; k++) {
        in[k].weights=new _Parameter[hiddennum];
        in[k].lastDelta=new _Parameter[hiddennum];
        for (i=0; i<hiddennum; i++) {
            in[k].weights[i]=randReal(.1);
            in[k].lastDelta[i]=0;
//          in[k].weights[i]=.1;
        }
    }
    for (k=0; k<=hiddennum; k++) {
        hidden[k].weights= new _Parameter[outnum];
        hidden[k].lastDelta=new _Parameter[outnum];
        for (i=0; i<outnum; i++) {
            hidden[k].weights[i]=randReal(.1);
            hidden[k].lastDelta[i]=0;
//          hidden[k].weights[i]=.1;
        }
    }
    hidden[hiddenNum].value=1;
    in[inNum].value=1;
    LR=learningRate= .1;
    epsilon= eps;
    momentum=.01;
//  cout << "INIT" << endl;
}

void Net::randomize()
{
    //if (verbose)
    //cout << "RANDING" << endl;
    int k,i;
    for (k=0; k<=inNum; k++) {
        for (i=0; i<hiddenNum; i++) {
            in[k].weights[i]=randReal(.1);
            in[k].lastDelta[i]=0;
        }
    }
    for (k=0; k<=hiddenNum; k++) {
        for (i=0; i<outNum; i++) {
            hidden[k].weights[i]=randReal(.1);
            hidden[k].lastDelta[i]=0;
        }
    }
}

void Net::destroy()
{
    int k;
    for (k=0; k<=inNum; k++) {
        delete in[k].weights;
        delete in[k].lastDelta;
    }
    delete in;
    for (k=0; k<hiddenNum; k++) {
        delete hidden[k].weights;
        delete hidden[k].lastDelta;
    }
    delete hidden;
    delete out;
    delete temp;
}

Net::~Net()
{
    destroy();
}

_Parameter Net::dOdW1(int cell, int w, int outN)
{
    // changed to make linear output
    _Parameter dOdH= hidden[w].weights[outN];
    _Parameter dHdW= (1-hidden[w].value)*hidden[w].value*in[cell].value;
    return dOdH*dHdW;
}

void Net::adjust()
{
    learningRate=1*learningRate;
}

_Parameter Net::dOdW2(int cell, int)
{
    // changed to make linear ouput
    return hidden[cell].value;
}

/*void Net::save(ofstream & fout)
{
    fout << endl;
    fout << inNum << " " << hiddenNum << " " << outNum << endl;
    int k, j;
    for (k=0; k<=inNum; k++)
    {
        for (j=0; j<hiddenNum; j++)
        {
            fout << in[k].weights[j] << " ";
        }
    }
    for (k=0; k<=hiddenNum; k++)
    {
        for (j=0; j<outNum; j++)
        {
            fout << hidden[k].weights[j] << " ";
        }
    }
    fout << endl;
}*/

void Net::save(FILE* f)
{
    fprintf (f, "\n{{%d,%d,%d}}\n{{",inNum, hiddenNum, outNum);
    int k, j;
    for (k=0; k<=inNum; k++)
        for (j=0; j<hiddenNum; j++)
            if (k+j == 0) {
                fprintf (f,"%g",in[k].weights[j]);
            } else {
                fprintf (f,",%g",in[k].weights[j]);
            }

    for (k=0; k<=hiddenNum; k++)
        for (j=0; j<outNum; j++) {
            fprintf (f,",%g",hidden[k].weights[j]);
        }

    fprintf (f,"\n}}");
}

void Net::load(FILE*)
{
    /*inNum       = dim[0];
    hiddenNum = dim[0];
    outNum    = (*dim)(2,0);

    hidden=new Node[hiddenNum+1];
    out= new _Parameter[outNum];
    temp= new _Parameter[outNum];
    int k, j, i;
    for (k=0; k<=inNum; k++)
    {
        in[k].weights=new _Parameter[hiddenNum];
        in[k].lastDelta=new _Parameter[hiddenNum];
        for (i=0; i<hiddenNum; i++)
        {
            in[k].lastDelta[i]=0;
        }
    }
    for (k=0; k<=hiddenNum; k++)
    {
        hidden[k].weights= new _Parameter[outNum];
        hidden[k].lastDelta=new _Parameter[outNum];
        for (i=0; i<outNum; i++)
        {
            hidden[k].lastDelta[i]=0;
        }
    }

    hidden[hiddenNum].value=1;
    in[inNum].value=1;
    LR=learningRate= .25;
    momentum=.25;

    long cntr = 0;
    for (k=0; k<=inNum; k++)
    {
        for (j=0; j<hiddenNum; j++)
        {
            in[k].weights[j] = (*coeff)(0,cntr++);
        }
    }
    for (k=0; k<=hiddenNum; k++)
    {
        for (j=0; j<outNum; j++)
        {
            hidden[k].weights[j] = (*coeff)(0,cntr++);
        }
    }*/
}

/*void Net::load(ifstream & fin)
{
    destroy();
    fin >> inNum >> hiddenNum >> outNum;
    in=new Node[inNum+1];
    hidden=new Node[hiddenNum+1];
    out= new _Parameter[outNum];
    temp= new _Parameter[outNum];
    int k, j, i;
    for (k=0; k<=inNum; k++)
    {
        in[k].weights=new _Parameter[hiddenNum];
        in[k].lastDelta=new _Parameter[hiddenNum];
        for (i=0; i<hiddenNum; i++)
        {
            in[k].lastDelta[i]=0;
        }
    }
    for (k=0; k<=hiddenNum; k++)
    {
        hidden[k].weights= new _Parameter[outNum];
        hidden[k].lastDelta=new _Parameter[outNum];
        for (i=0; i<outNum; i++)
        {
            hidden[k].lastDelta[i]=0;
        }
    }
    hidden[hiddenNum].value=1;
    in[inNum].value=1;
    LR=learningRate= .25;
    momentum=.25;


    for (k=0; k<=inNum; k++)
    {
        for (j=0; j<hiddenNum; j++)
        {
            fin >> in[k].weights[j];
        }
    }
    for (k=0; k<=hiddenNum; k++)
    {
        for (j=0; j<outNum; j++)
        {
            fin >> hidden[k].weights[j];
        }
    }
}*/

void Net::learn(_Parameter * input, _Parameter * output)
{
    int k/*, j*/;

    if (outNum == 1) {
        eval1(input);
    } else {
        eval(input);
    }


    for (k=0; k<outNum; k++) {
        temp[k]=2.*(out[k]-output[k]);
    }

    _Parameter *p1 = output,
                *p2 = out,
                 *p3 = out+outNum,
                  *p4 = temp;

    for (; p2!=p3; p1++,p2++,p4++) {
        *p4 = 2.*(*p2-*p1);
    }


    if (outNum == 1)
        for (k=0; k<=inNum; k++) {
            _Parameter* t1   = in[k].weights,
                        *   t2   = in[k].lastDelta,
                            iv   = in[k].value,
                            * t5   = t1,
                              * t5s  = t5+hiddenNum,
                                * t6   = t2;

            Node  *hn= hidden;

            for (; t5!=t5s; t5++,t6++,hn++) {
                // changed
                _Parameter  delta = hn->value;

                delta   = (1.-delta)*delta*iv* * hn->weights
                          *learningRate * *temp +momentum* *t6;
                *t5 -=  delta;
                *t6 =   delta;
            }
        }
    else
        for (k=0; k<=inNum; k++) {
            _Parameter* t1  = in[k].weights,
                        *   t2  = in[k].lastDelta,
                            * t3  = temp,
                              * t3s = temp+outNum,
                                iv  = in[k].value;

            long  outN = 0;

            for (; t3!=t3s; t3++, outN++) {
                _Parameter*  t5  = t1,
                             *  t5s = t5+hiddenNum,
                                *  t6  = t2;

                Node  *  hn  = hidden;

                for (; t5!=t5s; t5++,t6++,hn++) {
                    // changed
                    _Parameter  delta = hn->value;

                    delta   = (1.-delta)*delta*iv*hn->weights[outN]
                              *learningRate * *t3 +momentum* *t6;
                    *t5 -=  delta;
                    *t6 =   delta;
                }
            }
        }


    /*for (k=0; k<=inNum; k++)
    {
        for (j=0; j<outNum; j++)
        {
            for (int w=0; w<hiddenNum; w++)
            {
            // changed
                _Parameter delta=learningRate*dOdW1(k,w,j)*temp[j]+momentum*in[k].lastDelta[w];
                in[k].weights[w]-=delta;
                in[k].lastDelta[w]=delta;
            }
        }
    }*/

    Node   *hn  = hidden,
            *hns = hidden+hiddenNum+1;

    if (outNum > 1) {
        for (; hn!=hns; hn++) {
            _Parameter *x1  = temp,
                        *x1s = temp+outNum,
                         *x2    = hn->weights,
                          *x3  = hn->lastDelta,
                           hcv  = hn->value;

            for (; x1!=x1s; x1++,x2++,x3++) {
                _Parameter delta    =   learningRate* hcv * *x1 +momentum* *x3;
                *x2 -= delta;
                *x3  = delta;
            }
        }
    } else {
        for (; hn!=hns; hn++) {
            _Parameter delta    = learningRate* hn->value * *temp +momentum* *hn->lastDelta;
            *hn->weights   -= delta;
            *hn->lastDelta  = delta;
        }
    }



    /*for (k=0; k<=hiddenNum; k++)
    {
        for (j=0; j<outNum; j++)
        {
            _Parameter delta=learningRate*dOdW2(k,j)*temp[j]+momentum*hidden[k].lastDelta[j];
            hidden[k].weights[j]-=delta;
            hidden[k].lastDelta[j]=delta;
        }
    }*/
}

void Net::studyAll(_Parameter**input, _Parameter**output, int samp)
{
    _Parameter  max = 1.; // to make sure we enter for loop

    int     rep;

    for (rep=0 ; max>epsilon && rep<timeout; rep++) {
        max=-1;
        for (int count=0; count<density; count++) {
            for (int k=0; k<samp; k++) {
                learn(input[k],output[k]);
            }
        }

        for (int k=0; k<samp; k++) {
            _Parameter err=error();
            if (err>max) {
                max=err;
            }
            learn(input[k],output[k]);
        }

        learningRate=sigmaF(2*(max-.45))*coef;
        momentum=learningRate*mom;
    }

    /*if (verbose)
    {
        if (rep>=timeout)
            cout << "DIVERGE" << endl;
        cout << "Cycles taken: " << rep << endl << "Largest Error: " << max << endl;
    }*/

    cycles=rep;
//  cout << endl;
}

_Parameter Net::studyAll(_Parameter*input, _Parameter*output, int samp)
{
    _Parameter  max = 1.; // to make sure we enter for loop

    int     rep;

    for (rep=0 ; max>epsilon && rep<timeout; rep++) {
        max=-1;
        for (int count=0; count<density; count++) {
            for (int k=0; k<samp; k++) {
                learn(input + k*inNum ,output + k*outNum);
            }
        }
        for (int k=0; k<samp; k++) {
            _Parameter err=error();
            if (err>max) {
                max = err;
            }
            learn(input + k*inNum ,output + k*outNum);
        }

        learningRate=sigmaF(2*(max-.45))*coef;
        momentum=learningRate*mom;
    }

    /*if (verbose)
    {
        if (rep>=timeout)
            cout << "DIVERGE" << endl;
        cout << "Cycles taken: " << rep << endl << "Largest Error: " << max << endl;
    }*/

    cycles=rep;
    return max;
//  cout << endl;
}

_Parameter Net::sum(_Parameter * x, int p)
{
    _Parameter total=0;
    for (int k=0; k<p; k++) {
        total+=x[k];
    }
    //cout << total << endl;
    return total;
}

_Parameter Net::error()
{
    _Parameter tempTotal=0;
    for (int y=0; y<outNum; y++) {
        tempTotal+=fabs(temp[y] * .5);
    }
    // this is the Mean Squared Error
    return tempTotal;
}

bool Net::accurate(_Parameter **, _Parameter**output, int samp)
{
    for (int k=0; k<samp; k++) {
        if (!within(out,output[k])) {
            return false;
        }
    }
    return true;
}

bool Net::within(const _Parameter * pred, const _Parameter * act) const
{
    for (int k=0; k<outNum; k++) {
        if (fabs(pred[k]-act[k])>epsilon) {
            return false;
        }
    }
    return true;
}

const _Parameter * Net::eval(_Parameter * input)
{
    int     k,
            i;

    /*for (k=0; k<inNum; k++)
    {
        in[k].value=input[k];
    }*/

    _Parameter* in1 = input, *in2 = input+inNum;
    Node*   np  = in;

    while (in1<in2) {
        np->value = *in1;
        in1++;
        np++;
    }

    for (k=0; k<hiddenNum; k++) {
        _Parameter total = 0.;
        for (i=0; i<=inNum; i++) {
            total+=in[i].value*in[i].weights[k];
        }
        hidden[k].value= sigmaF(total);
    }

    for (k=0; k<outNum; k++) {
        _Parameter total = 0.;
        for (i=0; i<=hiddenNum; i++) {
            total+=hidden[i].value*hidden[i].weights[k];
        }
        out[k]= total;
    }
    return out;
}

const _Parameter * Net::eval1(_Parameter * input)
{
    int     k,
            i;

    _Parameter* in1 = input,
                *in2 = input+inNum,
                 temp [100]; // hack for serial access

    Node*   np  = in;

    while (in1<in2) {
        np->value = *in1;
        in1++;
        np++;
    }

    np = in;
    Node * nps = in+inNum+1;

    in2 = temp+hiddenNum;

    for (in1 = temp; in1 < in2; in1++) {
        *in1 = 0.0;
    }

    for (; np<nps; np++) {
        _Parameter iv  = np->value,
                   *w1  = np->weights;

        for (in1 = temp; in1<in2; in1++, w1++) {
            *in1 += iv * *w1;
        }
    }

    for (k=0; k<hiddenNum; k++) {
        hidden[k].value = sigmaF (temp[k]);
    }

    /*for (k=0; k<hiddenNum; k++)
    {
        _Parameter total = 0.;
        for (i=0; i<=inNum; i++)
        {
            total+=in[i].value*in[i].weights[k];
        }
        hidden[k].value= sigmaF(total);
    }*/

    {
        _Parameter total = 0.;
        for (i=0; i<=hiddenNum; i++) {
            total+=hidden[i].value**hidden[i].weights;
        }
        out[0]= total;
    }

    return out;
}

