// Oliver Serang
// 
// Net class implements a one hidden layer Neural Net with a
// variable number of inputs, outputs, and hidden Nodes.
// 
// 
// copyright Oliver Serang 2003

// Modified to 

#ifndef _NET_H
#define _NET_H

#include <math.h>
#include <stdlib.h>
//#include <fstream.h>

#include "baseobj.h"



struct Node
{
	_Parameter value;
	_Parameter*weights;
	_Parameter*lastDelta;
};

class Net
{
	public:
		Net () {}
		Net (int innum, int hiddennum, int outnum, _Parameter eps, _Parameter c, _Parameter m, int d, int t, bool v);
		Net (int innum, int hiddennum, int outnum, _Parameter eps);
		
		
		void learn(_Parameter*input, _Parameter*output);
		const _Parameter* eval	(_Parameter*input);
		const _Parameter* eval1	(_Parameter*input);
		bool within(const _Parameter *, const _Parameter *) const;
		bool accurate(_Parameter**input, _Parameter**output, int samp);
		_Parameter bruteDelta;
		void   studyAll(_Parameter**input,_Parameter**output,int samp);
		_Parameter studyAll(_Parameter*input,_Parameter*output,int samp);
		void destroy();
		void randomize();
		void init(int, int, int, _Parameter);
		~Net();
		//void save(ofstream & fout);
		//void load(ifstream & in);
		void save(FILE*);
		void load(FILE*);
		
		
		int cycles,
			timeout,
			density;
			
		bool verbose;
		
		_Parameter coef,
				   mom;
	private:
		_Parameter error();
		Node*in;	// input nodes of the network
		Node*hidden;	// hidden nodes of the network
		_Parameter*out;	// output of the network
		int inNum, hiddenNum, outNum;
		_Parameter*temp;
		_Parameter dOdW1(int , int, int);
		_Parameter dOdW2(int , int);
		_Parameter sum(_Parameter*x,int p);
		_Parameter LR;
		_Parameter learningRate;
		_Parameter momentum;
		_Parameter epsilon;
		void adjust();
		
	inline _Parameter sigmaF(const _Parameter x)
	{
		return 1./(1.+exp(-x));
	}
	
	inline _Parameter randReal(_Parameter magnitude)
	{
		// return a real # in [-magnitude, magnitude]
		//_Parameter x=rand()%1000;
		//x=x-500;
		// x is now an integer in [-500,500]
		//x/=500;
		// x is now a real [-1,1]
		//x*=magnitude;
		//return x;
		return (genrand_int32 () - 2147483648.0) * magnitude / 2147483648.0;
	}
};

#endif
