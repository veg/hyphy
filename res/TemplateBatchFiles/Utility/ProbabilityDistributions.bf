_piConst = 4*Arctan (1);

/*-------------------------------------------------*/

function sampleFromPoisson (lambda)
{
	_sampledValue = Random (0, Exp(lambda));
	_sum		  = 1;
	_term		  = 1;
	_k			  = 0;
	
	while (_sum < _sampledValue)
	{
		_k = _k+1;
		_term = _term * lambda/_k;
		_sum = _sum + _term;
	}
	return _k;
}

/*-------------------------------------------------*/

function sampleFromGamma (alpha,beta)
{
	_sampledValue = Random (0,1);
	FindRoot (_k,CGammaDist(_x_,alpha,beta)-_sampledValue,_x_,0,2500);
	return _k;
}

/*-------------------------------------------------*/

function sampleFromBeta (p,q)
{
	_sampledValue = Random (0,1);
	FindRoot (_k,IBeta(_x_,p,q)-_sampledValue,_x_,0,1);
	return _k;
}

/*-------------------------------------------------*/

function sampleFromNormal ()
{
	/* Box-Muller */
	
	return Sqrt (-2*Log (Random (1e-26,1)))*Cos (2 * _piConst * Random (1e-26,1));
}

/*-------------------------------------------------*/

function linearFit (xy)
{
    data_rows = Rows(xy);
    ss = data_rows;
    sx = +(xy[-1][0]);
    sy = +(xy[-1][1]);
    sxoss = sx/ss;
				
    fitB = 0;
    st2  = 0;
				
    sxx = 0;
    syy = 0;
    sxy = 0;
    ax  = sx/data_rows;
    ay  = sy/data_rows;
				
    for (count = 0; count<data_rows; count += 1)
    {
        xt = xy[count][0]-ax;
        yt = xy[count][1]-ay;
        sxx += (xt)^2;
        syy += (yt)^2;
        sxy += xt*yt;
        t    = xy[count][0]-sxoss;
        st2  = st2+t*t;
        fitB += t*xy[count][1];
    }

    fitB = fitB/st2;
    fitA = (sy-sx*fitB)/ss;
    varA = Sqrt ((1+sx*sx/(ss*st2))/ss);
    varB = Sqrt (1/st2);
    
    return {"Correlation": sxy/Sqrt(sxx*syy),
            "Intercept": fitA,
            "Slope"    : fitB,
            "Var(Intercept)": varA,
            "Var(Slope)": varB};
    
}


