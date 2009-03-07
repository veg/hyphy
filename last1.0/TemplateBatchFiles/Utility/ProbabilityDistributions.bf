
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