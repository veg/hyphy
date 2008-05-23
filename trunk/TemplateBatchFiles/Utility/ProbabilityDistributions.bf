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