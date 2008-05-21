function GatherDescriptiveStats (_dataVector)
{
	count = Rows (_dataVector);
	if (count == 1)
	{
		_dataVector = Transpose (_dataVector);
		count = Rows (_dataVector);
	}

	sum  = 0;
	sum2 = 0;
	sum3 = 0;
	sum4 = 0;
	
	_dataVector = _dataVector%0;
	
	data_min  = _dataVector [0];
	data_max  = _dataVector [count-1];
	data_25   = _dataVector [(count*0.025)$1];
	data_975  = _dataVector [Min((count*0.975+0.5)$1,count-1)];
	
	if (count%2)
	{
		median = _dataVector[count/2];
	}
	else
	{
		counter = count/2-1;
		median = (_dataVector[counter]+_dataVector[counter+1])/2;
	}

	for (counter=0; counter < count; counter = counter+1)
	{
		term = _dataVector [counter];
		sum  = sum+term;
		sum2 = sum2+term*term;
		sum3 = sum3+term^3;
		sum4 = sum4+term^4;
	}
	
	_dstats = {};

	_dstats ["Count"] 		= count;
	_dstats ["Mean"] 		= sum/count;
	_dstats ["Median"] 		= median;
	_dstats ["Min"] 		= data_min;
	_dstats ["Max"] 		= data_max;
	_dstats ["2.5%"] 		= data_25;
	_dstats ["97.5%"] 		= data_975;
	_dstats ["Sum"] 		= sum;
	_dstats ["Sq. sum"] 	= sum2;

	if (count > 1)
	{
		counter = (sum2-sum*sum/count)/(count-1);
		
		if (count > 3)
		{
			sum4 = (-6*sum^4+12*count*sum^2*sum2-3(count-1)*count*sum2^2-4*count(count+1)*sum*sum3+count^2*(count+1)*sum4)/(count*(count-1)*(count-2)*(count-3));
			sum4 = sum4/counter^2+3;
		}
		else
		{
			sum4 = 0;
		}	

		if (count > 2)
		{
			sum3 = (2*sum^3-3*count*sum*sum2+count^2*sum3)/(count*(count-1)*(count-2))/Sqrt(counter^2);
		}
		else
		{
			sum3 = 0;
		}	
		
		_dstats["Std.Dev"]   = Sqrt(counter);
		_dstats["Variance"]  = counter;
		_dstats["COV"] 		 = Sqrt (counter)*count/sum;
		_dstats["Skewness"] = sum3;
		_dstats["Kurtosis"] = sum4;
		
	}
	else
	{
		_dstats["Std.Dev"]   = "N/A";
		_dstats["Variance"]  = "N/A";
		_dstats["COV"] 		 = "N/A";
		_dstats["Skewness"]  = "N/A";
		_dstats["Kurtosis"]  = "N/A";
	}
	return _dstats;
}