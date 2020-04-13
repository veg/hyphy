LoadFunctionLibrary("libv3/UtilityFunctions.bf");

math.Infinity = 1/0;

/** @module math */

/**
* Converts a float to integer by rounding
* @name math.Int
* @param float
*/
lfunction math.Int (float) {
  return (float + 0.5)$1;
}

/**
* Computes small-sample AIC
* @name math.GetIC
* @param logl
* @param params
* @param samples
*/
lfunction math.GetIC(logl,params,samples) {
  return -2*logl + 2*samples/(samples-params-1)*params;
}

/**
* Computes the LRT and p-value
* @name math.DoLRT
* @param lognull
* @param logalt
* @param df
*/
lfunction math.DoLRT (lognull,logalt,df) {
  lrt = 2 * (logalt - lognull);
  return {
           utility.getGlobalValue("terms.LRT"): lrt__,
           utility.getGlobalValue("terms.p_value") : 1-CChi2 (lrt__, df__)
         }
}

/**
* Computes sum of 1xN vector
* @name math.Sum
* @param {Matrix} _data_vector - 1xN vector
*/
lfunction math.Sum(_data_vector) {

  return +_data_vector;

}

/**
* Returns median average of 1xN vector
* @name math.Median
* @param {Matrix} _data_vector - 1xN vector
*/
lfunction math.Median(_data_vector) {

  count = utility.Array1D (_data_vector);
     // this will also let you handle other non 1xN vectors

  // sort tmp_data_vector
  _tmp_data_vector = _data_vector % 0;

  if (count%2) {
    median = _tmp_data_vector[count$2];
    // $ is integer division, so you don't end up taking index 1.5
  } else {
    counter = count$2-1;
    median = (_tmp_data_vector[counter]+_tmp_data_vector[counter+1])/2;
  }

  return median;

}

/**
* Returns mean average of 1xN vector
* @name math.Mean
* @param {Matrix} _data_vector - 1xN vector
*/
lfunction math.Mean(_data_vector) {

  return math.Sum (_data_vector) / utility.Array1D (_data_vector);

}

/**
* Returns kurtosis of 1xN vector
* @name math.Kurtosis
* @param {Matrix} _data_vector - 1xN vector
*/
lfunction math.Kurtosis(_data_vector) {

  count = utility.Array1D (_data_vector);
  mean = math.Mean(_data_vector);

  /*
  for (_k = 0; _k < count; _k = _k+1) {
    diff += Abs(_data_vector[_k] - mean)^4;
  }*/

  // it is MUCH faster to do "iterator ops" like this
  diff = + _data_vector ["(_MATRIX_ELEMENT_VALUE_ - mean)^4"];
  _moment = diff/count;

  return _moment/math.Variance(_data_vector)^2;

}

/**
* Returns skewness of 1xN vector
* @name math.Skewness
* @param {Matrix} _data_vector - 1xN vector
*/
lfunction math.Skewness(_data_vector) {

  count = utility.Array1D (_data_vector);
  mean = math.Mean(_data_vector);
  diff = + _data_vector ["(_MATRIX_ELEMENT_VALUE_ - mean)^3"];

  _moment = diff/count;
  return _moment/math.Std(_data_vector)^3;

}

/**
* Returns variance of 1xN vector
* @name math.Variance
* @param {Matrix} _data_vector - 1xN vector
*/
lfunction math.Variance(_data_vector) {

  count = utility.Array1D (_data_vector);
  mean = math.Mean(_data_vector);

  diff = + _data_vector ["(_MATRIX_ELEMENT_VALUE_ - mean)^2"];
  return diff/count;

}

/**
* Returns standard deviation of 1xN vector
* @name math.Std
* @param {Matrix} _data_vector - 1xN vector
*/
lfunction math.Std(_data_vector) {
  return Sqrt (math.Variance (_data_vector));
};

/**
* Returns suite of statistical information for a 1xN vector
* @name math.GatherDescriptiveStats
* @param {Matrix} _data_vector - 1xN vector
*/
lfunction math.GatherDescriptiveStats(_data_vector) {


  _count = utility.Array1D (_data_vector);
  _sorted_data_vector = _data_vector % 0;

  _variance = math.Variance(_data_vector);
  _std = Sqrt(_variance);
  _sum = math.Sum(_data_vector);

  _nonnegatives = utility.Filter (_data_vector, "_value_", "_value_ >= 0");

  _dstats = {};
  _dstats[utility.getGlobalValue("terms.math.count")] = _count;
  _dstats[utility.getGlobalValue("terms.math.mean")] = math.Mean(_data_vector);
  _dstats[utility.getGlobalValue("terms.math.median")] = math.Median(_data_vector);
  _dstats[utility.getGlobalValue("terms.math.minimum")] = _sorted_data_vector[0];
  _dstats[utility.getGlobalValue("terms.math.maximum")] = _sorted_data_vector[_count - 1];
  _dstats[utility.getGlobalValue("terms.math._2.5")] = _sorted_data_vector[(_count*0.025)$1];
  _dstats[utility.getGlobalValue("terms.math._97.5")] = _sorted_data_vector[Min((_count*0.975+0.5)$1, _count-1)];
  _dstats[utility.getGlobalValue("terms.math.sum")] = _sum;
  //_dstats["Mean"] = _sum/_count;
  _dstats[utility.getGlobalValue("terms.math.stddev")] = _std;
  _dstats[utility.getGlobalValue("terms.math.variance")] = _variance;
  _dstats[utility.getGlobalValue("terms.math.cov")] = _std*_count/_sum;
  _dstats[utility.getGlobalValue("terms.math.skewness")] = math.Skewness(_data_vector);
  _dstats[utility.getGlobalValue("terms.math.kurtosis")] = math.Kurtosis(_data_vector);
  _dstats [utility.getGlobalValue("terms.math.square_sum")] = math.Sum(_data_vector*2);
  _dstats [utility.getGlobalValue("terms.math.non_negative")]= Abs(_nonnegatives);

  return _dstats;

}

/**
* Returns Holm-Bonferroni corrected p-values
* @name math.HolmBonferroniCorrection
* @param {Dict} ps - a key -> p-value (or null, for not done)

*/
lfunction math.HolmBonferroniCorrection(ps) {
  tested  = utility.Filter (ps, "_value_", "None!=_value_");
  count   = utility.Array1D (tested);
  names   = utility.Keys (tested);
  indexed = {count, 2};
  for (k = 0; k < count; k+=1) {
    indexed [k][0] = tested[names[k]];
    indexed [k][1] = k;
  }
  indexed = indexed % 0;


  for (k = 0; k < count; k+=1) {
    indexed[k][0] = indexed[k][0] * (count - k);
  }

  corrected = {};

  for (k = 0; k < count; k+=1) {
    corrected[names[indexed[k][1]]] = Min (1,indexed[k][0]);
  }
  utility.Extend (corrected, ps);
  return corrected;
};


/**
* Returns Benjamini-Hochberg corrected q-values (FDR)
* @name math.BenjaminiHochbergFDR 
* @param {Dict} ps - a key -> p-value (or null, for not done)

*/

lfunction math.BenjaminiHochbergFDR(ps) {
  tested  = utility.Filter (ps, "_value_", "None!=_value_");
  count   = utility.Array1D (tested);
  names   = utility.Keys (tested);
  indexed = {count, 2};
  for (k = 0; k < count; k+=1) {
    indexed [k][0] = tested[names[k]];
    indexed [k][1] = k;
  }
  indexed = indexed % 0;


  for (k = 0; k < count; k+=1) {
    indexed[k][0] = indexed[k][0] * count/(k+1);
  }

  corrected = {};

  for (k = 0; k < count; k+=1) {
    corrected[names[indexed[k][1]]] = Min (1,indexed[k][0]);
  }
  
  utility.Extend (corrected, ps);
  return corrected;
};


/**
* Returns the range normalized to the lowest value
* @name math.minNormalizedRange
* @param {Matrix || Dictonary} if Dictonary, the values will be used
* @returns {number}

*/
lfunction math.minNormalizedRange(object) {
  if (Type (object) == "AssociativeList") {
    object = utility.Values(object);
  }
  min_value = Min(object, 0);
  if (min_value == 0) {
    return ^"math.Infinity";
  }
  return (Max(object, 0) - min_value ) / min_value;
};


