LoadFunctionLibrary("libv3/UtilityFunctions.bf");

/** @module math */

/**
* Computes small-sample AIC
* @name math.getIC
* @param logl
* @param params
* @param samples
*/
lfunction math.getIC(logl,params,samples) {
  return -2*logl + 2*samples/(samples-params-1)*params;
}

/**
* Computes the LRT and p-value
* @name math.doLRT
* @param lognull
* @param logalt
* @param df
*/
lfunction math.doLRT (lognull,logalt,df) {
  lrt = 2 * (logalt - lognull);
  return {
           "LRT"     : lrt,
           "p-value" : 1-CChi2 (lrt, df)
         }
}

/**
* Computes sum of 1xN vector
* @name math.sum
* @param {Matrix} _data_vector - 1xN vector
*/
lfunction math.sum(_data_vector) {

  return +_data_vector;

}

/**
* Returns median average of 1xN vector
* @name math.median
* @param {Matrix} _data_vector - 1xN vector
*/
lfunction math.median(_data_vector) {

  count = utility.array1D (_data_vector);
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
* @name math.mean
* @param {Matrix} _data_vector - 1xN vector
*/
lfunction math.mean(_data_vector) {

  return math.sum (_data_vector) / utility.array1D (_data_vector);

}

/**
* Returns kurtosis of 1xN vector
* @name math.kurtosis
* @param {Matrix} _data_vector - 1xN vector
*/
lfunction math.kurtosis(_data_vector) {

  count = utility.array1D (_data_vector);
  mean = math.mean(_data_vector);

  /*
  for (_k = 0; _k < count; _k = _k+1) {
    diff += Abs(_data_vector[_k] - mean)^4;
  }*/

  // it is MUCH faster to do "iterator ops" like this
  diff = + _data_vector ["(_MATRIX_ELEMENT_VALUE_ - mean)^4"];
  _moment = diff/count;

  return _moment/math.variance(_data_vector)^2;

}

/**
* Returns skewness of 1xN vector
* @name math.skewness
* @param {Matrix} _data_vector - 1xN vector
*/
lfunction math.skewness(_data_vector) {

  count = utility.array1D (_data_vector);
  mean = math.mean(_data_vector);
  diff = + _data_vector ["(_MATRIX_ELEMENT_VALUE_ - mean)^3"];

  _moment = diff/count;
  return _moment/math.std(_data_vector)^3;

}

/**
* Returns variance of 1xN vector
* @name math.variance
* @param {Matrix} _data_vector - 1xN vector
*/
lfunction math.variance(_data_vector) {

  count = utility.array1D (_data_vector);
  mean = math.mean(_data_vector);

  diff = + _data_vector ["(_MATRIX_ELEMENT_VALUE_ - mean)^2"];
  return diff/count;

}

/**
* Returns standard deviation of 1xN vector
* @name math.std
* @param {Matrix} _data_vector - 1xN vector
*/
lfunction math.std(_data_vector) {
  return Sqrt (math.variance (_data_vector));
};

/**
* Returns suite of statistical information for a 1xN vector
* @name math.gather_descriptive_stats
* @param {Matrix} _data_vector - 1xN vector
*/
lfunction math.gather_descriptive_stats(_data_vector) {

  _count = utility.array1D (_data_vector);
  _sorted_data_vector = _data_vector % 0;

  _variance = math.variance(_data_vector);
  _std = Sqrt(_variance);
  _sum = math.sum(_data_vector);

  _nonnegatives = utility.filter (_data_vector, "_value_", "_value_ >= 0");

  _dstats = {};
  _dstats["Count"] = _count;
  _dstats["Mean"] = math.mean(_data_vector);
  _dstats["Median"] = math.median(_data_vector);
  _dstats["Min"] = _sorted_data_vector[0];
  _dstats["Max"] = _sorted_data_vector[_count - 1];
  _dstats["2.5%"] = _sorted_data_vector[(_count*0.025)$1];
  _dstats["97.5%"] = _sorted_data_vector[Min((_count*0.975+0.5)$1, _count-1)];
  _dstats["Sum"] = _sum;
  _dstats["Mean"] = _sum/_count;
  _dstats["Std.Dev"] = _std;
  _dstats["Variance"] = _variance;
  _dstats["COV"] = _std*_count/_sum;
  _dstats["Skewness"] = math.skewness(_data_vector);
  _dstats["Kurtosis"] = math.kurtosis(_data_vector);
  _dstats ["Sq. sum"] = math.sum(_data_vector*2);
  _dstats ["Non-negative"]= Abs(_nonnegatives);

  return _dstats;

}


