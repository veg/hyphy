LoadFunctionLibrary("libv3/convenience/math.bf");

_data_vector = {
  {2}
  {2}
  {4}
  {1}
};

// sum
assert(math.Sum(_data_vector) == 9,"sum == 9 failed");

// mean
assert(math.Mean(_data_vector) == 2.25,"mean == 2.25 failed");

// median
assert(math.Median(_data_vector) == 2,"median == 2 failed");

// standard deviation
assert(Format(math.Std(_data_vector), 4, 4) == "1.0897","std == 1.0897 failed");

// variance 
assert(Format(math.Variance(_data_vector), 4, 4) == "1.1875","variance == 1.1875 failed");

// kurtosis
assert(Format(math.Kurtosis(_data_vector), 3, 3) == "2.097","kurtosis == 2.097 failed");

// skewness 
assert(Format(math.Skewness(_data_vector), 3, 3) == "0.652","skewness == 0.652 failed");

// gather descriptive stats
assert(Abs(math.GatherDescriptiveStats(_data_vector)) == 15, "not all stats accounted for");

fprintf(stdout, math.GatherDescriptiveStats(_data_vector));
