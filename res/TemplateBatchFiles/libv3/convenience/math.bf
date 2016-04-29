LoadFunctionLibrary("libv3/UtilityFunctions.bf");

/**
 * Computes AIC
 * @param logl
 * @param params
 * @param samples
 */
lfunction math.getIC(logl, params, samples) {
    return -2 * logl + 2 * samples / (samples - params - 1) * params;
}

/**
 * Computes sum of 1xN vector
 * @param {Matrix} _data_vector - 1xN vector
 */
lfunction math.sum(_data_vector) {

    count = Rows(_data_vector);
    sum = 0;

    if (count == 1) {
        _data_vector = Transpose(_data_vector);
        count = Rows(_data_vector);
    }

    for (_k = 0; _k < count; _k = _k + 1) {
        sum += _data_vector[_k];
    }

    return sum;

}

/**
 * Returns median average of 1xN vector
 * @param {Matrix} _data_vector - 1xN vector
 */
lfunction math.median(_data_vector) {

    count = Rows(_data_vector);
    median = 0;

    // sort tmp_data_vector
    _tmp_data_vector = _data_vector % 0;

    if (count % 2) {
        median = _tmp_data_vector[count / 2];
    } else {
        counter = count / 2 - 1;
        median = (_tmp_data_vector[counter] + _tmp_data_vector[counter + 1]) / 2;
    }

    return median;

}

/**
 * Returns mean average of 1xN vector
 * @param {Matrix} _data_vector - 1xN vector
 */
lfunction math.mean(_data_vector) {

    count = Rows(_data_vector);
    mean = 0;

    if (count == 1) {
        _data_vector = Transpose(_data_vector);
        count = Rows(_data_vector);
    }

    sum = math.sum(_data_vector);
    mean = sum / count;
    return mean;

}

/**
 * Returns kurtosis of 1xN vector
 * @param {Matrix} _data_vector - 1xN vector
 */
lfunction math.kurtosis(_data_vector) {

    count = Rows(_data_vector);
    diff = 0;
    _moment = 0;
    _std = 0;

    if (count == 1) {
        _data_vector = Transpose(_data_vector);
        count = Rows(_data_vector);
    }

    mean = math.mean(_data_vector);

    for (_k = 0; _k < count; _k = _k + 1) {
        diff += Abs(_data_vector[_k] - mean) ^ 4;
    }

    _moment = diff / count;
    _std = math.variance(_data_vector) ^ 2;

    return _moment / _std;

}

/**
 * Returns skewness of 1xN vector
 * @param {Matrix} _data_vector - 1xN vector
 */
lfunction math.skewness(_data_vector) {

    count = Rows(_data_vector);
    diff = 0;
    _moment = 0;
    _std_cubed = 0;

    if (count == 1) {
        _data_vector = Transpose(_data_vector);
        count = Rows(_data_vector);
    }

    mean = math.mean(_data_vector);

    for (_k = 0; _k < count; _k = _k + 1) {
        diff += (_data_vector[_k] - mean) ^ 3;
    }

    _moment = diff / count;
    _std_cubed = math.std(_data_vector) ^ 3;

    return _moment / _std_cubed;

}

/**
 * Returns variance of 1xN vector
 * @param {Matrix} _data_vector - 1xN vector
 */
lfunction math.variance(_data_vector) {

    count = Rows(_data_vector);
    diff = 0;

    if (count == 1) {
        _data_vector = Transpose(_data_vector);
        count = Rows(_data_vector);
    }

    mean = math.mean(_data_vector);

    for (_k = 0; _k < count; _k = _k + 1) {
        diff += Abs(_data_vector[_k] - mean) ^ 2;
    }

    return diff / count;

}

/**
 * Returns standard deviation of 1xN vector
 * @param {Matrix} _data_vector - 1xN vector
 */
lfunction math.std(_data_vector) {

    // std = sqrt(mean(abs(x - x.mean())**2))

    count = Rows(_data_vector);
    diff = 0;

    if (count == 1) {
        _data_vector = Transpose(_data_vector);
        count = Rows(_data_vector);
    }

    mean = math.mean(_data_vector);

    for (_k = 0; _k < count; _k = _k + 1) {
        diff += Abs(_data_vector[_k] - mean) ^ 2;
    }

    return Sqrt(diff / count);

};

/**
 * Returns suite of statistical information of a 1xN vector
 * @param {Matrix} _data_vector - 1xN vector
 */
lfunction math.gather_descriptive_stats(_data_vector) {

    _count = Rows(_data_vector);

    if (count == 1) {
        _data_vector = Transpose(_data_vector);
        _count = Rows(_data_vector);
    }

    _sorted_data_vector = _data_vector % 0;

    _variance = math.variance(_data_vector);
    _std = Sqrt(_variance);
    _sum = math.sum(_data_vector);

    _nonnegatives = utility.filter(_data_vector, "_value_", "_value_ >= 0");

    _dstats = {};
    _dstats["Count"] = _count;
    _dstats["Mean"] = math.mean(_data_vector);
    _dstats["Median"] = math.median(_data_vector);
    _dstats["Min"] = _sorted_data_vector[0];
    _dstats["Max"] = _sorted_data_vector[_count - 1];
    _dstats["2.5%"] = _sorted_data_vector[(_count * 0.025) $1];
    _dstats["97.5%"] = _sorted_data_vector[Min((_count * 0.975 + 0.5) $1, _count - 1)];
    _dstats["Sum"] = _sum;
    _dstats["Mean"] = _sum / _count;
    _dstats["Std.Dev"] = _std;
    _dstats["Variance"] = _variance;
    _dstats["COV"] = _std * _count / _sum;
    _dstats["Skewness"] = math.skewness(_data_vector);
    _dstats["Kurtosis"] = math.kurtosis(_data_vector);
    _dstats["Sq. sum"] = math.sum(_data_vector * 2);
    _dstats["Non-negative"] = Abs(_nonnegatives);

    return _dstats;

}