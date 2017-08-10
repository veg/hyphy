LoadFunctionLibrary("libv3/UtilityFunctions.bf");

absrel.stats = {
 "Count":4,
 "Mean":2.25,
 "Median":2,
 "Min":1,
 "Max":4,
 "2.5%":1,
 "97.5%":4,
 "Sum":9,
 "Std.Dev":1.089724735885168,
 "Variance":1.1875,
 "COV":0.4843221048378527,
 "Skewness":0.6520236646847543,
 "Kurtosis":2.096952908587257,
 "Sq. sum":18,
 "Non-negative":4
};

function test_utility_for_each_pair() {
    test = {};
    utility.ForEachPair(absrel.stats, "_key_", "_value_",
    "test[_key_] = _value_;");
    assert(test["Count"] == 4, "dictionary did not iterate appropriately");
}

function test_utility.Keys() {
    test = {};
    keys = utility.Keys(absrel.stats);
    expected_keys =  {terms.math.count, terms.math.mean, terms.math.median, terms.math.min, terms.math.max, terms.math._2.5, terms.mat._97.5,
    terms.math.sum, terms.math.stddev, terms.math.variance, terms.math.cov, terms.math.skewness, terms.math.kurtosis, terms.math.square_sum,
    terms.math.non_negative};
    //TODO: write function to test string matrix equality
    assert(Abs(keys) == 10, "did not return all keys");
}

function test_set_intersection() {

    set1 = {};
    first = {{"a","b","c"}};
    set2 = {};
    second = {{"a","b","c","d"}};


    utility.AddToSet(set1, first);
    utility.AddToSet(set2, second);

    set3 = {};
    utility.Intersect(set3, set1, set2);

    //expected_set = {
    //    "a" : 1,
    //    "b" : 1,
    //    "c" : 1
    //};
    assert(set3["a"] == 1);
    assert(set3["b"] == 1);
    assert(set3["c"] == 1);
    assert(set3["d"] != 1);

    utility.Intersect(set3, set2, set1);
    assert(set3["a"] == 1);
    assert(set3["b"] == 1);
    assert(set3["c"] == 1);
    assert(set3["d"] != 1);

}

test_utility_for_each_pair();
test_utility.Keys();
test_set_intersection();


