LoadFunctionLibrary("libv3/IOFunctions.bf");

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

function test_reportStatsMD () {

    io.ReportStatsMD("aBSREL", absrel.stats);
}


// TODO: Doesn't work due to not being able to parse arrays
function test_parse_json() {
    file_path = "./data/CD2.nex.slac.json";
    parsed_json = io.ParseJSON(file_path);
    fprintf(stdout, parsed_json);
}

function test_check_key() {
    assert(utility.KeyExists(absrel.stats, "Count"), "Count not found");
}

expected_loglikelihood = -3467.072352344857; 
actual = -3467.082352344857;

fprintf(stdout,  Abs(actual - expected_loglikelihood));
assert(Abs(actual - expected_loglikelihood) <= 1, "whoops");

test_reportStatsMD();
test_check_key();
