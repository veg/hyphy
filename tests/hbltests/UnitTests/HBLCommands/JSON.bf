// Tests associativeListtoJSON function.
LoadFunctionLibrary("lib2014/UtilityFunctions.bf");

test_json    =    {"fits"           : {},
                  "timers"          : {},
                  "profiles"        : "[[inf,0.9779047623694379],[inf,0.02209523763056209]]","GOOSE_HONGKONG_W355_97" : "[[1,-nan]]","DUCK_HONGKONG_Y283_97" : "[[-nan,inf]]",
                  "evidence ratios" : {}
                  };


new_test_json = utility.associativeListToJSON(test_json);

new_test_string = " " + new_test_json + " ";

// inf and -nan are not parseable JSON
inf_exists = new_test_string/"*inf*";
nan_exists = new_test_string/"*-nan*";
e999_exists = new_test_string/"*1e999*";
null_exists = new_test_string/"*null*";
assert(!inf_exists, "inf still exists in JSON");
assert(!nan_exists, "-nan still exists in JSON");
assert(e999_exists, "1e999 does not exist in JSON");
assert(null_exists, "null does not exist in JSON");
