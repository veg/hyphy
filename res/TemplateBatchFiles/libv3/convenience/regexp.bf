
/** @module regexp */


namespace terms {

    namespace regexp {
    
        strings    = "strings";
        separators = "separators";
   
    }

}
 
/**
* @name regexp.Replace
* @param {String} string
* @param {String} re - search for this expression
* @param {String} repl - replace each occurence of re with repl
* @returns {Dictionary} a dictionary containing split strings
*/
 
lfunction regexp.Replace(string, re, repl) {
    return string ^ {{re,repl}};
}

/**
 * @name regexp.Split
 * @param {String} string
 * @param {String} re - the regular expression
 * @returns {Dictionary} a dictionary containing split strings
 */
lfunction regexp.Split(string, re) {
    coordinates = string || re;
    if (coordinates[0] >= 0) {
        result = {};

        current_end = 0;

        for (i = 0; i < Rows(coordinates); i += 2) {
            from = coordinates[i];

            if (current_end < from) {
                result + string[current_end][from - 1];
            } else {
                result + "";
            }
            current_end = coordinates[i + 1] + 1;
        }

        if (current_end < Abs(string)) {
            result + string[current_end][Abs(string) - 1];
        } else {
            result + "";
        }

        return result;
    }
    return {
        "0": string
    };
}

/**
 * @name regexp.SplitWithMatches
 * @param {String} string
 * @param {String} re - the regular expression
 * @returns {Dictionary} a dictionary containing split strings under the key "strings" and separators under the key "separators"
 */
lfunction regexp.SplitWithMatches (string, re) {
    coordinates = string || re;
    if (coordinates[0] >= 0) {
        strings = {};
        separators = {};
        
        current_end = 0;

        for (i = 0; i < Rows(coordinates); i += 2) {
            from = coordinates[i];
            separators + string[coordinates[i]][coordinates[i+1]];

            if (current_end < from) {
                strings + string[current_end][from - 1];
            } else {
                strings + "";
            }
            current_end = coordinates[i + 1] + 1;
        }

        if (current_end < Abs(string)) {
            strings + string[current_end][Abs(string) - 1];
        } else {
            strings + "";
        }

        return {
            ^("terms.regexp.strings") :strings,
            ^("terms.regexp.separators") : separators
        };
    }
    return {
        ^("terms.regexp.strings") : {
            "0": string
        },
        ^("terms.regexp.separators") : {
        }
    };
}

/**
 * @name regexp.Find
 * @param {String} string
 * @param {String} re - the regular expression
 * @returns {String} the portion of the string that is matching
 */
lfunction regexp.Find(string, re) {
    //console.log (string);
    coordinates = string $ re;
    if (coordinates[0] >= 0) {
        return string[coordinates[0]][coordinates[1]];
    }
    return None;
}

/**
 * @name regexp.FindSubexpressions
 * @param {String} string
 * @param {String} re - the regular expression
 * @returns {AssociativeList} all matched RE subsexpressions
 */
lfunction regexp.FindSubexpressions(string, re) {
    coordinates = string $ re;

    if (coordinates[0] >= 0) {
        matched = {};
        upto = utility.Array1D (coordinates);
        for (k = 0; k < upto; k += 2) {
            matched + string[coordinates[k]][coordinates[k+1]];
        }
        return matched;
    }
    return None;
}

/**
 * @name regexp.PartitionByRegularExpressions
 * @param {Dict/Matrix} strings set of strings to partition (if Dict, will use the set of VALUES)
 * @param {Dict/Matrix} rex - string matrix of regular expressions ((if Dict, will use the set of VALUES)
 * @returns {Dict} "reg-exp" : "set of matched strings" (as dict) PLUS a special key ("") which did not match any regular expression
 */
lfunction regexp.PartitionByRegularExpressions(strings, rex) {
    result = {};
    
    for (_value_; in; rex) {
        result[_value_] = {};
    }

    result[""] = {};
    matched_regexp = None;
    
    for (_value_; in; strings) {
        matched_regexp = None;
        for (re; in; rex) {
            matched_regexp = regexp.Find (_value_, re);
            if (matched_regexp != None) {
                break;
            }
        }
        if (None == matched_regexp) {
            result[""] + _value_;
        } else {
            result[re] + _value_;
        }
    }

    
    return result;
}

/**
 * @name regexp.PartitionByRegularExpressionsWithSub
 * partition the set of strings by regular expressions; if a regular expression 
 * contains subexpressions, then these will be added to the matched pattern 
 * @param {Dict/Matrix} strings set of strings to partition (if Dict, will use the set of VALUES)
 * @param {Dict/Matrix} rex - string matrix of regular expressions (if Dict, will use the set of VALUES)
 * @returns {Dict} "reg-exp" : "set of matched strings" (as dict) PLUS a special key ("") which did not match any regular expression
 */
lfunction regexp.PartitionByRegularExpressionsWithSub(strings, rex) {
    result = {};
    result[""] = {};
    matched_regexp = None;
    subs = None;

    for (_value_; in; strings ) {
        matched_regexp = utility.First (rex, "_regex_", "None!=regexp.Find ('" + _value_ +"', _regex_)");
        if (None == matched_regexp) {
            result[""] + _value_;
        } else {
            subs = Join ("|",regexp.FindSubexpressions (_value_, matched_regexp));
            if (utility.Has (result,subs,"AssociativeList") == FALSE) {
                result[subs] = {};
            }
            result[subs] + _value_;        
        }
    }
    
    return result;
}



