//LoadFunctionLibrary ("libv3/UtilityFunctions.bf");


//#profile START;

array_in = {};
N = 500000;
M = 100;
for (k = 0; k < N; k+=1) {
    array_in + k;
}

array_out = {};

/*for (m = 0; m < M; m+=1) {
    ForEachPair (array_in, '_k_', '_v_', 'array_out[_k_]=_v_^2');
}*/

/*#profile _hyphy_profile_dump;
stats  			= _hyphy_profile_dump["STATS"];
_profile_summer = ({1,Rows(stats)}["1"]) * stats;
_instructions   = _hyphy_profile_dump["INSTRUCTION"];
_indices	    = _hyphy_profile_dump["INSTRUCTION INDEX"];

fprintf (stdout, "\nTotal run time (seconds)      : ", Format(_profile_summer[1],15,6),
                 "\nTotal number of steps         : ", Format(_profile_summer[0],15,0), "\n\n");

to_sort        =  stats["-_MATRIX_ELEMENT_VALUE_*_MATRIX_ELEMENT_COLUMN_+(_MATRIX_ELEMENT_COLUMN_==0)*_MATRIX_ELEMENT_ROW_"] % 1;

for (k=0; k<Columns(_instructions); k=k+1)
{
    k2 = to_sort[k][0];
    fprintf (stdout, Format (_indices[k2],6,0), " : ", _instructions[k2], "\n\tCall count: ", stats[k2][0],
                                                   "\n\tTime (seconds): ", stats[k2][1], "\n");
}*/


function ForEachPair(object, key_name, value_name, transform) {

    io.CheckAssertion ("!utility.ForEachPair.warn_non_reentrant", "utility.ForEachPair is non re-entrant");
    utility.ForEachPair.warn_non_reentrant = TRUE;

    ExecuteCommands ("function  utility.ForEachPair.CB (`key_name`, `value_name`) {`transform`}", enclosing_namespace);

    if (Type (object) == "AssociativeList") {
        utility.ForEachPair.keys = Rows (object);
        utility.ForEachPair.ub =  Abs (object);

         ExecuteCommands ('
         for (utility.ForEachPair.k = 0; utility.ForEachPair.k < utility.ForEachPair.ub; utility.ForEachPair.k += 1) {
            utility.ForEachPair.key = utility.ForEachPair.keys[utility.ForEachPair.k];
            Call ("utility.ForEachPair.CB",utility.ForEachPair.key , object [utility.ForEachPair.key]);
        }', compiled);
    } else {
        if (Type (object) == "Matrix") {
            utility.ForEachPair.rows = Rows (object);
            utility.ForEachPair.columns = Columns (object);

            if (utility.ForEachPair.rows && utility.ForEachPair.columns) {

                utility.ForEachPair.key = {{utility.ForEachPair.r,utility.ForEachPair.c}};
                for (utility.ForEachPair.r = 0; utility.ForEachPair.r < utility.ForEachPair.rows; utility.ForEachPair.r += 1) {
                    for (utility.ForEachPair.c = 0; utility.ForEachPair.c < utility.ForEachPair.columns; utility.ForEachPair.c += 1) {
                        Call ("utility.ForEachPair.CB",utility.ForEachPair.key , object [utility.ForEachPair.r][utility.ForEachPair.c]);
                    }
                }
            }
        }
    }

    Eval ("`key_name` = None");
    Eval ("`value_name` = None");
    // reset bindings here to avoid calls on stale formula in subsequent invocations

    utility.ForEachPair.warn_non_reentrant = FALSE;

}
