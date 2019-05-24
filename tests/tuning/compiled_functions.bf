function regular_sum (N) {
    sum = 0;
    for (i = 0; i < N; i+=1) {
        sum += i;
    }
    return sum;
}

cfunction compiled_sum (N) {
    sum = 0;
    for (i = 0; i < N; i+=1) {
        sum += i;
    }
    return sum;
}

#profile START;

N = 1000000;

fprintf (stdout, "Regular = ", regular_sum (N), "\n");
fprintf (stdout, "Compiled sum = ", compiled_sum (N), "\n");

#profile _hyphy_profile_dump;
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
}

