function check_value (expected, observed, rel_tolerance) {
    if (expected == observed && observed == 0) {
        return TRUE;
    }
    return Abs (expected - observed) / Abs (Max (observed,expected)) <= rel_tolerance;
}
