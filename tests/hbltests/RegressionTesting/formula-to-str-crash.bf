function _profileFit(_xxv_,_variableIndex) {
    SetParameter(lf,_variableIndex,_xxv_);
    LFCompute(lf,_xxres);
    return _xxres;
}

x := _profileFit (_xx,0) - 0.05;

GetString (s, x, -2);

fprintf (stdout, s, "\n");
