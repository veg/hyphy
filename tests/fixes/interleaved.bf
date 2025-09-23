DataSet ds = ReadDataFile (PATH_TO_CURRENT_BF + "interleaved.nex");

GetDataInfo (di, ds, "SITE_MAP");

fprintf (stdout, di);