ExecuteAFile (PATH_TO_CURRENT_BF + "TestTools.ibf");
runATest ();


function getTestName () {
	return "Differentiate";
}

function getTestedFunctions () {
	return {{"_ElementaryCommand::HandleDifferentiate"}};
}

lfunction checkExpression (expression2check, expected_value, errorMessage) {
	result = 0;
	GetString (expr, ^expression2check, -2);
	assert (expr == expected_value, errorMessage + " (wrong result = " + expr + ")");
	result = 1;
	return result;
}


function runTest () {
	ASSERTION_BEHAVIOR = 1; /* print warning to console and go to the end of the execution list */
	testResult  	   = 0;

 
    assert (runCommandWithSoftErrors ("Differentiate (2 invalid, x^2,x)", "is not a valid variable identifier in call to Differentiate"), "Invalid variable identifier in call to Differentiate.");
    assert (runCommandWithSoftErrors ("Differentiate (dfx, x^2,x,\"beavis\")", "\\(the number of times to differentiate\\) must be a non-negative integer"), "Invalid order option in call to Differentiate.");
    assert (runCommandWithSoftErrors ("Differentiate (dfx, x^2,x-y)", "is not a valid variable identifier"), "Invalid order option in call to Differentiate.");
    assert (runCommandWithSoftErrors ("Differentiate (dfx, x^2,x,z#4)", "<ERROR HERE>"), "Unparseable order option in call to Differentiate.");
    assert (runCommandWithSoftErrors ("Differentiate (dfx, Time(x),x)","Differentiation of .+ failed"), "Unparseable order option in call to Differentiate.");


    Differentiate (dlog, Log (x*y),x);    
    if (!checkExpression ("dlog", "x^(-1)","Failed to correctly differentiate and simplify Log (x*y)")) {
    	return 0;
    }
    
    Differentiate (mlog, -y + x*y + y*z*w + 5*y, y);
    if (!checkExpression ("mlog", "4+z*w+x","Failed to correctly differentiate (wrt y) and simplify -y + x*y + y*z + 5*y")) {
    	return 0;
    }
    

    Differentiate (res, (x*y)^2, x);
    if (!checkExpression ("res", "2*x*y^2","Failed to correctly differentiate (wrt x) and simplify (x*y)^2")) {
    	return 0;
    }

    Differentiate (res, Exp (Log (2*x)), x); //  2
    if (!checkExpression ("res", "2","Failed to correctly differentiate and simplify Exp (Log (2*x))")) {
    	return 0;
    }


    Differentiate (res, x^3,x,2);
    if (!checkExpression ("res", "6*x","Failed to correctly 2x differentiate and simplify x^3")) {
    	return 0;
    }
   //assert (Abs(dfx+12)<1e-10, "Checking (x^3)'' == 6x derivative");

    Differentiate (res, 0.00154446*t*XbpPYxBh.user.theta_WY+0.0074036*t*XbpPYxBh.user.theta_VY+0.00263518*t*XbpPYxBh.user.theta_VW+0.008335179999999999*t*XbpPYxBh.user.theta_TY+0.00296676*t*XbpPYxBh.user.theta_TW+0.0142216*t*XbpPYxBh.user.theta_TV+0.00367729*t*XbpPYxBh.user.theta_SY+0.00130886*t*XbpPYxBh.user.theta_SW+0.00627424*t*XbpPYxBh.user.theta_SV+0.00706371*t*XbpPYxBh.user.theta_ST+0.00375083*t*XbpPYxBh.user.theta_RY+0.00133504*t*XbpPYxBh.user.theta_RW+0.00639972*t*XbpPYxBh.user.theta_RV+0.00720499*t*XbpPYxBh.user.theta_RT+0.00317867*t*XbpPYxBh.user.theta_RS+0.00387341*t*XbpPYxBh.user.theta_QY+0.00137867*t*XbpPYxBh.user.theta_QW+0.00660886*t*XbpPYxBh.user.theta_QV+0.00744044*t*XbpPYxBh.user.theta_QT+0.00328255*t*XbpPYxBh.user.theta_QS+0.0033482*t*XbpPYxBh.user.theta_QR+0.00291731*t*XbpPYxBh.user.theta_PY+0.00103837*t*XbpPYxBh.user.theta_PW+0.00497756*t*XbpPYxBh.user.theta_PV+0.00560388*t*XbpPYxBh.user.theta_PT+0.0024723*t*XbpPYxBh.user.theta_PS+0.00252175*t*XbpPYxBh.user.theta_PR+0.00260416*t*XbpPYxBh.user.theta_PQ+0.00492756*t*XbpPYxBh.user.theta_NY+0.00175388*t*XbpPYxBh.user.theta_NW+0.00840748*t*XbpPYxBh.user.theta_NV+0.009465370000000001*t*XbpPYxBh.user.theta_NT+0.0041759*t*XbpPYxBh.user.theta_NS+0.00425942*t*XbpPYxBh.user.theta_NR+0.00439861*t*XbpPYxBh.user.theta_NQ+0.00331288*t*XbpPYxBh.user.theta_NP+0.00181413*t*XbpPYxBh.user.theta_MY+0.000645706*t*XbpPYxBh.user.theta_MW+0.00309529*t*XbpPYxBh.user.theta_MV+0.00348476*t*XbpPYxBh.user.theta_MT+0.0015374*t*XbpPYxBh.user.theta_MS+0.00156814*t*XbpPYxBh.user.theta_MR+0.00161939*t*XbpPYxBh.user.theta_MQ+0.00121967*t*XbpPYxBh.user.theta_MP+0.00206011*t*XbpPYxBh.user.theta_MN+0.00558947*t*XbpPYxBh.user.theta_LY+0.00198947*t*XbpPYxBh.user.theta_LW+0.00953684*t*XbpPYxBh.user.theta_LV+0.0107368*t*XbpPYxBh.user.theta_LT+0.00473684*t*XbpPYxBh.user.theta_LS+0.00483158*t*XbpPYxBh.user.theta_LR+0.00498947*t*XbpPYxBh.user.theta_LQ+0.00375789*t*XbpPYxBh.user.theta_LP+0.00634737*t*XbpPYxBh.user.theta_LN+0.00233684*t*XbpPYxBh.user.theta_LM+0.00637396*t*XbpPYxBh.user.theta_KY+0.0022687*t*XbpPYxBh.user.theta_KW+0.0108753*t*XbpPYxBh.user.theta_KV+0.0122438*t*XbpPYxBh.user.theta_KT+0.00540166*t*XbpPYxBh.user.theta_KS+0.0055097*t*XbpPYxBh.user.theta_KR+0.00568975*t*XbpPYxBh.user.theta_KQ+0.00428532*t*XbpPYxBh.user.theta_KP+0.00723823*t*XbpPYxBh.user.theta_KN+0.00266482*t*XbpPYxBh.user.theta_KM+0.008210530000000001*t*XbpPYxBh.user.theta_KL+0.0070849*t*XbpPYxBh.user.theta_IY+0.00252175*t*XbpPYxBh.user.theta_IW+0.0120884*t*XbpPYxBh.user.theta_IV+0.0136094*t*XbpPYxBh.user.theta_IT+0.00600416*t*XbpPYxBh.user.theta_IS+0.00612424*t*XbpPYxBh.user.theta_IR+0.00632438*t*XbpPYxBh.user.theta_IQ+0.0047633*t*XbpPYxBh.user.theta_IP+0.00804557*t*XbpPYxBh.user.theta_IN+0.00296205*t*XbpPYxBh.user.theta_IM+0.00912632*t*XbpPYxBh.user.theta_IL+0.0104072*t*XbpPYxBh.user.theta_IK+0.00188767*t*XbpPYxBh.user.theta_HY+0.000671884*t*XbpPYxBh.user.theta_HW+0.00322078*t*XbpPYxBh.user.theta_HV+0.00362604*t*XbpPYxBh.user.theta_HT+0.00159972*t*XbpPYxBh.user.theta_HS+0.00163172*t*XbpPYxBh.user.theta_HR+0.00168504*t*XbpPYxBh.user.theta_HQ+0.00126911*t*XbpPYxBh.user.theta_HP+0.00214363*t*XbpPYxBh.user.theta_HN+0.0007891970000000001*t*XbpPYxBh.user.theta_HM+0.00243158*t*XbpPYxBh.user.theta_HL+0.00277285*t*XbpPYxBh.user.theta_HK+0.00308213*t*XbpPYxBh.user.theta_HI+0.00661911*t*XbpPYxBh.user.theta_GY+0.00235596*t*XbpPYxBh.user.theta_GW+0.0112936*t*XbpPYxBh.user.theta_GV+0.0127147*t*XbpPYxBh.user.theta_GT+0.00560942*t*XbpPYxBh.user.theta_GS+0.00572161*t*XbpPYxBh.user.theta_GR+0.00590859*t*XbpPYxBh.user.theta_GQ+0.00445014*t*XbpPYxBh.user.theta_GP+0.00751662*t*XbpPYxBh.user.theta_GN+0.00276731*t*XbpPYxBh.user.theta_GM+0.00852632*t*XbpPYxBh.user.theta_GL+0.009722990000000001*t*XbpPYxBh.user.theta_GK+0.0108075*t*XbpPYxBh.user.theta_GI+0.0028795*t*XbpPYxBh.user.theta_GH+0.00406953*t*XbpPYxBh.user.theta_FY+0.00144848*t*XbpPYxBh.user.theta_FW+0.00694349*t*XbpPYxBh.user.theta_FV+0.00781717*t*XbpPYxBh.user.theta_FT+0.00344875*t*XbpPYxBh.user.theta_FS+0.00351773*t*XbpPYxBh.user.theta_FR+0.00363269*t*XbpPYxBh.user.theta_FQ+0.00273601*t*XbpPYxBh.user.theta_FP+0.00462133*t*XbpPYxBh.user.theta_FN+0.00170139*t*XbpPYxBh.user.theta_FM+0.00524211*t*XbpPYxBh.user.theta_FL+0.00597784*t*XbpPYxBh.user.theta_FK+0.0066446*t*XbpPYxBh.user.theta_FI+0.00177036*t*XbpPYxBh.user.theta_FH+0.00620776*t*XbpPYxBh.user.theta_FG+0.00436371*t*XbpPYxBh.user.theta_EY+0.00155319*t*XbpPYxBh.user.theta_EW+0.00744543*t*XbpPYxBh.user.theta_EV+0.008382270000000001*t*XbpPYxBh.user.theta_ET+0.00369806*t*XbpPYxBh.user.theta_ES+0.00377202*t*XbpPYxBh.user.theta_ER+0.00389529*t*XbpPYxBh.user.theta_EQ+0.0029338*t*XbpPYxBh.user.theta_EP+0.0049554*t*XbpPYxBh.user.theta_EN+0.00182438*t*XbpPYxBh.user.theta_EM+0.00562105*t*XbpPYxBh.user.theta_EL+0.00640997*t*XbpPYxBh.user.theta_EK+0.00712493*t*XbpPYxBh.user.theta_EI+0.00189834*t*XbpPYxBh.user.theta_EH+0.00665651*t*XbpPYxBh.user.theta_EG+0.00409252*t*XbpPYxBh.user.theta_EF+0.00573657*t*XbpPYxBh.user.theta_DY+0.00204183*t*XbpPYxBh.user.theta_DW+0.009787809999999999*t*XbpPYxBh.user.theta_DV+0.0110194*t*XbpPYxBh.user.theta_DT+0.0048615*t*XbpPYxBh.user.theta_DS+0.00495873*t*XbpPYxBh.user.theta_DR+0.00512078*t*XbpPYxBh.user.theta_DQ+0.00385679*t*XbpPYxBh.user.theta_DP+0.0065144*t*XbpPYxBh.user.theta_DN+0.00239834*t*XbpPYxBh.user.theta_DM+0.00738947*t*XbpPYxBh.user.theta_DL+0.008426589999999999*t*XbpPYxBh.user.theta_DK+0.00936648*t*XbpPYxBh.user.theta_DI+0.00249557*t*XbpPYxBh.user.theta_DH+0.00875069*t*XbpPYxBh.user.theta_DG+0.00538006*t*XbpPYxBh.user.theta_DF+0.00576898*t*XbpPYxBh.user.theta_DE+0.00220637*t*XbpPYxBh.user.theta_CY+0.000785319*t*XbpPYxBh.user.theta_CW+0.00376454*t*XbpPYxBh.user.theta_CV+0.00423823*t*XbpPYxBh.user.theta_CT+0.00186981*t*XbpPYxBh.user.theta_CS+0.0019072*t*XbpPYxBh.user.theta_CR+0.00337812*t*XbpPYxBh.user.theta_AC+0.0087831*t*XbpPYxBh.user.theta_AD+0.00668116*t*XbpPYxBh.user.theta_AE+0.00623075*t*XbpPYxBh.user.theta_AF+0.0101343*t*XbpPYxBh.user.theta_AG+0.00289017*t*XbpPYxBh.user.theta_AH+0.0108475*t*XbpPYxBh.user.theta_AI+0.009759*t*XbpPYxBh.user.theta_AK+0.00855789*t*XbpPYxBh.user.theta_AL+0.00277756*t*XbpPYxBh.user.theta_AM+0.00754446*t*XbpPYxBh.user.theta_AN+0.00446662*t*XbpPYxBh.user.theta_AP+0.00593047*t*XbpPYxBh.user.theta_AQ+0.0057428*t*XbpPYxBh.user.theta_AR+0.00563019*t*XbpPYxBh.user.theta_AS+0.0127618*t*XbpPYxBh.user.theta_AT+0.0113355*t*XbpPYxBh.user.theta_AV+0.00236468*t*XbpPYxBh.user.theta_AW+0.00664363*t*XbpPYxBh.user.theta_AY+0.0029169*t*XbpPYxBh.user.theta_CD+0.00221884*t*XbpPYxBh.user.theta_CE+0.00206925*t*XbpPYxBh.user.theta_CF+0.00336565*t*XbpPYxBh.user.theta_CG+0.000959834*t*XbpPYxBh.user.theta_CH+0.00360249*t*XbpPYxBh.user.theta_CI+0.003241*t*XbpPYxBh.user.theta_CK+0.00284211*t*XbpPYxBh.user.theta_CL+0.000922438*t*XbpPYxBh.user.theta_CM+0.00250554*t*XbpPYxBh.user.theta_CN+0.00148338*t*XbpPYxBh.user.theta_CP+0.00196953*XbpPYxBh.user.theta_CQ*t,XbpPYxBh.user.theta_WY,1);
    if (!checkExpression ("res", "0.00154446*t","Failed to correctly differentiate branch length expression")) {
    	return 0;
    }


	testResult = 1;

	return testResult;
}
