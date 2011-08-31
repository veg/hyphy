/*

A collection of HBL utility functions shared by various test 
scripts in the HyPhy distribution. 

Written by SL Kosakovsky Pond (spond@ucsd.edu), March 2009

For linencsing term see
http://www.gnu.org/licenses/quick-guide-gplv3.html

*/

RequireVersion ("2");

/*###########################################################################

Call startTestTimer to set up timing variables and (optionally) print out 
the description of a test (passed as an argument). Nothing is printed if
an something other than a string is passed (e.g. a 0)

----------------------------------------------------------------------------*/

function startTestTimer (descriptiveString)
{
	_hyphyTestTimerInWall = Time (1);
	_hyphyTestTimerInUser = Time (0);
	
	if (Type(descriptiveString) == "String")
	{
		fprintf (stdout, "[STARTING TEST  : ", descriptiveString, "]\n");
	}
	
	return 0;
}

/*###########################################################################

Call endTestTimer following a call to startTestTimer to finish timing a test 
and optionally print out a diagnostic message, run time (both wall clock and
user space) and their ratio (cpu utilization). Nothing is printed if
an something other than a string is passed (e.g. a 0). Returns a 1x2 matrix
of wall clock time and user space time elapsed (in seconds).

----------------------------------------------------------------------------*/

function endTestTimer (descriptiveString)
{
	_hyphyTestRangeWall = Time(1) - _hyphyTestTimerInWall;
	_hyphyTestRangeUser = Time(0) - _hyphyTestTimerInUser;
	
	if (Type(descriptiveString) == "String")
	{
		fprintf (stdout, "[FINISHED TEST  : ", descriptiveString, "]\n");
		fprintf (stdout, "[WALL CLOCK TIME: ", reportTimeUsed (_hyphyTestRangeWall), "]\n");
		fprintf (stdout, "[CPU  TIME      : ", reportTimeUsed (_hyphyTestRangeUser), "]\n");
		fprintf (stdout, "[CPU UTILIZATION: ", Format(_hyphyTestRangeUser/_hyphyTestRangeWall*100,5,2), "%]\n");
	}
	
	return {{_hyphyTestRangeWall__, _hyphyTestRangeUser__}};
}

/*###########################################################################

reportTimeUsed takes a positive numeric argument (seconds) and returns it as a string
formatted as hh:mm:ss 
An error string is returned if a negative value is passed in.

----------------------------------------------------------------------------*/

function reportTimeUsed (secondCount)
{
	if (secondCount >= 0)
	{
		_hrs  = secondCount$3600;
		_mins = (secondCount%3600)$60;
		_secs = secondCount - secondCount$60*60;
		
		if (_hrs < 10)
		{
			_hrs = "0" + _hrs;
		}
		else
		{
			_hrs = "" + _hrs;
		}
		
		if (_mins < 10)
		{
			_mins = "0" + _mins;
		}
		else
		{
			_mins = "" + _mins;
		}

		if (_secs < 10)
		{
			_secs = "0" + _secs;
		}
		else
		{
			_secs = "" + _secs;
		}
		return _hrs + ":" + _mins + ":" + _secs;

	}
	return "[ERROR: NEGATIVE TIME ELAPSED VALUE]\n";
}

/*###########################################################################

logTestResults takes a boolean (numeric argument); reports test results
and logs them to a test database (if one is open)

----------------------------------------------------------------------------*/

function logTestResult (testResult)
{
	if (testResult)
	{
		fprintf (stdout, "[TEST PASSED]\n");
	}
	else
	{
		fprintf (stdout, "[TEST FAILED]\n");
	}
	return testResult != 0; 
}

