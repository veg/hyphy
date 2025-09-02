#!/bin/bash
# A simple script to run all unit tests against the currently checked out version


#export ASAN_OPTIONS=detect_leaks=1

currentDir=$PWD
export ASAN_OPTIONS=detect_leaks=1
#HYPHYMP=$PWD/HYPHYDEBUG
HYPHYMP=$PWD/hyphy

testsRun=0
testsFailed=0
failedTests=()

for filename in ./tests/hbltests/UnitTests/HBLCommands/*.bf; do
  echo $filename

  # Run the test checking to see if it faield
  $HYPHYMP $filename 2>&1 >/dev/null
  if [ $? -ne 0 ]; then
    ((testFailed++))
    failedTests+=($filename)
  fi

  ((testsRun++))

done

if [ $testFailed ]
then
  echo "\n\n------------------------SUMMARY (Failed Tests)----------------------------\n"
  echo "\n The following tests failed:"

  for failedTest in "${failedTests[@]}"; do
    echo $failedTest
  done

  echo "\n The output of the failed tests is below: \n"
  for failedTest in "${failedTests[@]}"; do
    echo "--------------------------------------------------------------"
    $HYPHYMP $failedTest
    echo "\n"
  done

  echo "\n\n------------------------SUMMARY (Failed Tests)----------------------------\n"
  echo "$testFailed Tests Failed"
  echo "of $testsRun Tests Run"
  exit 1
else
  echo "\n\n------------------------ALL TEST PASSED------------------------------------\n"
fi
