# A simple script to run all unit tests against the currently checked out version

currentDir=$PWD
HYPHYMP=$PWD/HYPHYMP

testsRun=0
testsFailed=0
failedTests=()

for filename in ./tests/hbltests/UnitTests/HBLCommands/*.bf; do

  # Run the test checking to see if it faield
  if $HYPHYMP $filename | grep -q "TEST FAILED"; then
    ((testFailed++))
    failedTests+=($filename)
  fi

  ((testsRun++))
  
done

echo "\n\n------------------------SUMMARY (Failed Tests)----------------------------\n"
echo "\n The following tests failed:"

for failedTest in "${failedTests[@]}"; do
  $failedTest
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
