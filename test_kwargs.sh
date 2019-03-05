#!/bin/bash
# Run the commands with different key word arguments and ensure that the batch file doesn't error uppon startup
# by looking for a known output string in the output file


#*******************************************************************************
#   Global variables
#*******************************************************************************
# File type combinations
declare -a fileTypes=(  "--alignment ./tests/hbltests/data/CD2.nex" # nexus with tree
                        "--alignment ./tests/hbltests/data/CD2_noTree.nex --tree ./tests/hbltests/data/CD2.newick" # nexus with tree separate
                        "--alignment ./tests/hbltests/data/CD2_reduced.fna" # fasta with tree
                        "--alignment ./tests/hbltests/data/CD2_reduced.fasta --tree ./tests/hbltests/data/CD2_reduced.nhx" # fasta with tree separate
                        "--alignment ./tests/hbltests/data/CD2.phylip --tree ./tests/hbltests/data/CD2.newick" # phylip
                        )

declare -a universalArgs=(  "--code Universal"
                            "--branches All"
                            "--branches Internal"
                            "--output ./testMethodOutput"
                            )


#*******************************************************************************
#   General Function
#*******************************************************************************
# Define the general function for evaluating a method
# the function "evaluateMethod" takes three args: 
#   1. the method name (i.e. busted)
#   2. the text to look for in the output
#   3. the array of method specific arguments to test

evaluateMethod() {
    
    methodName=$1
    expectedOutput=$2
    shift # shift function args
    shift # shift function args again so only the methodSpecific arguments are left
    methodSpecificArgs=("$@")

    echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "Testing $methodName"
    echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

    for fileType in "${fileTypes[@]}"
    do
        # Run the command
        hyphyCall="./HYPHYMP $methodName $fileType"
        echo "      ***************************************************************************************************"
        echo "      calling: $hyphyCall"
        eval $hyphyCall > tempOutput.txt &
        sleep 5s
        killall HYPHYMP

        # Confirm that the expected output is there
        if grep -q "$expectedOutput" ./tempOutput.txt; then
            echo "      +++ PASSED +++"
        else
            echo "      --- FAILED ---"
            cat ./tempOutput.txt
            exit 1
        fi    
    done

    for arguments in "${methodSpecificArgs[@]}"
    do
        # Run the command
        hyphyCall="./HYPHYMP $methodName ${fileTypes[0]} $arguments"
        echo "      ***************************************************************************************************"
        echo "      calling: $hyphyCall"
        eval $hyphyCall > tempOutput.txt &
        sleep 5s
        killall HYPHYMP

        # Confirm that the expected output is there
        if grep -q "$expectedOutput" ./tempOutput.txt; then
            echo "      +++ PASSED +++"
        else
            echo "      --- FAILED ---"
            cat ./tempOutput.txt
            exit 1
        fi    
    done
}

#*******************************************************************************
#   Individual Methods
#*******************************************************************************

# BUSTED
declare -a bustedArgs=( "${universalArgs[@]}"
                        "--srv Yes"
                        "--srv No"
                        )                     
evaluateMethod "busted" "Obtaining branch lengths" "${bustedArgs[@]}"

# MEME
declare -a memeArgs=( "${universalArgs[@]}"
                        "--pvalue 0.1"
                        "--pvalue 0.4"
                        )                     
evaluateMethod "meme" "Obtaining branch lengths" "${memeArgs[@]}"