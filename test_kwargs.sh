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
        sleep 2s
        killall HYPHYMP

        # Confirm that the expected output is there
        if grep -q "$expectedOutput" ./tempOutput.txt; then
            echo "      +++ PASSED +++"
        else
            echo "      --- FAILED ---"
            cat ./tempOutput.txt
            exit 1
        fi
        removeTempFiles
    done

    for arguments in "${methodSpecificArgs[@]}"
    do
        # Run the command
        hyphyCall="./HYPHYMP $methodName ${fileTypes[0]} $arguments"
        echo "      ***************************************************************************************************"
        echo "      calling: $hyphyCall"
        eval $hyphyCall > tempOutput.txt &
        sleep 2s
        killall HYPHYMP

        # Confirm that the expected output is there
        if grep -q "$expectedOutput" ./tempOutput.txt; then
            echo "      +++ PASSED +++"
        else
            echo "      --- FAILED ---"
            cat ./tempOutput.txt
            exit 1
        fi
        removeTempFiles
    done
}

removeTempFiles() {
    if [ -f ./res/TemplateBatchFiles/SelectionAnalyses/tempCache.txt ]; then
        rm ./res/TemplateBatchFiles/SelectionAnalyses/tempCache.txt
    fi 
    if [ -f ./tests/hbltests/data/CD2_AA.fna.FADE.cache ]; then
        rm ./tests/hbltests/data/CD2_AA.fna.FADE.cache
    fi
    sleep .5s
}

#*******************************************************************************
#  Nucleotide Methods
#*******************************************************************************

# aBSREL
declare -a absrelArgs=( "${universalArgs[@]}"
                        )                     
evaluateMethod "absrel" "Obtaining branch lengths" "${absrelArgs[@]}"

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

# SLAC
declare -a slacArgs=( "${universalArgs[@]}"
                        "--samples 0"
                        "--samples 100"
                        "--samples 100000"
                        "--pvalue 0.1"
                        "--pvalue 0.4"
                        )                     
evaluateMethod "slac" "Obtaining branch lengths" "${slacArgs[@]}"


#*******************************************************************************
#  Amino Acid Methods
#*******************************************************************************
# Amino acid methods require different files and args
# Redefining the two array variables to be able to use same "evaluateMethod" function defined previously

declare -a fileTypes=(  "--alignment ./tests/hbltests/data/CD2_AA.fna" # amino acid fna with rooted tree
                        )
declare -a universalArgs=(  "--branches All"
                            "--branches Internal"
                            "--output ./testMethodOutput"
                            )

# FADE
declare -a fadeArgs=( "${universalArgs[@]}"
                        "--cache tempCache.txt"
                        "--grid 25"
                        "--model LG"
                        "--method Metropolis-Hastings"
                        "--method Collapsed-Gibbs"
                        "--method Variational-Bayes"
                        "--concentration_parameter 0.5"
                        "--method Variational-Bayes --grid 50 --model WAG --concentration_parameter 0.01"
                        "--method Collapsed-Gibbs --grid 20 --model WAG --concentration_parameter 0.5 --chains 6 --chain-length 1900000 --burn-in 1200000 --samples 125"
                        )                     
evaluateMethod "fade" "Fitting the baseline" "${fadeArgs[@]}"