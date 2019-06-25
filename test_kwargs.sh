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
        sleep 8s
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
        sleep 8s
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
    if [ -f ./res/TemplateBatchFiles/SelectionAnalyses/tempCache.cache ]; then
        rm ./res/TemplateBatchFiles/SelectionAnalyses/tempCache.cache
    fi 
    if [ -f ./tests/hbltests/data/*.cache ]; then
        rm ./tests/hbltests/data/*.cache
    fi
    sleep 1s
}

#*******************************************************************************
#  Nucleotide Methods
#*******************************************************************************

# aBSREL
declare -a absrelArgs=( "${universalArgs[@]}"
                        )                     
evaluateMethod "absrel" "Obtaining branch lengths" "${absrelArgs[@]}"

# GM
declare -a bgmArgs=( "${universalArgs[@]}"
                        "--steps 200000"
                        "--burn-in 20000"
                        "--burn-in 0"
                        "--samples 0"
                        "--max-parents 3"
                        "--min-subs 5"
                        )                     
evaluateMethod "bgm" "Performing initial model fit" "${bgmArgs[@]}"

# BUSTED
declare -a bustedArgs=( "${universalArgs[@]}"
                        "--srv Yes"
                        "--srv No"
                        )                     
evaluateMethod "busted" "Obtaining branch lengths" "${bustedArgs[@]}"

# FEL
declare -a felArgs=( "${universalArgs[@]}"
                        "--srv Yes"
                        "--srv No"
                        )                     
evaluateMethod "fel" "Obtaining branch lengths" "${felArgs[@]}"

# FUBAR
declare -a fubarArgs=( "${universalArgs[@]}"
                        "--cache tempCache.cache"
                        "--grid 25"
                        "--model LG"
                        "--method Metropolis-Hastings"
                        "--method Collapsed-Gibbs"
                        "--method Variational-Bayes"
                        "--concentration_parameter 0.5"
                        "--method Variational-Bayes --grid 50 --model WAG --concentration_parameter 0.01"
                        "--method Collapsed-Gibbs --grid 20 --model WAG --concentration_parameter 0.5 --chains 6 --chain-length 1900000 --burn-in 1200000 --samples 120"
                        )                    
# Disabling the FUBAR and FADE tests for now as they only fail on travis
#evaluateMethod "fubar" "Obtaining branch lengths" "${fubarArgs[@]}"

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
#  GARD
#*******************************************************************************
# Redefining the two array variables to be able to use same "evaluateMethod" function defined previously

declare -a fileTypes=(  "--alignment ./tests/hbltests/data/CD2.nex" # nexus with tree
                        "--alignment ./tests/hbltests/data/CD2_noTree.nex " # nexus with no tree
                        "--alignment ./tests/hbltests/data/CD2.phylip" # phylip
                        )

# GARD
declare -a gardArgs=(   ""
                        "--type nucleotide"
                        "--type codon"
                        "--type codon --code Universal"
                        "--type codon --model JTT"
                        "--rv GDD"
                        "--rv GDD --rate-classes 2"
                        "--output ./testMethodOutput"
                        "--output-lf ./testMethodOutput"
                        "--method Collapsed-Gibbs --grid 20 --model WAG --concentration_parameter 0.5 --chains 6 --chain-length 1900000 --burn-in 1200000 --samples 125"
                        )                     

evaluateMethod "gard" "Fitting the baseline" "${gardArgs[@]}"


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
                        "--cache tempCache.cache"
                        "--grid 25"
                        "--model LG"
                        "--method Metropolis-Hastings"
                        "--method Collapsed-Gibbs"
                        "--method Variational-Bayes"
                        "--concentration_parameter 0.5"
                        "--method Variational-Bayes --grid 50 --model WAG --concentration_parameter 0.01"
                        "--method Collapsed-Gibbs --grid 20 --model WAG --concentration_parameter 0.5 --chains 6 --chain-length 1900000 --burn-in 1200000 --samples 125"
                        )                     
# Disabling the FUBAR and FADE tests for now as they only fail on travis
#evaluateMethod "fade" "Fitting the baseline" "${fadeArgs[@]}"