# Running MSS GA

The MSS Genetic Algorithm (GA) estimates an MSS model by dividing synonymous codon pairs into either NEUTRAL or SELECTED classes. 

#### Requirements

1.  The `develop` branch of [HyPhy](https://github.com/veg/hyphy)

##### What you need

1. Alignment files you wish to run the GA on. Each alignment file must include a tree
2. A line list of relative paths to these files

### Example

#### Make line list of file names

```
ls -1 ./files/* > files.lst
```

#### Run by codon GA using the BIC selection criterion

```
mpirun -np 6 HYPHYMPI MSS-GA \
--ic BIC \
--filelist files.lst \
--output BIC.json
```

Partial output...

* Loaded a list with **40** files
code: Universal
ic: BIC
mss-type: SynREVCodon
mss-reference-rate: synonymous

### Individual file statistics and simple model fits

|                                                        Filepath                                                        |   Sequences    |     Sites      |   TreeLength   |     Log(L)     |
|------------------------------------------------------------------------------------------------------------------------|----------------|----------------|----------------|----------------|
|          /Users/hverdonk/./data/FBgn0000014_NT_033777.3_12_species_12_seqs.msa.raxml.bestTree.replicate.1.nex          |       12       |      500       |       0.280    |    -4497.3412  |
|           /Users/hverdonk/./data/FBgn0000015_NT_033777.3_4_species_4_seqs.msa.raxml.bestTree.replicate.1.nex           |       4        |      500       |       0.039    |    -2405.3639  |
|           /Users/hverdonk/./data/FBgn0000017_NT_037436.4_5_species_5_seqs.msa.raxml.bestTree.replicate.1.nex           |       5        |      500       |       0.187    |    -3212.6213  |
|          /Users/hverdonk/./data/FBgn0000018_NT_033779.5_11_species_11_seqs.msa.raxml.bestTree.replicate.1.nex          |       11       |      500       |       0.616    |    -6370.0707  |
|          /Users/hverdonk/./data/FBgn0000100_NT_037436.4_12_species_13_seqs.msa.raxml.bestTree.replicate.1.nex          |       12       |      500       |       0.171    |    -3628.1248  |

....

```
[MSS-GA]   2 rates [ generation      1, total models       65,    0% converged,  0.16/sec]. Min (c-AIC) =  392084.2262 [no improvement], rate profile 0000000000111001100001010010101111101100111100110001011001101000111       0000000001111000100011011110111110101000111100110001111101101000110 IC =     392096.8403, delta =   -20.06
        0000001010110011100001011110101111001100101000110001001010101001110 IC =     392087.8570, delta =   -11.07
        0000100000111001000001100111110100101100100101010101101100110000000 IC =     392095.9454, delta =   -19.16
        0001010100001000100001010011101111100111111000100011011001001000111 IC =     392092.3239, delta =   -15.54
        0001110010110000010001111110100110011100100101110111101000111111010 IC =     392097.8549, delta =   -21.07
        0010000100101001101100111111100000001110111001010001111000111000110 IC =     392100.8079, delta =   -24.02
        0010000101100000000001011011001110101100110010110000101101100110101 IC =     392101.1409, delta =   -24.36
        0010000101100011000001011011010000001000110001010001001101100000010 IC =     392101.1555, delta =   -24.37
        0010010000111111010000110010101110110101111100110101111001011100111 IC =     392096.2922, delta =   -19.51
        0010010110000101100000000011111111010101110010000101010011111100111 IC =     392087.8727, delta =   -11.09
        0010110000101010010000110010101111100100111111110111111010011101011 IC =     392099.5169, delta =   -22.73
        0010110100001111000010001011100111011101111001100011101011111000110 IC =     392088.1948, delta =   -11.41
        0010111110001101010001111011110110110100111110111011111010011011010 IC =     392101.1305, delta =   -24.35
        0011000001101000110001010010101111100110111001110001011100001110010 IC =     392100.7144, delta =   -23.93
        0011001011100010001110011111111101001110111011010101001100001110010 IC =     392091.9536, delta =   -15.17
```

#### Process the results

```
hyphy MSS-GA-processor --json BIC.json --tsv BIC.tsv
```

#### Repeat the run with c-AIC 

```
mpirun -np 6 HYPHYMPI MSS-GA \
--ic AIC-c \
--filelist files.lst \
--output AIC.json
```

#### Process the results

```
hyphy MSS-GA-processor --json AIC.json --tsv AIC.tsv
```

#### Compare simulated truth with results



