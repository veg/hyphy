RequireVersion ("2.5.7");

LoadFunctionLibrary("../codon.bf");
LoadFunctionLibrary("../DNA.bf");
LoadFunctionLibrary("../parameters.bf");
LoadFunctionLibrary("../frequencies.bf");
LoadFunctionLibrary("../../UtilityFunctions.bf");
LoadFunctionLibrary("MG_REV.bf");
LoadFunctionLibrary("../protein.bf");

/** @module models.codon.MG_REV_PROPERTIES */

/**
 * @name models.codon.MG_REV_PROPERTIES.ModelDescription
 * @param {String} type
 * @param {String} code
 * @param {Dict}   properties
 *                 'name' -> dict {
                        'one letter amino-acid' -> property
                    }
 */



terms.model.residue_name_map = "residue_to_id_map";
terms.model.MG_REV_PROPERTIES.mean_prop = "mean_property_value";


terms.model.MG_REV_PROPERTIES.predefined = {
    "Atchley" : {
             "Factor I bipolar":{
               "A":-0.591,
               "C":-1.343,
               "D":1.05,
               "E":1.357,
               "F":-1.006,
               "G":-0.384,
               "H":0.336,
               "I":-1.239,
               "K":1.831,
               "L":-1.019,
               "M":-0.663,
               "N":0.945,
               "P":0.189,
               "Q":0.931,
               "R":1.538,
               "S":-0.228,
               "T":-0.032,
               "V":-1.337,
               "W":-0.595,
               "Y":0.26
              },
             "Factor II secondary structure":{
               "A":-1.302,
               "C":0.465,
               "D":0.302,
               "E":-1.453,
               "F":-0.59,
               "G":1.652,
               "H":-0.417,
               "I":-0.547,
               "K":-0.5610000000000001,
               "L":-0.987,
               "M":-1.524,
               "N":0.828,
               "P":2.081,
               "Q":-0.179,
               "R":-0.055,
               "S":1.399,
               "T":0.326,
               "V":-0.279,
               "W":0.008999999999999999,
               "Y":0.83
              },
             "Factor III volume":{
               "A":-0.733,
               "C":-0.862,
               "D":-3.656,
               "E":1.477,
               "F":1.891,
               "G":1.33,
               "H":-1.673,
               "I":2.131,
               "K":0.533,
               "L":-1.505,
               "M":2.219,
               "N":1.299,
               "P":-1.628,
               "Q":-3.005,
               "R":1.502,
               "S":-4.76,
               "T":2.213,
               "V":-0.544,
               "W":0.672,
               "Y":3.097
              },
             "Factor IV composition":{
               "A":1.57,
               "C":-1.02,
               "D":-0.259,
               "E":0.113,
               "F":-0.397,
               "G":1.045,
               "H":-1.474,
               "I":0.393,
               "K":-0.277,
               "L":1.266,
               "M":-1.005,
               "N":-0.169,
               "P":0.421,
               "Q":-0.503,
               "R":0.44,
               "S":0.67,
               "T":0.908,
               "V":1.242,
               "W":-2.128,
               "Y":-0.838
              },
             "Factor V charge":{
               "A":-0.146,
               "C":-0.255,
               "D":-3.242,
               "E":-0.837,
               "F":0.412,
               "G":2.064,
               "H":-0.078,
               "I":0.8159999999999999,
               "K":1.648,
               "L":-0.912,
               "M":1.212,
               "N":0.9330000000000001,
               "P":-1.392,
               "Q":-1.853,
               "R":2.897,
               "S":-2.647,
               "T":1.313,
               "V":-1.262,
               "W":-0.184,
               "Y":1.512
              }
     },
     "3PROP": {
    "Hydrophobicity_KyteDoolittle": {
      "A": 1.8,
      "C": 2.5,
      "D": -3.5,
      "E": -3.5,
      "F": 2.8,
      "G": -0.4,
      "H": -3.2,
      "I": 4.5,
      "K": -3.9,
      "L": 3.8,
      "M": 1.9,
      "N": -3.5,
      "P": -1.6,
      "Q": -3.5,
      "R": -4.5,
      "S": -0.8,
      "T": -0.7,
      "V": 4.2,
      "W": -0.9,
      "Y": -1.3
    },
    "Volume_Angstrom3": {
      "A": 88.6,
      "C": 108.5,
      "D": 111.1,
      "E": 138.4,
      "F": 189.9,
      "G": 60.1,
      "H": 153.2,
      "I": 166.7,
      "K": 168.6,
      "L": 166.7,
      "M": 162.9,
      "N": 114.1,
      "P": 112.7,
      "Q": 143.8,
      "R": 173.4,
      "S": 89.0,
      "T": 116.1,
      "V": 140.0,
      "W": 227.8,
      "Y": 193.6
    },
    "Isoelectric_Point_pI": {
      "A": 6.0,
      "C": 5.07,
      "D": 2.77,
      "E": 3.22,
      "F": 5.48,
      "G": 5.97,
      "H": 7.59,
      "I": 6.02,
      "K": 9.74,
      "L": 5.98,
      "M": 5.74,
      "N": 5.41,
      "P": 6.3,
      "Q": 5.65,
      "R": 10.76,
      "S": 5.68,
      "T": 5.6,
      "V": 5.96,
      "W": 5.89,
      "Y": 5.66
    }
  },
    "4PROP": {
    "Hydrophobicity_KyteDoolittle": {
      "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
      "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
      "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
      "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3
    },
    "Volume_Angstrom3": {
      "A": 88.6, "C": 108.5, "D": 111.1, "E": 138.4, "F": 189.9,
      "G": 60.1, "H": 153.2, "I": 166.7, "K": 168.6, "L": 166.7,
      "M": 162.9, "N": 114.1, "P": 112.7, "Q": 143.8, "R": 173.4,
      "S": 89.0, "T": 116.1, "V": 140.0, "W": 227.8, "Y": 193.6
    },
    "Isoelectric_Point_pI": {
      "A": 6.00, "C": 5.07, "D": 2.77, "E": 3.22, "F": 5.48,
      "G": 5.97, "H": 7.59, "I": 6.02, "K": 9.74, "L": 5.98,
      "M": 5.74, "N": 5.41, "P": 6.30, "Q": 5.65, "R": 10.76,
      "S": 5.68, "T": 5.60, "V": 5.96, "W": 5.89, "Y": 5.66
    },
    "AlphaHelix_Propensity_ChouFasman": {
      "A": 1.42, "C": 0.70, "D": 1.01, "E": 1.51, "F": 1.13,
      "G": 0.57, "H": 1.00, "I": 1.08, "K": 1.16, "L": 1.21,
      "M": 1.45, "N": 0.67, "P": 0.57, "Q": 1.11, "R": 0.98,
      "S": 0.77, "T": 0.83, "V": 1.06, "W": 1.08, "Y": 0.69
    }
    },
    "5PROP": {
    "Hydrophobicity_KyteDoolittle": {
      "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
      "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
      "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
      "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3
    },
    "Volume_Angstrom3": {
      "A": 88.6, "C": 108.5, "D": 111.1, "E": 138.4, "F": 189.9,
      "G": 60.1, "H": 153.2, "I": 166.7, "K": 168.6, "L": 166.7,
      "M": 162.9, "N": 114.1, "P": 112.7, "Q": 143.8, "R": 173.4,
      "S": 89.0, "T": 116.1, "V": 140.0, "W": 227.8, "Y": 193.6
    },
    "Isoelectric_Point_pI": {
      "A": 6.00, "C": 5.07, "D": 2.77, "E": 3.22, "F": 5.48,
      "G": 5.97, "H": 7.59, "I": 6.02, "K": 9.74, "L": 5.98,
      "M": 5.74, "N": 5.41, "P": 6.30, "Q": 5.65, "R": 10.76,
      "S": 5.68, "T": 5.60, "V": 5.96, "W": 5.89, "Y": 5.66
    },
    "AlphaHelix_Propensity_ChouFasman": {
      "A": 1.42, "C": 0.70, "D": 1.01, "E": 1.51, "F": 1.13,
      "G": 0.57, "H": 1.00, "I": 1.08, "K": 1.16, "L": 1.21,
      "M": 1.45, "N": 0.67, "P": 0.57, "Q": 1.11, "R": 0.98,
      "S": 0.77, "T": 0.83, "V": 1.06, "W": 1.08, "Y": 0.69
    },
    "BetaSheet_Propensity_ChouFasman": {
      "A": 0.83, "C": 1.19, "D": 0.54, "E": 0.37, "F": 1.38,
      "G": 0.75, "H": 0.87, "I": 1.60, "K": 0.74, "L": 1.30,
      "M": 1.05, "N": 0.89, "P": 0.55, "Q": 1.10, "R": 0.93,
      "S": 0.75, "T": 1.19, "V": 1.70, "W": 1.37, "Y": 1.47
    }
  },
     "Random-2" : {
             "Random Factor 1":{
               "A":Random(-2,2),
               "C":Random(-2,2),
               "D":Random(-2,2),
               "E":Random(-2,2),
               "F":Random(-2,2),
               "G":Random(-2,2),
               "H":Random(-2,2),
               "I":Random(-2,2),
               "K":Random(-2,2),
               "L":Random(-2,2),
               "M":Random(-2,2),
               "N":Random(-2,2),
               "P":Random(-2,2),
               "Q":Random(-2,2),
               "R":Random(-2,2),
               "S":Random(-2,2),
               "T":Random(-2,2),
               "V":Random(-2,2),
               "W":Random(-2,2),
               "Y":Random(-2,2)
              },
             "Random Factor 2":{
               "A":Random(-2,2),
               "C":Random(-2,2),
               "D":Random(-2,2),
               "E":Random(-2,2),
               "F":Random(-2,2),
               "G":Random(-2,2),
               "H":Random(-2,2),
               "I":Random(-2,2),
               "K":Random(-2,2),
               "L":Random(-2,2),
               "M":Random(-2,2),
               "N":Random(-2,2),
               "P":Random(-2,2),
               "Q":Random(-2,2),
               "R":Random(-2,2),
               "S":Random(-2,2),
               "T":Random(-2,2),
               "V":Random(-2,2),
               "W":Random(-2,2),
               "Y":Random(-2,2)
              }
     },
     "Random-3" : {
             "Random Factor 1":{
               "A":Random(-2,2),
               "C":Random(-2,2),
               "D":Random(-2,2),
               "E":Random(-2,2),
               "F":Random(-2,2),
               "G":Random(-2,2),
               "H":Random(-2,2),
               "I":Random(-2,2),
               "K":Random(-2,2),
               "L":Random(-2,2),
               "M":Random(-2,2),
               "N":Random(-2,2),
               "P":Random(-2,2),
               "Q":Random(-2,2),
               "R":Random(-2,2),
               "S":Random(-2,2),
               "T":Random(-2,2),
               "V":Random(-2,2),
               "W":Random(-2,2),
               "Y":Random(-2,2)
              },
             "Random Factor 2":{
               "A":Random(-2,2),
               "C":Random(-2,2),
               "D":Random(-2,2),
               "E":Random(-2,2),
               "F":Random(-2,2),
               "G":Random(-2,2),
               "H":Random(-2,2),
               "I":Random(-2,2),
               "K":Random(-2,2),
               "L":Random(-2,2),
               "M":Random(-2,2),
               "N":Random(-2,2),
               "P":Random(-2,2),
               "Q":Random(-2,2),
               "R":Random(-2,2),
               "S":Random(-2,2),
               "T":Random(-2,2),
               "V":Random(-2,2),
               "W":Random(-2,2),
               "Y":Random(-2,2)
              },
             "Random Factor 3":{
               "A":Random(-2,2),
               "C":Random(-2,2),
               "D":Random(-2,2),
               "E":Random(-2,2),
               "F":Random(-2,2),
               "G":Random(-2,2),
               "H":Random(-2,2),
               "I":Random(-2,2),
               "K":Random(-2,2),
               "L":Random(-2,2),
               "M":Random(-2,2),
               "N":Random(-2,2),
               "P":Random(-2,2),
               "Q":Random(-2,2),
               "R":Random(-2,2),
               "S":Random(-2,2),
               "T":Random(-2,2),
               "V":Random(-2,2),
               "W":Random(-2,2),
               "Y":Random(-2,2)
              }
     },
     "Random-4" : {
             "Random Factor 1":{
               "A":Random(-2,2),
               "C":Random(-2,2),
               "D":Random(-2,2),
               "E":Random(-2,2),
               "F":Random(-2,2),
               "G":Random(-2,2),
               "H":Random(-2,2),
               "I":Random(-2,2),
               "K":Random(-2,2),
               "L":Random(-2,2),
               "M":Random(-2,2),
               "N":Random(-2,2),
               "P":Random(-2,2),
               "Q":Random(-2,2),
               "R":Random(-2,2),
               "S":Random(-2,2),
               "T":Random(-2,2),
               "V":Random(-2,2),
               "W":Random(-2,2),
               "Y":Random(-2,2)
              },
             "Random Factor 2":{
               "A":Random(-2,2),
               "C":Random(-2,2),
               "D":Random(-2,2),
               "E":Random(-2,2),
               "F":Random(-2,2),
               "G":Random(-2,2),
               "H":Random(-2,2),
               "I":Random(-2,2),
               "K":Random(-2,2),
               "L":Random(-2,2),
               "M":Random(-2,2),
               "N":Random(-2,2),
               "P":Random(-2,2),
               "Q":Random(-2,2),
               "R":Random(-2,2),
               "S":Random(-2,2),
               "T":Random(-2,2),
               "V":Random(-2,2),
               "W":Random(-2,2),
               "Y":Random(-2,2)
              },
             "Random Factor 3":{
               "A":Random(-2,2),
               "C":Random(-2,2),
               "D":Random(-2,2),
               "E":Random(-2,2),
               "F":Random(-2,2),
               "G":Random(-2,2),
               "H":Random(-2,2),
               "I":Random(-2,2),
               "K":Random(-2,2),
               "L":Random(-2,2),
               "M":Random(-2,2),
               "N":Random(-2,2),
               "P":Random(-2,2),
               "Q":Random(-2,2),
               "R":Random(-2,2),
               "S":Random(-2,2),
               "T":Random(-2,2),
               "V":Random(-2,2),
               "W":Random(-2,2),
               "Y":Random(-2,2)
              },
             "Random Factor 4":{
               "A":Random(-2,2),
               "C":Random(-2,2),
               "D":Random(-2,2),
               "E":Random(-2,2),
               "F":Random(-2,2),
               "G":Random(-2,2),
               "H":Random(-2,2),
               "I":Random(-2,2),
               "K":Random(-2,2),
               "L":Random(-2,2),
               "M":Random(-2,2),
               "N":Random(-2,2),
               "P":Random(-2,2),
               "Q":Random(-2,2),
               "R":Random(-2,2),
               "S":Random(-2,2),
               "T":Random(-2,2),
               "V":Random(-2,2),
               "W":Random(-2,2),
               "Y":Random(-2,2)
              }
     },
     "Random-5" : {
             "Random Factor 1":{
               "A":Random(-2,2),
               "C":Random(-2,2),
               "D":Random(-2,2),
               "E":Random(-2,2),
               "F":Random(-2,2),
               "G":Random(-2,2),
               "H":Random(-2,2),
               "I":Random(-2,2),
               "K":Random(-2,2),
               "L":Random(-2,2),
               "M":Random(-2,2),
               "N":Random(-2,2),
               "P":Random(-2,2),
               "Q":Random(-2,2),
               "R":Random(-2,2),
               "S":Random(-2,2),
               "T":Random(-2,2),
               "V":Random(-2,2),
               "W":Random(-2,2),
               "Y":Random(-2,2)
              },
             "Random Factor 2":{
               "A":Random(-2,2),
               "C":Random(-2,2),
               "D":Random(-2,2),
               "E":Random(-2,2),
               "F":Random(-2,2),
               "G":Random(-2,2),
               "H":Random(-2,2),
               "I":Random(-2,2),
               "K":Random(-2,2),
               "L":Random(-2,2),
               "M":Random(-2,2),
               "N":Random(-2,2),
               "P":Random(-2,2),
               "Q":Random(-2,2),
               "R":Random(-2,2),
               "S":Random(-2,2),
               "T":Random(-2,2),
               "V":Random(-2,2),
               "W":Random(-2,2),
               "Y":Random(-2,2)
              },
             "Random Factor 3":{
               "A":Random(-2,2),
               "C":Random(-2,2),
               "D":Random(-2,2),
               "E":Random(-2,2),
               "F":Random(-2,2),
               "G":Random(-2,2),
               "H":Random(-2,2),
               "I":Random(-2,2),
               "K":Random(-2,2),
               "L":Random(-2,2),
               "M":Random(-2,2),
               "N":Random(-2,2),
               "P":Random(-2,2),
               "Q":Random(-2,2),
               "R":Random(-2,2),
               "S":Random(-2,2),
               "T":Random(-2,2),
               "V":Random(-2,2),
               "W":Random(-2,2),
               "Y":Random(-2,2)
              },
             "Random Factor 4":{
               "A":Random(-2,2),
               "C":Random(-2,2),
               "D":Random(-2,2),
               "E":Random(-2,2),
               "F":Random(-2,2),
               "G":Random(-2,2),
               "H":Random(-2,2),
               "I":Random(-2,2),
               "K":Random(-2,2),
               "L":Random(-2,2),
               "M":Random(-2,2),
               "N":Random(-2,2),
               "P":Random(-2,2),
               "Q":Random(-2,2),
               "R":Random(-2,2),
               "S":Random(-2,2),
               "T":Random(-2,2),
               "V":Random(-2,2),
               "W":Random(-2,2),
               "Y":Random(-2,2)
              },
             "Random Factor 5":{
               "A":Random(-2,2),
               "C":Random(-2,2),
               "D":Random(-2,2),
               "E":Random(-2,2),
               "F":Random(-2,2),
               "G":Random(-2,2),
               "H":Random(-2,2),
               "I":Random(-2,2),
               "K":Random(-2,2),
               "L":Random(-2,2),
               "M":Random(-2,2),
               "N":Random(-2,2),
               "P":Random(-2,2),
               "Q":Random(-2,2),
               "R":Random(-2,2),
               "S":Random(-2,2),
               "T":Random(-2,2),
               "V":Random(-2,2),
               "W":Random(-2,2),
               "Y":Random(-2,2)
              }
     },
     "2PROP": {
            "Hydrophobicity_KyteDoolittle": {
              "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
              "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
              "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
              "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3
            },
            "Volume_Angstrom3": {
              "A": 88.6, "C": 108.5, "D": 111.1, "E": 138.4, "F": 189.9,
              "G": 60.1, "H": 153.2, "I": 166.7, "K": 168.6, "L": 166.7,
              "M": 162.9, "N": 114.1, "P": 112.7, "Q": 143.8, "R": 173.4,
              "S": 89.0, "T": 116.1, "V": 140.0, "W": 227.8, "Y": 193.6
            }
    }
};

//----------------------------------------------------------------------------------------------------------------


lfunction model.codon.MG_REV_PROPERTIES.prompt_and_define (type, code) {
    KeywordArgument ("property-set", "How to partition synonymous codons into classes", "3PROP");
    
    property_set = io.SelectAnOption (
            {
                "Atchley":"Use the five properties derived from a factor analysis of 500 amino-acid properties [Table 2 in PNAS (2005) 102(18) 6395-6400 doi: 10.1073/pnas.0408677102]",
                "LCAP":"Use the five properties defined in the Conant and Stadler LCAP model [Mol Biol Evol (2009) 26 (5): 1155-1161. doi: 10.1093/molbev/msp031]",
                 "2PROP" : "Use two primary properties: Hydrophobicity (Kyle-Doolitle), Volume (cubic angstroms)",
                 "3PROP" : "Use three primary properties: Hydrophobicity (Kyle-Doolitle), Volume (cubic angstroms), Isoelectric Point (pI)",
                "4PROP" : "Use four primary properties: Hydrophobicity (Kyle-Doolitle), Volume (cubic angstroms), Isoelectric Point (pI), Alpha-Helix Propensity (Chou-Fasman)",
                "5PROP" : "Use five primary properties: Hydrophobicity (Kyle-Doolitle), Volume (cubic angstroms), Isoelectric Point (pI), Alpha-Helix Propensity (Chou-Fasman), Beta-Sheet Propensity (Chou-Fasman)",
                "Random-2" : "Two random properties (for null hypothesis testing)",
                "Random-3" : "Three random properties (for null hypothesis testing)",
                "Random-4" : "Four random properties (for null hypothesis testing)",
                "Random-5" : "Five random properties (for null hypothesis testing)",
                "Custom":    "Load the set of properties from a file"
            }, 
            "The set of properties to use in the model");

    
    
    if (property_set == "Custom") {
        KeywordArgument ("property-file", "JSON file which defines amino-acid properties");
        property_set = io.PromptUserForFilePathRead ("JSON file which defines amino-acid properties");
        property_set = io.ParseJSON(property_set);
        console.log (">Loaded a set of `Abs(property_set)` properties");
     }
    
     
    
     return models.codon.MG_REV_PROPERTIES.ModelDescription(type, code, property_set);
}


lfunction models.codon.MG_REV_PROPERTIES.ModelDescription(type, code, properties) {


    // piggyback on the standard MG_REV model for most of the code

    mg_base = models.codon.MG_REV.ModelDescription (type, code);
    mg_base[utility.getGlobalValue("terms.description")] = "The Muse-Gaut 94 codon-substitution model coupled with the general time reversible (GTR) model of nucleotide substitution, which incorporates amino-acid residue properties into the non-synonymous rates";
    mg_base[utility.getGlobalValue("terms.model.q_ij")] = "models.codon.MG_REV_PROPERTIES._GenerateRate";
    mg_base[utility.getGlobalValue("terms.model.residue_properties")] = models.codon.MG_REV_PROPERTIES._munge_properties(properties);
    mg_base[utility.getGlobalValue("terms.model.residue_name_map")] = parameters.ValidateIDs (utility.Keys (mg_base[utility.getGlobalValue("terms.model.residue_properties")]));
    mg_base[utility.getGlobalValue("terms.model.post_definition")] = "models.codon.MG_REV_PROPERTIES.post_definition";
    mg_base[utility.getGlobalValue("terms.model.set_branch_length")] = "models.codon.MG_REV_PROPERTIES.set_branch_length";
        
    return mg_base;
}

/**
 * @name models.codon.MG_REV_PROPERTIES._munge_properties
 * @param {String or AssociativeList} properties
 * @return {AssociativeList}
 * Convert a set of amino-acid properties into the expected format
 * This will also recenter the properties (mean 0, variance of 1)
 */


lfunction models.codon.MG_REV_PROPERTIES._munge_properties (properties) {

    if (Type (properties) == "String") {
        assert (^"terms.model.MG_REV_PROPERTIES.predefined" / properties, properties + " is not a valid predefined property set");
        return models.codon.MG_REV_PROPERTIES._munge_properties ((^"terms.model.MG_REV_PROPERTIES.predefined")[properties]);
    }

    assert (Type (properties) == "AssociativeList", "Properties definition must be an AssociativeArray or a String for predefined properties");

        /*
             each property must map to either a dictionary with '1 letter AA' -> value
             or to an array with 20 values
        */


    valid_aa = utility.MatrixToDict (^"models.protein.alphabet");

    for (key, value; in; properties) {
        if (Type (value) == "AssociativeList") {
            io.CheckAssertion ("utility.Array1D(`&value`)==20", "A dictionary of amino-acid properties must have 20 entries for " + key);
            for (_i = 0; _i < 20; _i += 1) {
                assert (value / (^"models.protein.alphabet")[_i],
                    "Property " + key + " is not defined for amino-acid " + (^"models.protein.alphabet")[_i]);
            }
        } else {
            if (Type (value) == "Matrix") {
                io.CheckAssertion ("utility.Array1D(`&value`)==20", "A matrix of amino-acid properties must have 20 entries for " + key);
                value.dict = {};
                for (_i = 0; _i < 20; _i += 1) {
                    value.dict [(^"models.protein.alphabet")[_i]] = value[_i];
                }
                properties[key] = value.dict;
            } else {
                assert (0, "Invalid entry for properties in class " + key);
            }
        }
    }
        
        
    for (key, value; in; properties) {
        sum  = 0;
        sum2 = 0;
        
        for (aa, prop; in; value) {
            sum += prop;
            sum2 += prop*prop;
        }
        
        sum = sum / 20;
        norm = Sqrt (sum2 / 20 - sum*sum);
        
        for (aa, prop; in; value) {
            value [aa] = (prop-sum) / norm;
        }       
    }
        
    assert (utility.Array1D(properties) > 0, "At least one valid property set must be defined");
    return properties;
}

lfunction models.codon.MG_REV_PROPERTIES._GenerateRate(fromChar, toChar, namespace, model_type, model) {
    return models.codon.MG_REV_PROPERTIES._GenerateRate_generic (fromChar, toChar, namespace, model_type,
        model[utility.getGlobalValue("terms.translation_table")],
        "alpha", utility.getGlobalValue("terms.parameters.synonymous_rate"),
        "beta",  utility.getGlobalValue("terms.parameters.nonsynonymous_rate"),
        "lambda", "",
        "omega", utility.getGlobalValue("terms.parameters.log_omega_ratio"),
        model[utility.getGlobalValue("terms.model.residue_properties")],
        model[utility.getGlobalValue("terms.model.residue_name_map")]
        );
}

/**
 * @name models.codon.MG_REV_PROPERTIES._GenerateRate
 * @param {Number} fromChar
 * @param {Number} toChar
 * @param {String} namespace
 * @param {String} model_type
 * @param {Matrix} _tt - translation table
 */


lfunction models.codon.MG_REV_PROPERTIES._GenerateRate_generic (fromChar, toChar, namespace, model_type, _tt, alpha, alpha_term, beta, beta_term, lambda, lambda_term, omega, omega_term, properties, property_id_map) {

    _GenerateRate.p = {};
    _GenerateRate.diff = models.codon.diff.complete(fromChar, toChar);
    diff_count = utility.Array1D (_GenerateRate.diff);

    if (diff_count == 1) {

        _GenerateRate.p[model_type] = {};
        _GenerateRate.p[utility.getGlobalValue("terms.global")] = {};

        nuc_rate = "";

        for (i = 0; i < diff_count; i += 1) {
            if ((_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.from")] > (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.to")]) {
                nuc_p = "theta_" + (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.to")] + (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.from")];
            } else {
                nuc_p = "theta_" + (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.from")] +(_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.to")];
            }
            nuc_p = parameters.ApplyNameSpace(nuc_p, namespace);
            (_GenerateRate.p[utility.getGlobalValue("terms.global")])[terms.nucleotideRateReversible((_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.from")], (_GenerateRate.diff[i])[utility.getGlobalValue("terms.diff.to")])] = nuc_p;

            nuc_rate = parameters.AppendMultiplicativeTerm (nuc_rate, nuc_p);
       }

        rate_entry = nuc_rate;

        if (_tt[fromChar] != _tt[toChar]) {
            if (model_type == utility.getGlobalValue("terms.global")) {
                aa_rate = {};
                prop_count = Abs (properties);
                prop_names = utility.Keys (properties);
                (_GenerateRate.p[model_type]) [omega_term] = parameters.ApplyNameSpace (omega, namespace);

                for (prop_name, prop_values; in; properties) {
                    prop_diff = Abs (prop_values[_tt[fromChar]] - prop_values[_tt[toChar]]);
                    term_rate = parameters.ApplyNameSpace(lambda + "_" + property_id_map[prop_name], namespace);
                    (_GenerateRate.p[model_type])[terms.propertyImportance (prop_name, lambda_term)] = term_rate;
                    aa_rate + ("`term_rate`*" + prop_diff);
                }

                //aa_rate = parameters.ApplyNameSpace(lambda, namespace);
                rate_entry += "*Min(10000,Exp(" + (_GenerateRate.p[model_type]) [omega_term] + "-(" + Join("+",aa_rate) + ")))";
             } else {
                aa_rate = {};
                prop_count = Abs (properties);
                prop_names = utility.Keys (properties);
                (_GenerateRate.p[model_type]) [beta_term] = beta;

                for (prop_name, prop_values; in; properties) {
                    prop_diff = Abs (prop_values[_tt[fromChar]] - prop_values[_tt[toChar]]);
                    term_rate = lambda + "_" + property_id_map[prop_name];
                    (_GenerateRate.p[model_type])[terms.propertyImportance (prop_name, lambda_term)] = term_rate;
                    aa_rate + ("`term_rate`*" + prop_diff);
                }
                rate_entry += "*Exp(-(" + Join("+",aa_rate) + "))";
                (_GenerateRate.p[model_type])[alpha_term] = alpha;
                rate_entry = "Min(10000," + rate_entry  + "*" + beta + ")";
 
             }
        } else {
            if (model_type == utility.getGlobalValue("terms.local")) {
                 (_GenerateRate.p[model_type])[alpha_term] = alpha;
                rate_entry += "*" + alpha;
           } else {
                _GenerateRate.p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate;
            }
        }

        _GenerateRate.p[utility.getGlobalValue("terms.model.rate_entry")] = rate_entry;
    }
    return _GenerateRate.p;
}

lfunction  models.codon.MG_REV_PROPERTIES.post_definition (model) {
    prop_range = {
        ^"terms.lower_bound": "-10",
        ^"terms.upper_bound": "10"
    };


    for (id; in ; model.GetParameters_RegExp (model, terms.propertyImportance ('', '') + "|" + ^"terms.parameters.log_omega_ratio")) {
        parameters.SetRange(id, prop_range);
        parameters.SetValue(id, 0.1);
    }

    models.generic.post.definition (model);
}

lfunction models.codon.MG_REV_PROPERTIES.set_branch_length(model, value, parameter) {
    
    if (utility.Has (model, ^"terms.model.MG_REV_PROPERTIES.mean_prop", "Number")) {
      properties = model.GetLocalParameters_RegExp(model, terms.propertyImportance ("",""));
      for (tag, id; in; properties) {
        parameters.SetValue (id, model[^"terms.model.MG_REV_PROPERTIES.mean_prop"]);
      }
      if (utility.Has (model, "fraction_same", "Number")) {
        fs = model["fraction_same"];
        models.codon.MG_REV.set_branch_length(model,value,parameter);
        weighted = model ["models.codon.MG_REV.set_branch_length"];
        if (Type (weighted) == "AssociativeList") {
          for (tag, id; in; properties) {
            parameters.SetValue (id, 0.);
          }
          models.codon.MG_REV.set_branch_length(model,value,parameter);
          unweighted = model ["models.codon.MG_REV.set_branch_length"];
          for (tag, value; in; weighted) {
            //console.log (fs * unweighted[tag] + (1-fs) * value);
            parameters.SetValue (tag, fs * unweighted[tag] + (1-fs) * value);
          }
          return 0;
        }
      }
    }
    return models.codon.MG_REV.set_branch_length(model,value,parameter);

}
