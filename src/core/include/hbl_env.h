/*
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon42@uwo.ca)
 Steven Weaver (sweaver@temple.edu)
 
 Module Developers:
 Lance Hepler (nlhepler@gmail.com)
 Martin Smith (martin.audacis@gmail.com)
 
 Significant contributions from:
 Spencer V Muse (muse@stat.ncsu.edu)
 Simon DW Frost (sdf22@cam.ac.uk)
 
 Permission is hereby granted, free of charge, to any person obtaining a
 copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject to
 the following conditions:
 
 The above copyright notice and this permission notice shall be included
 in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

#include "hy_strings.h"
#include "avllistxl.h"
#include "parser.h"
#include "defines.h"



namespace hy_env {
 
 
    bool       EnvVariableTrue (_String const& name);
    /**
         Check if the value of the environment variable is "true"

         @param name the name of the env variable (in HBL)
         
         @return TRUE if the variable is defined and true-ike (not zero)
     
     */

    
    hyFloat       EnvVariableGetDefaultNumber (_String const& name);
    /**
     Look up the default value _numeric_ for a given env variable,
     
     @param name the name of the env variable (in HBL)
     
     @return default value if available and of the correct type, otherwise HY_INVALID_RETURN_VALUE
     
     */

    
    _PMathObj       EnvVariableGetDefault (_String const& name, unsigned long type);
    /**
         Look up the default value for a given env variable,
              checking that it is of a paricular type (could be any)
         
         @param name the name of the env variable (in HBL)
         @param type check for return type
         
         @return default value if available and of the correct type, otherwise nil
     
     */
    
    _PMathObj       EnvVariableGet (_String const& name, unsigned long type);
    /**
     Look up the value for a given env variable,
     checking that it is of a paricular type (could be any)
     
     @param name the name of the env variable (in HBL)
     @param type check for return type
     
     @return current value if set and of the correct type, otherwise default value if available and of the correct type, otherwise nil
     
     */

    void            EnvVariableSet (_String const& name, _PMathObj value, bool copy );
    /**
     Set the value for a env variable (global scope)
     copying the value if requested
     
     @param name the name of the env variable (in HBL)
     @param value the value to set the variable to
     @param copy is set to true, make a deep copy of 'value' prior to copying
     
     
     */

    void            EnvVariableSetNamespace (_String const& name, _PMathObj value, _String* nmspace, bool copy );
    /**
     Set the value for a env variable (in a given namespace if provided)
     copying the value if requested
     
     @param name the name of the env variable (in HBL)
     @param value the value to set the variable to
     @param copy is set to true, make a deep copy of 'value' prior to copying
     
     
     */


    void       SetupEnvDefaults (void);
    /**
        Populate default value arrays
     */

    extern const _String
    
          data_file_tree,
          data_file_tree_string,
          nexus_file_tree_matrix      ,
          data_file_partition_matrix  ,
          use_traversal_heuristic    ,
          normalize_sequence_names   ,
          dataset_save_memory_size,
          sitewise_matrix,
          blockwise_matrix,
          covariance_parameter,
          selection_strings,
          random_seed,
          assigned_seed,
          base_directory,
          lib_directory,
          directory_separator_char,
          path_to_current_bf,
          true_const,
          false_const,
          harvest_frequencies_gap_options,
          last_file_path,
          matrix_element_row,
          matrix_element_column,
          matrix_element_value,
          message_logging,
          mpi_node_id,
          mpi_node_count,
          mpi_last_sent_message,
          error_report_format_expression,
          error_report_format_expression_string,
          error_report_format_expression_stack,
          error_report_format_expression_stdin,
          status_bar_update_string,
          use_last_model,
          last_model_parameter_list,
          kGetStringFromUser,
          kPairwiseCountAmbiguitiesResolve,
          kPairwiseCountAmbiguitiesAverage,
          kPairwiseCountAmbiguitiesSkip,
          get_data_info_returns_only_the_index;

    ;
  
    
};
