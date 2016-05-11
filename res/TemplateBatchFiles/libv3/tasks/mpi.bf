LoadFunctionLibrary ("../terms-json.bf");

namespace mpi {
    job_id = 0;

    function get_job_id () {
        job_id += 1;
        return job_id;
    }


    lfunction create_queue (nodesetup) {
        /** create and return an empty FIFO queue for MPI jobs */
        mpi_node_count = utility.GetEnvVariable ("MPI_NODE_COUNT");
    
        queue = {};
        send_to_nodes = ""; 
        if (mpi_node_count > 1) {
    
    
            if (None != nodesetup) {
                if (Abs (nodesetup)) {
                
                    utility.SetEnvVariable ("LF_NEXUS_EXPORT_EXTRA", 
                                            'PRESERVE_SLAVE_NODE_STATE = TRUE; MPI_NEXUS_FILE_RETURN = "None";');
                
                    send_to_nodes * 128;
                
                    utility.ForEach (nodesetup["LikelihoodFunctions"], "_value_", 
                                     '
                                        ExecuteCommands ("Export (create_queue.temp, " + _value_ + ")");
                                        for (`&k` = 1; `&k` < `mpi_node_count`; `&k` += 1) {
                                    
                                            MPISend (`&k`, create_queue.temp);
                                            MPIReceive (-1, ignore, ignore);
                                        }
                                    
                                 
                                     ');
                
                    send_to_nodes * 0;
                }
            }
    
            for (k = 1; k < mpi_node_count; k += 1) {
                queue [k] = {"job_id" : None, "callback" : None, "arguments": None};
            }
        }
        return queue;
    }


    lfunction queue_job (queue, job, arguments, result_callback) {
        /** 
            send the job function with provided arguments to 
            the first available node.
        
            When the job is finished; call the "result_callback" function 
        
        */
    
        mpi_node_count = utility.GetEnvVariable ("MPI_NODE_COUNT");
    
        if (mpi_node_count > 1) {    
            for (node = 1; node < mpi_node_count; node += 1) {
                if (None == (queue [node])["job_id"]) {
                    break;
                }
            }
    
            if (node == mpi_node_count) {
                node = aux._handle_receieve (queue);        
            }

        
            complete_function_dump = aux.queue_export_function (job);
            job_id = get_job_id();   
            queue [node] = {"job_id" : job_id, "callback" : result_callback, "arguments" : arguments};
            MPISend (node, complete_function_dump + "; return " + job + '(' + Join (",",utility.Map (arguments,"_value_", "utility.convertToArgumentString (_value_)")) + ')');    
        } else {
            Call (result_callback, 0, Eval (job + '(' + Join (",",utility.Map (arguments,"_value_", "utility.convertToArgumentString (_value_)")) + ')'), arguments);
        }
    }

    lfunction queue_complete (queue) {
       
        mpi_node_count = utility.GetEnvVariable ("MPI_NODE_COUNT");
    
        if (mpi_node_count > 1) {
            do {
    
                for (node = 1; node < mpi_node_count; node += 1) {
                    if (None != (queue [node])["job_id"]) {
                        break;
                    }
                }
    
                if (node < mpi_node_count) {
                    node = aux._handle_receieve (queue);        
                }
            } while (node < mpi_node_count);
        }
    }

    namespace aux {
        function queue_export_function (func_id) {
            Export (complete_function_dump, ^func_id);
            return complete_function_dump;
        }

        lfunction _handle_receieve (queue) {
            MPIReceive (-1,from,result);
            Call ((queue [from])["callback"], from, Eval(result), (queue [from])["arguments"]);
            queue [from] = {"job_id" : None, "callback" : None};
            return from; 
        }
    }
}
