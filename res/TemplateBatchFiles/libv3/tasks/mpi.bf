LoadFunctionLibrary ("../terms-json.bf");

/** @module mpi
        Functions for creating, populating, and manipulating
        MPI job queues. In the absence of an MPI environment
        the jobs are executed serially
 */



namespace mpi {
    job_id = 0;

    function get_next_job_id () {
        job_id += 1;
        return job_id;
    }


    lfunction CreateQueue (nodesetup) {
        /** create and return an empty FIFO queue for MPI jobs
         * @name mpi.CreateQueue
         * @param  {Dict} nodesetup
         *      controls what gets passed to slave nodes
         *      "Headers" -> iterable (matrix/dict) of string paths of header files to load
         *      "Models" ->  matrix of model names to make available to slave nodes
         *      "LikelihoodFunctions" -> iterable (matrix/dict) of LikelihoodFunction IDs to export to slave nodes
         * @return {Dict} an "opaque" queue structure
         */

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
                                        }
                                        for (`&k` = 1; `&k` < `mpi_node_count`; `&k` += 1) {
                                            MPIReceive (-1, ignore, ignore);
                                        }

                                     ');


                    if (utility.Has (nodesetup, "Headers", None)) {
                        send_to_nodes * "PRESERVE_SLAVE_NODE_STATE = TRUE;\n";
                        send_to_nodes * (Join (";\n",utility.Map (nodesetup["Headers"], "_value_", "'LoadFunctionLibrary(\"' + _value_ +'\")'")) + ";");
                    }

                    if (utility.Has (nodesetup, "Functions", None)) {
                        send_to_nodes * "PRESERVE_SLAVE_NODE_STATE = TRUE;\n";
                        utility.ForEach (nodesetup["Functions"], "_value_",
                            '
                                ExecuteCommands ("Export (_test_id_," + _value_ + ")");
                                `&send_to_nodes` * _test_id_;
                            '
                        );
                    }

                    if (utility.Has (nodesetup, "Variables", None)) {
                        utility.ForEach (nodesetup["Variables"], "_value_",
                            '
                                `&send_to_nodes` * ("\n" + _value_ + " = " +  (parameters.Quote(^_value_)) + ";\n") ;
                            '
                        );
                    }

                    model_count = utility.Array1D (nodesetup["Models"]);

                    if (model_count) {
                        send_to_nodes * "PRESERVE_SLAVE_NODE_STATE = TRUE;\n";

                        globals_to_export = {};
                        functions_to_export = {};

                        for (m = 0; m < model_count; m+=1) {
                            model_name = (nodesetup["Models"])[m];
                            model_globals = utility.Values(((^model_name)["parameters"])[^"terms.global"]);
                            model_global_count = utility.Array1D (model_globals);
                            for (v = 0; v < model_global_count; v+=1) {
                                globals_to_export [model_globals[v]] = 1;
                            }

                            utility.ForEach ({{"get-branch-length","set-branch-length"}}, "_value_",
                            '
                                _test_id_ = (^(`&model_name`))[_value_];
                                if (Type (_test_id_) == "String" && Abs (_test_id_) > 0) {
                                    `&functions_to_export` [_test_id_] = 1;
                                }
                            ');
                        }

                        utility.ForEach (utility.Keys(globals_to_export), "_value_",
                            '
                                `&send_to_nodes` * parameters.ExportParameterDefinition (_value_);
                            '
                        );

                        utility.ForEach (utility.Keys(functions_to_export), "_value_",
                            '
                                ExecuteCommands ("Export (_test_id_," + _value_ + ")");
                                `&send_to_nodes` * _test_id_;
                            '
                        );
                    }


                    send_to_nodes * 0;
                }
            }

            if (Abs (send_to_nodes)) {
            	
            	
                for (k = 1; k < mpi_node_count; k += 1) {
                   MPISend (k, send_to_nodes);
                }
                for (k = 1; k < mpi_node_count; k += 1) {
                   MPIReceive (-1, ignore, ignore);
                }
            }

            //assert (0);
            for (k = 1; k < mpi_node_count; k += 1) {
                queue [k] = {"job_id" : None, "callback" : None, "arguments": None};
            }
        }
        return queue;
    }


    lfunction QueueJob (queue, job, arguments, result_callback) {
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
            //console.log (complete_function_dump);
            job_id = get_next_job_id();
            //fprintf (stdout, "Sending to node ", node, "\n");
            queue [node] = {"job_id" : job_id, "callback" : result_callback, "arguments" : arguments};
             MPISend (node, complete_function_dump + "; return " + job + '(' + Join (",",utility.Map (arguments,"_value_", "utility.convertToArgumentString (_value_)")) + ')');

        } else {
            Call (result_callback, 0, Eval (job + '(' + Join (",",utility.Map (arguments,"_value_", "utility.convertToArgumentString (_value_)")) + ')'), arguments);
        }
    }

    lfunction QueueComplete (queue) {

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
