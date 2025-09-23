# HyPhy Batch Language (HBL) Commands

This document provides a comprehensive reference for the commands available in the HyPhy Batch Language (HBL). These commands are used to control various aspects of phylogenetic analysis, from data manipulation to model execution.

## Command Reference

### `function`

* **Syntax:** `function <name>(<args>) { ... }`
* **Description:** Defines a user function. This allows you to encapsulate a block of code that can be called by name.
* **Arguments:**
    * `name`: The name of the function.
    * `args`: A comma-separated list of arguments for the function.

### `ffunction`

* **Syntax:** `ffunction <name>(<args>) { ... }`
* **Description:** Defines a user function that returns a value. This is similar to `function`, but the function is expected to return a value using the `return` command.
* **Arguments:**
    * `name`: The name of the function.
    * `args`: A comma-separated list of arguments for the function.

### `return`

* **Syntax:** `return <value>;` or `return(<value>);`
* **Description:** Returns a value from a function defined with `ffunction`.
* **Arguments:**
    * `value`: The value to be returned.

### `if`

* **Syntax:** `if (<condition>) { ... }`
* **Description:** Executes a block of code if a condition is true.
* **Arguments:**
    * `condition`: An expression that evaluates to a boolean value.

### `else`

* **Syntax:** `else { ... }`
* **Description:** Executes a block of code if the preceding `if` condition is false.

### `do`

* **Syntax:** `do { ... } while (<condition>);`
* **Description:** Executes a block of code at least once, and then repeatedly as long as a condition is true.
* **Arguments:**
    * `condition`: An expression that evaluates to a boolean value.

### `break`

* **Syntax:** `break;`
* **Description:** Exits a loop (e.g., `for`, `while`, `do-while`).

### `continue`

* **Syntax:** `continue;`
* **Description:** Skips the current iteration of a loop and proceeds to the next one.

### `#include`

* **Syntax:** `#include "<file_path>"`
* **Description:** Includes the contents of another HBL file. This is useful for organizing code into multiple files.
* **Arguments:**
    * `file_path`: The path to the HBL file to include.

### `DataSet`

* **Syntax:** `DataSet <name> = <data>;`
* **Description:** Defines a new dataset. The data can be provided in various formats.
* **Arguments:**
    * `name`: The name of the dataset.
    * `data`: The dataset content.

### `DataSetFilter`

* **Syntax:** `DataSetFilter <name> = <filter_definition>;`
* **Description:** Defines a new dataset filter, which can be used to select a subset of a dataset.
* **Arguments:**
    * `name`: The name of the dataset filter.
    * `filter_definition`: The definition of the filter.

### `ConstructCategoryMatrix`

* **Syntax:** `ConstructCategoryMatrix(<receptacle>, <Likelihood Function|Tree>, [optional <COMPLETE|SHORT|WEIGHTS|CLASSES (default = COMPLETE)> , matrix argument with partitions to include (default = all)>])`
* **Description:** Constructs a category matrix, which is used to define different rate categories for sites in a sequence alignment.
* **Arguments:**
    * `receptacle`: The variable where the resulting matrix will be stored.
    * `Likelihood Function|Tree`: The likelihood function or tree to use for constructing the matrix.
    * `COMPLETE|SHORT|WEIGHTS|CLASSES`: The type of matrix to construct.
    * `matrix argument`: An optional argument specifying which partitions to include.

### `Tree`

* **Syntax:** `Tree <name> = <newick_string>;`
* **Description:** Defines a new tree from a Newick-formatted string.
* **Arguments:**
    * `name`: The name of the tree.
    * `newick_string`: The tree in Newick format.

### `LikelihoodFunction`

* **Syntax:** `LikelihoodFunction <name> = (<DataSet>, <Tree>);`
* **Description:** Defines a new likelihood function, which is the core of many phylogenetic analyses.
* **Arguments:**
    * `name`: The name of the likelihood function.
    * `DataSet`: The dataset to be used.
    * `Tree`: The tree to be used.

### `LikelihoodFunction3`

* **Syntax:** `LikelihoodFunction3 <name> = (<DataSet>, <Tree>, <Model>);`
* **Description:** Defines a new likelihood function with three components: a dataset, a tree, and a model.
* **Arguments:**
    * `name`: The name of the likelihood function.
    * `DataSet`: The dataset to be used.
    * `Tree`: The tree to be used.
    * `Model`: The substitution model to be used.

### `MolecularClock`

* **Syntax:** `MolecularClock(tree or tree node, local variable 1 [optional ,<local variable 2>, ..., <local variable N>])`
* **Description:** Enforces a molecular clock on a tree, which assumes that the rate of evolution is constant across all branches.
* **Arguments:**
    * `tree or tree node`: The tree or a specific node in the tree to apply the clock to.
    * `local variable`: One or more local variables to be used in the clock model.

### `fscanf`

* **Syntax:** `fscanf(file path (string),<optional 'REWIND'>,'type 1 (s|S|n|N|m|M|h|H|f|F)')[optional , <type 2>, ... <type N>], var1 [optional , <var 2>, ... <var N>])`
* **Description:** Reads formatted data from a file.
* **Arguments:**
    * `file path`: The path to the file to read from.
    * `REWIND`: An optional keyword to rewind the file pointer to the beginning of the file before reading.
    * `type`: The type of data to read (e.g., `s` for string, `n` for number).
    * `var`: The variable to store the read data in.

### `sscanf`

* **Syntax:** `sscanf(string,<optional 'REWIND'>,'type 1 (s|S|n|N|m|M|h|H|f|F)')[optional , <type 2>, ... <type N>], var1 [optional , <var 2>, ... <var N>])`
* **Description:** Reads formatted data from a string.
* **Arguments:**
    * `string`: The string to read from.
    * `REWIND`: An optional keyword to rewind the string pointer to the beginning of the string before reading.
    * `type`: The type of data to read (e.g., `s` for string, `n` for number).
    * `var`: The variable to store the read data in.

### `ReplicateConstraint`

* **Syntax:** `ReplicateConstraint(<constraint pattern in terms of 'this1', 'this2',...>, <an argument to replace 'this*', for each 'thisN' in the pattern>);`
* **Description:** Replicates a constraint, which is a way to enforce relationships between parameters in a model.
* **Arguments:**
    * `constraint pattern`: The pattern for the constraint.
    * `argument`: The arguments to replace the placeholders in the pattern.

### `category`

* **Syntax:** `category <name> = <definition>;`
* **Description:** Defines a new category, which can be used to group sites with similar evolutionary characteristics.
* **Arguments:**
    * `name`: The name of the category.
    * `definition`: The definition of the category.

### `Model`

* **Syntax:** `Model <name> = <definition>;`
* **Description:** Defines a new substitution model.
* **Arguments:**
    * `name`: The name of the model.
    * `definition`: The definition of the model.

### `ChoiceList`

* **Syntax:** `ChoiceList(store_here {ID}, title {String}, how many choices (0 for any number >= 1) {Integer}, [NO_SKIP or SKIP_NONE | list of indices not to show as options], [source object | comma separated list of 'key', 'description' pairs])`
* **Description:** Creates a choice list for user interaction, allowing the user to select from a list of options.
* **Arguments:**
    * `store_here`: The variable to store the user's choice in.
    * `title`: The title of the choice list.
    * `how many choices`: The number of choices to allow.
    * `NO_SKIP or SKIP_NONE`: Options for skipping choices.
    * `source object`: The object to get the choices from.

### `GetInformation`

* **Syntax:** `GetInformation(<receptacle>, <DataSet or DataSetFilter or LikelihoodFunction or Model or Variable or Regexp or String>)`
* **Description:** Retrieves information about an object, such as a dataset, model, or variable.
* **Arguments:**
    * `receptacle`: The variable to store the information in.
    * `object`: The object to get information about.

### `ExecuteCommands`

* **Syntax:** `ExecuteCommands(<source code>, [optional <'compiled' | (input redirect , [optional <namespace>]) ])`
* **Description:** Executes a string of HBL commands.
* **Arguments:**
    * `source code`: The string of HBL commands to execute.
    * `compiled`: An optional keyword to indicate that the source code is pre-compiled.
    * `input redirect`: An optional argument to redirect input.
    * `namespace`: An optional argument to specify a namespace.

### `ExecuteAFile`

* **Syntax:** `ExecuteAFile(<file path>, [optional <'compiled' | (input redirect , [optional <namespace>]) ])`
* **Description:** Executes an HBL script from a file.
* **Arguments:**
    * `file path`: The path to the HBL file to execute.
    * `compiled`: An optional keyword to indicate that the source code is pre-compiled.
    * `input redirect`: An optional argument to redirect input.
    * `namespace`: An optional argument to specify a namespace.

### `LoadFunctionLibrary`

* **Syntax:** `LoadFunctionLibrary(<file path | library name>, [optional <'compiled' | (input redirect , [optional <namespace>]) ])`
* **Description:** Loads a function library, which is a collection of HBL functions.
* **Arguments:**
    * `file path | library name`: The path to the library file or the name of the library.
    * `compiled`: An optional keyword to indicate that the source code is pre-compiled.
    * `input redirect`: An optional argument to redirect input.
    * `namespace`: An optional argument to specify a namespace.

### `FindRoot`

* **Syntax:** `FindRoot (<receptacle>, <expression>, <variable to solve for>,<left bound>,<right bound>)`
* **Description:** Finds the root of an equation using a numerical solver.
* **Arguments:**
    * `receptacle`: The variable to store the root in.
    * `expression`: The expression to find the root of.
    * `variable to solve for`: The variable to solve for.
    * `left bound`: The left bound of the search interval.
    * `right bound`: The right bound of the search interval.

### `MPIReceive`

* **Syntax:** `MPIReceive (<from node; or -1 to receive from any>, <message storage>, <sender index storage>)`
* **Description:** Receives a message in an MPI (Message Passing Interface) environment, which is used for parallel computing.
* **Arguments:**
    * `from node`: The ID of the node to receive from, or -1 to receive from any node.
    * `message storage`: The variable to store the received message in.
    * `sender index storage`: The variable to store the sender's ID in.

### `MPISend`

* **Syntax:** `MPISend(<node id>, <string | likelihood function ID | filename [in conjunction with argument 3]>, [if specified, treat the second argument as a script path, and use the dict supplied here as input options to the script])`
* **Description:** Sends a message in an MPI environment.
* **Arguments:**
    * `node id`: The ID of the node to send to.
    * `string | likelihood function ID | filename`: The message to send.
    * `dict`: An optional dictionary of input options.

### `GetDataInfo`

* **Syntax:** `GetDataInfo(<receptacle>, <DataSet or DataSetFilter>, [optional <sequence ref, site ref | sequence 1 , sequence 2, DISTANCES>])`
* **Description:** Retrieves information about a dataset, such as the number of sequences and sites.
* **Arguments:**
    * `receptacle`: The variable to store the information in.
    * `DataSet or DataSetFilter`: The dataset or dataset filter to get information about.
    * `sequence ref, site ref`: Optional arguments to specify a subset of the data.

### `StateCounter`

* **Syntax:** `StateCounter(<receptacle>, <DataSet or DataSetFilter>, <state>)`
* **Description:** Counts the occurrences of a specific state (e.g., a nucleotide or amino acid) in a dataset.
* **Arguments:**
    * `receptacle`: The variable to store the count in.
    * `DataSet or DataSetFilter`: The dataset or dataset filter to count the state in.
    * `state`: The state to count.

### `Integrate`

* **Syntax:** `Integrate (<receptacle>, <expression>, <variable to integrate over for>,<left bound>,<right bound>)`
* **Description:** Performs numerical integration of an expression over a given interval.
* **Arguments:**
    * `receptacle`: The variable to store the result of the integration in.
    * `expression`: The expression to integrate.
    * `variable to integrate over for`: The variable to integrate with respect to.
    * `left bound`: The lower bound of the integration interval.
    * `right bound`: The upper bound of the integration interval.

### `DoSQL`

* **Syntax:** `DoSQL (<dbID | SQL_OPEN | SQL_CLOSE>, <transaction string | file name>, <ID here | result here>)`
* **Description:** Executes a SQL query on a database.
* **Arguments:**
    * `dbID | SQL_OPEN | SQL_CLOSE`: The database ID, or a command to open or close a database connection.
    * `transaction string | file name`: The SQL query to execute or a file containing the query.
    * `ID here | result here`: The variable to store the result in.

### `Topology`

* **Syntax:** `Topology <name> = <tree_string>;`
* **Description:** Defines a new tree topology from a Newick-formatted string.
* **Arguments:**
    * `name`: The name of the topology.
    * `tree_string`: The tree in Newick format.

### `AlignSequences`

* **Syntax:** `AlignSequences (result, sequences, options)`
* **Description:** Aligns a set of sequences using a specified alignment algorithm.
* **Arguments:**
    * `result`: The variable to store the alignment in.
    * `sequences`: The sequences to align.
    * `options`: The alignment options.

### `#profile`

* **Syntax:** `#profile START|PAUSE|RESUME|where to store`
* **Description:** Controls the code profiler, which is used to measure the performance of HBL code.
* **Arguments:**
    * `START|PAUSE|RESUME|where to store`: A command to start, pause, resume, or store the profiling results.

### `SCFG`

* **Syntax:** `SCFG ident = (Rules1, Rules2 <,start>)`
* **Description:** Defines a Stochastic Context-Free Grammar (SCFG), which is a probabilistic model used for sequence analysis.
* **Arguments:**
    * `ident`: The identifier for the SCFG.
    * `Rules1, Rules2`: The rules of the grammar.
    * `start`: The start symbol of the grammar.

### `NeuralNet`

* **Syntax:** `NeuralNet <name> = <definition>;`
* **Description:** Defines a neural network, which can be used for various machine learning tasks.
* **Arguments:**
    * `name`: The name of the neural network.
    * `definition`: The definition of the neural network.

### `BGM`

* **Syntax:** `BGM ident = (<nodes>)`
* **Description:** Defines a Bayesian Graphical Model (BGM), which is a probabilistic model that represents the dependencies between variables.
* **Arguments:**
    * `ident`: The identifier for the BGM.
    * `nodes`: The nodes of the graph.

### `SimulateDataSet`

* **Syntax:** `SimulateDataSet (<DataSetFilter>, <Tree>, <Model>)`
* **Description:** Simulates a dataset based on a tree, a model, and a dataset filter.
* **Arguments:**
    * `DataSetFilter`: The dataset filter to use.
    * `Tree`: The tree to simulate the data on.
    * `Model`: The substitution model to use.

### `KeywordArgument`

* **Syntax:** `KeywordArgument (keyword, description, [default value, [dialog reference]])`
* **Description:** Defines a keyword argument for a function, which allows you to pass arguments by name instead of by position.
* **Arguments:**
    * `keyword`: The name of the keyword.
    * `description`: A description of the argument.
    * `default value`: An optional default value for the argument.
    * `dialog reference`: An optional reference to a dialog.

### `ConvertBranchLength`

* **Syntax:** `ConvertBranchLength(<receptacle>, <expression>, <variable to solve for>, <target value>)`
* **Description:** Converts the branch lengths of a tree based on an expression.
* **Arguments:**
    * `receptacle`: The variable to store the new tree in.
    * `expression`: The expression to use for converting the branch lengths.
    * `variable to solve for`: The variable to solve for in the expression.
    * `target value`: The target value for the conversion.

### `for`

* **Syntax:** `for (<initialization>;<condition>;<increment>) {loop body}`
* **Description:** A standard for loop, which allows you to execute a block of code a specific number of times.
* **Arguments:**
    * `initialization`: An expression to be executed before the loop starts.
    * `condition`: An expression to be evaluated before each iteration of the loop.
    * `increment`: An expression to be executed after each iteration of the loop.

### `while`

* **Syntax:** `while (<condition>) {loop body}`
* **Description:** A standard while loop, which allows you to execute a block of code as long as a condition is true.
* **Arguments:**
    * `condition`: An expression to be evaluated before each iteration of the loop.

### `SetDialogPrompt`

* **Syntax:** `SetDialogPrompt(<prompt string>);`
* **Description:** Sets the prompt for a dialog, which is a window that can be used to get input from the user.
* **Arguments:**
    * `prompt string`: The string to be used as the prompt.

### `HarvestFrequencies`

* **Syntax:** `HarvestFrequencies(<receptacle>, <DataSet or DataSetFilter>, <atom INTEGER>, <unit INTEGER <= atom>, <position aware 0 or 1>, [optional site partition], [optional sequence partition] (only for DataSetArguments))`
* **Description:** Harvests frequencies of states (e.g., nucleotides or amino acids) from a dataset.
* **Arguments:**
    * `receptacle`: The variable to store the frequencies in.
    * `DataSet or DataSetFilter`: The dataset or dataset filter to harvest the frequencies from.
    * `atom`: The size of the atom to be used for harvesting.
    * `unit`: The size of the unit to be used for harvesting.
    * `position aware`: A flag to indicate whether the harvesting should be position-aware.
    * `site partition`: An optional argument to specify a site partition.
    * `sequence partition`: An optional argument to specify a sequence partition.

### `Optimize`

* **Syntax:** `Optimize (<receptacle>, <likelihood function/scfg/bgm>, [optional dictionary of arguments])`
* **Description:** Optimizes a likelihood function, SCFG, or BGM to find the parameter values that maximize the likelihood of the data.
* **Arguments:**
    * `receptacle`: The variable to store the optimization results in.
    * `likelihood function/scfg/bgm`: The object to be optimized.
    * `dictionary of arguments`: An optional dictionary of arguments for the optimizer.

### `LFCompute`

* **Syntax:** `LFCompute (<likelihood function/scfg/bgm>,<LF_START_COMPUTE|LF_DONE_COMPUTE|receptacle>)`
* **Description:** Computes a likelihood function.
* **Arguments:**
    * `likelihood function/scfg/bgm`: The object to be computed.
    * `LF_START_COMPUTE|LF_DONE_COMPUTE|receptacle`: A command to start or finish the computation, or a variable to store the result in.

### `CovarianceMatrix`

* **Syntax:** `CovarianceMatrix (<receptacle>, <likelihood function/scfg/bgm>)`
* **Description:** Computes the covariance matrix of a model, which is a measure of the uncertainty in the parameter estimates.
* **Arguments:**
    * `receptacle`: The variable to store the covariance matrix in.
    * `likelihood function/scfg/bgm`: The object to compute the covariance matrix for.

### `SelectTemplateModel`

* **Syntax:** `SelectTemplateModel(<DataSetFilter>);`
* **Description:** Selects a template model, which is a pre-defined model that can be used as a starting point for building a new model.
* **Arguments:**
    * `DataSetFilter`: The dataset filter to use for selecting the model.

### `UseModel`

* **Syntax:** `UseModel (<model ID>)`
* **Description:** Specifies which model to use for a particular analysis.
* **Arguments:**
    * `model ID`: The ID of the model to use.

### `SetParameter`

* **Syntax:** `SetParameter(<object>, <parameter index>, <value>)`
* **Description:** Sets the value of a parameter in a model or other object.
* **Arguments:**
    * `object`: The object to set the parameter in.
    * `parameter index`: The index of the parameter to set.
    * `value`: The new value for the parameter.

### `assert`

* **Syntax:** `assert (<statement>,[optional message on failure]>)`
* **Description:** Asserts that a statement is true. If the statement is false, the program will terminate with an error message.
* **Arguments:**
    * `statement`: The statement to be asserted.
    * `message on failure`: An optional message to be displayed if the assertion fails.

### `RequireVersion`

* **Syntax:** `RequireVersion (<version string>)`
* **Description:** Requires a specific version of HyPhy to be running. If the version is not correct, the program will terminate.
* **Arguments:**
    * `version string`: The required version of HyPhy.

### `DeleteObject`

* **Syntax:** `DeleteObject(<object 1> [optional ,<object 2>, <object 3>, ..., <object N>])`
* **Description:** Deletes one or more objects from memory.
* **Arguments:**
    * `object`: The object(s) to be deleted.

### `ClearConstraints`

* **Syntax:** `ClearConstraints(<object 1> [optional ,<object 2>, <object 3>, ..., <object N>])`
* **Description:** Clears constraints from one or more objects.
* **Arguments:**
    * `object`: The object(s) to clear the constraints from.

### `fprintf`

* **Syntax:** `fprintf(stdout|MESSAGE_LOG|TEMP_FILE_NAME|PROMPT_FOR_FILE|file path, object 1 [optional ,<object 2>, ..., <object N>])`
* **Description:** Prints formatted output to a file or stream.
* **Arguments:**
    * `stdout|MESSAGE_LOG|TEMP_FILE_NAME|PROMPT_FOR_FILE|file path`: The destination for the output.
    * `object`: The object(s) to be printed.

### `Export`

* **Syntax:** `Export (<string variable ID>, <object ID>, [export options as a dict])`
* **Description:** Exports an object to a string representation.
* **Arguments:**
    * `string variable ID`: The variable to store the exported string in.
    * `object ID`: The ID of the object to be exported.
    * `export options`: An optional dictionary of export options.

### `GetURL`

* **Syntax:** `GetURL (<receptacle>,<URL>[, SAVE_TO_FILE])`
* **Description:** Fetches data from a URL.
* **Arguments:**
    * `receptacle`: The variable to store the fetched data in.
    * `URL`: The URL to fetch the data from.
    * `SAVE_TO_FILE`: An optional keyword to save the data to a file.

### `GetString`

* **Syntax:** `GetString(<receptacle>,<object>,<index>,[optional <second index>])`
* **Description:** Retrieves a string from an object.
* **Arguments:**
    * `receptacle`: The variable to store the string in.
    * `object`: The object to retrieve the string from.
    * `index`: The index of the string to retrieve.
    * `second index`: An optional second index for multi-dimensional objects.

### `Differentiate`

* **Syntax:** `Differentiate(<receptacle>, <the expression to differentiate>, <variable to differentiate>[, number of times, default = 1])`
* **Description:** Differentiates an expression with respect to a variable.
* **Arguments:**
    * `receptacle`: The variable to store the result in.
    * `the expression to differentiate`: The expression to be differentiated.
    * `variable to differentiate`: The variable to differentiate with respect to.
    * `number of times`: The number of times to differentiate.