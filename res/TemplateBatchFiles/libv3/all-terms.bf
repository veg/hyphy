/* Terms used throughout HyPhy */

// These terms must be defined in this manner in order to avoid conflicts w/ built-in functions global and Gamma
terms.global               = "global";
terms.json.global          = "Global model fit"; // TODO: Kill
terms.rate_variation.Gamma = "Gamma";


namespace terms{

    /* Generic terms which are used in a variety of contexts */
    alphabet                = "alphabet";
    code                    = "code";
    bases                   = "bases";
    amino_acid              = "amino-acid";
    codons                  = "codons";
    codon                   = "codon";
    sense_codons            = "sense";
    nucleotide              = "nucleotide";
    dinucleotide            = "dinucleotide";
    binary                  = "binary";
    stop_codons             = "stop";
    terminal_stop           = "terminal_stop";
    translation_table       = "translation-table";
    synonymous_sub_count    = "synonymous substitution count";
    nonsynonymous_sub_count = "nonsynonymous substitution count";
    original_name           = "original name";
    replicates              = "replicates";
    data_type               = "datatype";
    devnull                 = "/dev/null";
    _namespace              = "namespace";

    category            = "category";
    mixture             = "mixture";
    //timers              = "timers";
    constraint          = "constraint";
    fix                 = "fix-me";
    likelihood_function = "LF";
    model               = "model";
    number_precision    = "number-precision";
    three_way           = "three-way";
    reduced             = "reduced";
    substitutions       = "substitutions";
    search_grid         = "search_grid";
    search_restarts     = "search_restarts";

    parameters     = "parameters";
    local          = "local";
   // global         = "global"; // Defined at the top of file
    globals_only   = "globals only";
    initial_values = "Initial values";
    default        = "default";
    id             = "ID"; // everything now
    description    = "description";

    efv_estimate   = "EFV";
    branch_length  = "branch length";
    alternative    = "alternative";
    Null           = "null";
    LRT            = "LRT";
    p_value        = "p-value";

    branch_selection_attributes   = "Branch selection attributes";
    empirical_bayes_factor  = "BF";
    posterior               = "posterior";

    global_mg94xrev = "Global MG94xREV";


    lower_bound = "LB";
    upper_bound = "UB";
    range01 = {
        lower_bound: "0",
        upper_bound: "1"
    };

    range_almost_01 = {
        lower_bound: "1e-8",
        upper_bound: "1"
    };

    range_gte1 = {
        lower_bound: "1",
        upper_bound: "1e10"
    };

    range_any = {
        lower_bound: "0",
        upper_bound: "1e25"
    };

    range_positive = {
        lower_bound: "0",
        upper_bound: "1e5"
    };

    range_clamp_locals = {
        lower_bound: "0",
        upper_bound: "100"
    };






    /* Term functions */
    function characterFrequency (character) {
        return "Equilibrium frequency for " + character;
    }
    function nucleotideRate(fromC, toC) {
        return "Substitution rate from nucleotide " + fromC + " to nucleotide " + toC;
    }

    function nucleotideRateReversible (fromC, toC) {
        if (fromC < toC) {
            return nucleotideRate (fromC, toC);
        }
        return nucleotideRate (toC, fromC);
    }

    function aminoacidRate(fromA, toA) {
        return "Substitution rate from amino-acid " + fromA + " to amino-acid " + toA;
    }

    function propertyImportance  (id) {
        return "Importance factor for property " + id;
    }

    function binaryRate(fromX, toX) {
        return "Substitution rate from character " + fromX + " to character " + toX;
    }
    function timeParameter() {
        return "Evolutionary time parameter";
    }
    function AddCategory (term, categoryID) {
        return term + " for category " + categoryID;
    }


    /* Terms accompanying category defintions */
    //category = "category" // Defined above under genetic terms. Left here as a comment for future reminders.
    namespace category {
        category_parameters = "category parameters";
        bins                = "bins";
        weights             = "weights";
        represent           = "represent";
        PDF                 = "PDF";
        CDF                 = "CDF";
        dCDF                = "dCDF";
        HMM                 = "HMM";
    }


    /* Terms for code structures */
    //code = "code";// Defined above under general terms. Left here as a comment for future reminders.
    namespace code{
        stops    = "stops";
        ordering = "ordering";
        mapping  = "mapping";
    }


    /* Terms associated with data structures */
    namespace data {
        sites          = "sites";
        sequences      = "sequences";
        sequence       = "sequence";
        sample_size    = "sample size";
        composition    = "composition";
        file           = "file";
        cache          = "cache";
        name           = "name";
        name_mapping   = "name-mapping";
        partitions     = "partitions";
        tree           = "tree";
        dataset        = "dataset";
        datafilter     = "datafilter";
        coverage       = "coverage";
        filter_string  = "filter-string";
        is_constant    = "is_constant";
        pattern_id     = "pattern id";
        filename_to_index = "filename-to-index";
    }


    /* Terms for evaluating diff's b/w states */
    namespace diff{
        to       = "to";
        from     = "from";
        position = "position";
    }


    /* Terms specific to model fits.*/
    // Note that some terms used in model fits are also under the genetic terms
    namespace fit{
        MLE                 = "MLE";
        trees               = "Trees";
        nonsynonymous_trees = "non-synonymous-trees";
        synonymous_trees    = "synonymous-trees";
        log_likelihood      = "LogL";
        filters             = "Filters";
        phase               = "phase";
    }


    /* Terms accompanying models/frequencies.bf */
    namespace frequencies {
        _4x1       = "Nucleotide 4x1 estimator";
        equal      = "Equal frequencies";
        CF3x4      = "Corrected 3x4 frequency estimator";
        F3x4      = "Standard F3x4 codon frequency estimator";
        F1x4      = "Standard F1x4 codon frequency estimator";
        _20x1      = "Protein 20x1 estimator";
        MLE        = "Maximum likelihood frequency estimator";
        predefined = "Based on a training set";
        binary     = "Binary character frequency estimator";
        run_time   = "Run time frequency estimator via matrix inversion"
    }

    /* Terms accompanying tasks/genetic_code.bf */
    namespace genetic_code {
        synonymous        = "synonymous";
        nonsynonymous     = "nonsynonymous";
        weighting_matrix  = "weighting-matrix";
        count_stop_codons = "count-stop-codons";
        EPS               = "EPS";
        EPN               = "EPN";
        OPS               = "OPS";
        OPN               = "OPN";
        NTP               = "NTP";
        SS                = "SS";
        NS                = "NS";
    }

    /* Terms used in I/O */
    namespace io {
        //Analysis banner
        info         = "info";
        requirements = "requirements";
        reference    = "citation";
        authors      = "authors";
        contact      = "contact";
        version      = "version";
        help         = "help";

        //ReadDelimitedFile
        header       = "header";
        rows         = "rows";
    }

    /* Terms used to write to JSONs */
    namespace json{
        // For any json
        json                  = "json";
        analysis              = "analysis";
        input                 = "input";
        file                  = "file name";
        sequences             = "number of sequences";
        sites                 = "number of sites";
        fits                  = "fits";
        timers                = "timers";
        trees                 = "trees";
        MLE                   = "MLE";
        parameters            = "estimated parameters";
        PMID                  = "PMID";
 //       PMCID                 = "PMCID";
        test_results          = "test results";
        tree_string           = "tree";
        tree_length           = "tree length";
        rate_distribution     = "Rate Distributions";
        log_likelihood        = "Log Likelihood";
        AICc                  = "AIC-c";
        global_mg94xrev       = "Global MG94xREV";
        mg94xrev_sep_rates    = "MG94xREV with separate rates for branch sets";
        nucleotide_gtr        = "Nucleotide GTR";
        baseline              = "Baseline";
        frequencies           = "Equilibrium frequencies";
        model                 = "model"; // TODO: change string to "model name"
       // global              = "Global model fit"; // Defined at the top of file
        attribute             = "attributes";
        display_order         = "display order";
        attribute_type        = "attribute type";
        node_label           = "node label";
        branch_label          = "branch label";
        branch_attributes     = "branch attributes";
        branch_annotations    = "branch annotations";
        annotation_tag        = "annotation tag";
        branch_lengths        = "branch lengths";
  //      background            = "background";

        headers               = "headers";
        content               = "content";
        partition_count       = "partition count";
        partitions            = "data partitions";

        tested                = "tested";
        uncorrected_pvalue    = "Uncorrected P-value";
        corrected_pvalue      = "Corrected P-value";
        pvalue_threshold      = "P-value threshold";
        relative_site_rates   = "Relative site rate estimates";
      //  site_log_likelihood   = "site log likelihood";
     //   evidence_ratios       = "evidence ratios";
        options               = "options";
        runtime               = "runtime";
        version               = "version";
        convergence_failures  = "convergence failures";
        omega_ratio           = "omega";
        rate                  = "rate";
        proportion            = "proportion";
        positive              = "positive test results";
    }


    /* Terms accompanying convenience/math.bf */
    namespace math {
        count        = "Count";
        mean         = "Mean";
        median       = "Median";
        stddev       = "Std.Dev";
        minimum      = "Min";
        maximum      = "Max";
        _2.5         = "2.5%";
        _97.5        = "97.5%";
        sum          = "Sum";
        variance     = "Variance";
        cov          = "COV";
        skewness     = "Skewness";
        kurtosis     = "Kurtosis";
        square_sum   = "Sq. sum";
        non_negative = "Non-negative";
    }



    /* Terms for mixture models */
    //mixture  = "mixture"; // Defined above under general terms. Left here as a comment for future reminders.
    namespace mixture {
        mixture_components = "mixture components";
        mixture_weight     = "Mixture weight";
        mixture_aux_weight = "Mixture auxiliary weight";
    }


    /* Terms accompanying models/* ;Primary terms used in model definitions */
    namespace model {

        efv_estimate_name       = "Equilibrium frequency estimator";
        frequency_estimator     = "frequency-estimator";
        efv_matrix              = "pi";
        efv_id                  = "efv-id";

        defineQ                 = "defineQ";
        q_ij                    = "q_ij";
        rate_matrix             = "Q";
        matrix_id               = "matrix-id";
        rate_entry              = "rate-entry";
        empirical_rates         = "empirical-rates";

        rate_variation          = "Rate variation";

        time                    = "time";

        type                    = "type";
        components              = "components";
        canonical               = "canonical";
        reversible              = "reversible";
        empirical               = "empirical";

        branch_length_constrain = "branch length constrain";// TODO
        get_branch_length       = "get-branch-length";
        set_branch_length       = "set-branch-length";
        constrain_branch_length = "constrain-branch-length";
        branch_length_string    = "branch-length-string";
        branch_length_string_conditional = "branch-length-string-conditional";
        branch_length_scaler    = "branch length scaler";
        post_definition         = "post-definition";
        length                  = "length";
        length_parameter        = "length parameter";
        length_expression       = "length-expression";

        data                    = "data";
        residue_properties      = "residue_properties";

    }



    /* Terms accompanying tasks/mpi.bf */
    namespace mpi {
        job_id              = "job_id";
        callback            = "callback";
        arguments           = "arguments";
        Models              = "Models";
        LikelihoodFunctions = "LikelihoodFunctions";
        Headers             = "Headers";
        Variables           = "Variables";
        Functions           = "Functions";
        DataSetFilters      = "DataSetFilters";
    }



    /* Terms accompanying models/parameters.bf or any terms specifically referring to a model parameter */
    namespace parameters {
        local_constrained             = "Local Constrained";// from terms.lf.local.constrained
        global_constrained            = "Global Constrained";// from terms.lf.global.constrained
        local_independent             = "Local Independent";
        global_independent            = "Global Independent";
        transition                    = "transition";
        transversion                  = "transversion";
        transition_transversion_ratio = "Transitition/transversion ratio";
        kappa                         = "kappa";
        synonymous_rate               = "synonymous rate";
        nonsynonymous_rate            = "non-synonymous rate";
        omega_ratio                   = "non-synonymous/synonymous rate ratio";
        log_omega_ratio               = "log of non-synonymous/synonymous rate ratio";
        multiple_hit_rate             = "rate at which 2 nucleotides are changed instantly within a single codon";
        triple_hit_rate               = "rate at which 3 nucleotides are changed instantly within a single codon";
        triple_hit_rate_syn           = "rate at which 3 nucleotides are changed instantly within a single codon between synonymous codon islands";

        one                           = "1";
        theta                         = "theta";
        default_time                  = "t";
        omegas                        = "omegas";
        omega                         = "omega";
        delta                         = "delta";
        weights                       = "weights";
        weight                        = "weight";
        rates                         = "rates";
        freqs                         = "f";
    }


    /* Terms for an AssociativeList prefix */
    prefix = "prefix";
    settings = "settings";
    namespace settings {
        branch_selector = "branch-selector";
    }


    /* Terms accompanying rate_variation.bf rate variation */
    //Note that rate_variation itself is a **model** term
    namespace rate_variation {
        distribution  = "distribution";
        options       = "options";
        rate_modifier = "rate_modifier";
        bins          = "Rate variation bins";
        //Gamma         = "Gamma"; // Defined at the top of file
        GammaI        = "Gamma+I";
        gamma_alpha   = "Shape parameter for the gamma distribution (alpha)";
        gamma_beta    = "Variance parameter for the gamma distribution (beta)";
        gamma_p_inv   = "Estimated proportion of invariant sites";
        hmm_lambda    = "HMM rate switching parameter";
        before        = "before";
        after         = "after";
        HMM           = "HMM";
    }



    /* Terms used to specify runtime options for model fitting */
    namespace run_options {
        model_type                        = "model-type";
        proportional_branch_length_scaler = "proportional-branch-length-scaler";
        retain_lf_object                  = "retain-lf-object";
        retain_model_object               = "retain-model-object";
        partitioned_omega                 = "partitioned-omega";
        apply_user_constraints            = "apply-user-constraints";
        optimization_settings             = "optimization-settings";
    }


    /* Terms formatting table output */
    namespace table_options{
        header               = "header";
        minimum_column_width = "min-column-width";
        align                = "align";
        column_widths        = "column-widths";
    }

    /* Terms used for runtime tracking */
    //timers = "timers";// Defined above under genetic terms. Left here as a comment for future reminders.
    namespace timers {
        timer = "timer";
        order = "order";
    }


    /* Terms associated with tree structures */
    namespace trees {
        newick              = "string";
        newick_with_lengths = "string_with_lengths";
        newick_annotated    = "annotated_string";
        model_map           = "model_map";
        partitioned         = "partitioned";
        model_list          = "model_list";
        rooted              = "rooted";
        root                = "root";
        branches            = "branches";
        meta                = "meta";
        comment             = "comment";
        neighbor_joining    = "neighbor-joining";
        data_for_neighbor_joining
                            = "FILTER_FOR_NEIGHBOR_JOINING";

        //node_name = "Name";
        //children = "Children";
        //parent = "Parent";
    }

    /* Terms associated with tree labeling */
    namespace tree_attributes {
        internal    = "internal";
        leaf        = "leaf";
        test        = "test";
        background  = "background";
    }

    /* Terms associated with BGMs */

    namespace bgm {
        namespace node {
            id = "NodeID";
            type = "NodeType";
            max_parents = "MaxParents";
            prior_size = "PriorSize";
            levels = "NumLevels";
        }
    }

}
