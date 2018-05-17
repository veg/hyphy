
/**
 * The R-8 method https://en.wikipedia.org/wiki/Quantile#Estimating_quantiles_from_a_sample
 * @name stats.Quantile
 * @param {Matrix} column - the vector to operate on
 * @param {Number} p
 * @returns {Number} the quantile
 */
lfunction stats.Quantile(column, p) {

    N = Rows(column);

    h = (N + 1 / 3) * p - 2 / 3;

    if (h < 0) {
        return column[0];
    }

    if (h >= N - 1) {
        return column[N - 1];
    }

    h_floor = h$1;

    return column[h_floor] + (h - h_floor) * (column[h_floor + 1] - column[h_floor]);
}

/**
 * Compute a simple Bayes factor given priors and posteriors
 * @name stats.BayesFactor
 * @param {Number} prior - prior for event E
 * @param {Number} posterior - posterior for event E
 * @returns {Number} the BF
 */

lfunction stats.BayesFactor(prior, posterior) {
    if (prior == 1) { // infinite prior odds
        if (posterior == 1) {
            return 1;
        } 
        return 0;
    }
    return posterior * (1-prior) / (1-posterior) / prior;
}
