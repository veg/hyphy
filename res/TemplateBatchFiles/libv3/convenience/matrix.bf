

/** @module matrix */

/** 
 * upper triangle gets copied to lower triangle
 * for non-square matrices, the largest square minor is symmetrized        
 * @name matrix.Symmetrize
 * @param {Matrix} mx - matrix to perform transform on
 * @returns nothing
 */
lfunction matrix.Symmetrize(mx) {

    copy_dim = Min(Rows(mx), Columns(mx));

    for (r = 0; r < copy_dim; r += 1) {
        for (c = r + 1; c < copy_dim; c += 1) {
            mx[c][r] = mx[r][c];
        }
    }

}
