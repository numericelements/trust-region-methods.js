// www.numericelements.com
// <nu>meric elemen<ts> = nuts
/*global nuts*/


nuts.CholeskyDecomposition = function (symmetricMatrix) {
    'use strict';

    // see Matrix Computation, Golub and Van Loan, p. 164

    this.g = symmetricMatrix.squareMatrix();
    this.success = false;
    this.CLOSE_TO_ZERO = 10e-8;
    this.firstNonPositiveDefiniteLeadingSubmatrixSize = undefined;
    //this.firstNonPositiveDefiniteLeadingSubmatrixSize = 1;


    if (this.g.get(0, 0) < this.CLOSE_TO_ZERO) {
        return;
    }

    var sqrtGjj = Math.sqrt(this.g.get(0, 0)),
        i,
        j,
        k,
        sum;
    for (i = 0; i < this.g.size; i += 1) {
        this.g.divideAt(i, 0, sqrtGjj);
    }

    for (j = 1; j < this.g.size; j += 1) {
        for (i = j; i < this.g.size; i += 1) {
            sum = 0;
            for (k = 0; k < j; k += 1) {
                sum += this.g.get(i, k) * this.g.get(j, k);
            }
            this.g.substractAt(i, j, sum);
        }
        if (this.g.get(j, j) < this.CLOSE_TO_ZERO) {
            this.firstNonPositiveDefiniteLeadingSubmatrixSize = j + 1;
            return;
        }
        sqrtGjj = Math.sqrt(this.g.get(j, j));
        for (i = j; i < this.g.size; i += 1) {
            this.g.divideAt(i, j, sqrtGjj);
        }

    }

    for (j = 0; j < this.g.size; j += 1) {
        for (i = 0; i < j; i += 1) {
            this.g.set(i, j, 0);
        }
    }

    this.success = true;


};


nuts.CholeskyDecomposition.prototype = {
    solve : function (b) {
        'use strict';
        // See Numerical Recipes Third Edition p. 101

        //if (!this.success) {
        //    throw "CholeskyDecomposistion.success === false";
        //}

        var n = this.g.size,
            x = b.slice(),
            i,
            k,
            sum;
        // Ly = b
        for (i = 0; i < n; i += 1) {
            sum = b[i];
            for (k = i - 1; k >= 0; k -= 1) {
                sum -= this.g.get(i, k) * x[k];
            }
            x[i] = sum / this.g.get(i, i);
        }
        // LT x = Y
        for (i = n - 1; i >= 0; i -= 1) {
            sum = x[i];
            for (k = i + 1; k < n; k += 1) {
                sum -= this.g.get(k, i) * x[k];
            }
            x[i] = sum / this.g.get(i, i);
        }

        return x;
    },

    solve_LT_result_equal_b : function (b) {
        'use strict';
        var n = this.g.size,
            x = b.slice(),
            i,
            k,
            sum;
        for (i = 0; i < n; i += 1) {
            sum = b[i];
            for (k = i - 1; k >= 0; k -= 1) {
                sum -= this.g.get(i, k) * x[k];
            }
            x[i] = sum / this.g.get(i, i);
        }
        return x;
    },

    getG : function () {
        'use strict';
        return this.g;
    }
};
