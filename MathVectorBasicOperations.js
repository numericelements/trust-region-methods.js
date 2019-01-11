// www.numericelements.com
// <nu>meric elemen<ts> = nuts
/*global nuts*/



nuts.multiplyVectorByScalar = function (vector, value) {
    'use strict';
    var result = [],
        i;
    for (i = 0; i < vector.length; i += 1) {
        result.push(vector[i] * value);
    }
    return result;
};

nuts.divideVectorByScalar = function (vector, value) {
    'use strict';
    var result = [],
        i;
    for (i = 0; i < vector.length; i += 1) {
        result.push(vector[i] / value);
    }
    return result;
};

nuts.saxpy2 = function (a, x, y) {
    'use strict';
    if (x.length !== y.length) {
        throw "nut.saxpy2, adding two vectors of different length";
    }
    var result = [],
        i;
    for (i = 0; i < x.length; i += 1) {
        result.push(a * x[i] + y[i]);
    }
    return result;
};

nuts.saxpy = function (a, x, y) {
    'use strict';
    // y = y + ax;
    // See Matrix Computations 4th edition, Golub and Van Loan, p. 4
    if (x.length !== y.length) {
        throw "nut.saxpy, adding two vectors of different length";
    }
    var i;
    for (i = 0; i < x.length; i += 1) {
        y[i] += a * x[i];
    }
};

nuts.dotProduct = function (x, y) {
    'use strict';
    var result = 0,
        i;
    for (i = 0; i < x.length; i += 1) {
        result += x[i] * y[i];
    }
    return result;
};

nuts.addTwoVectors = function (x, y) {
    'use strict';
    if (x.length !== y.length) {
        throw "nut.addTwoVectors, adding two vectors of different length";
    }
    var result = [],
        i;
    for (i = 0; i < x.length; i += 1) {
        result.push(x[i] + y[i]);
    }
    return result;
};

nuts.squaredNorm = function (v) {
    'use strict';
    var result = 0,
        i;
    for (i = 0; i < v.length; i += 1) {
        result += v[i] * v[i];
    }
    return result;
};

nuts.norm = function (v) {
    'use strict';
    return Math.sqrt(nuts.squaredNorm(v));
};

nuts.zeroVector = function (n) {
    'use strict';
    var result = [],
        i;
    for (i = 0; i < n; i += 1) {
        result.push(0);
    }
    return result;
};
