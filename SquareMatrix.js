// www.numericelements.com
// <nu>meric elemen<ts> = nuts
/*global nuts*/


nuts.SquareMatrix = function (size, data) {
    'use strict';
    this.size = size;
    this.data = [];
    var i;
    if (data === undefined) {
        this.data = [];
        for (i = 0; i < this.size * this.size; i += 1) {
            this.data.push(0);
        }
    } else {
        this.data = data.slice();
    }
};

nuts.SquareMatrix.prototype =  {
    get : function (row, column) {
        'use strict';
        return this.data[row * this.size + column];
    },

    set : function (row, column, value) {
        'use strict';
        this.data[row * this.size + column] = value;
    },

    substractAt : function (row, column, value) {
        'use strict';
        this.data[row * this.size + column] -= value;
    },

    multiplyAt : function (row, column, value) {
        'use strict';
        this.data[row * this.size + column] *= value;
    },

    divideAt : function (row, column, value) {
        'use strict';
        this.data[row * this.size + column] /= value;
    },

    multiplyByMatrix : function (that) {
        'use strict';
        var result = new nuts.SquareMatrix(this.size),
            i,
            j,
            k;
        for (i = 0; i < this.size; i += 1) {
            for (j = 0; j < this.size; j += 1) {
                for (k = 0; k < this.size; k += 1) {
                    result.set(i, j, this.get(i, k) * that.get(k, j));
                }
            }
        }
        return result;
    }

};
