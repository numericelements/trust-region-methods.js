// www.numericelements.com
// <nu>meric elemen<ts> = nuts
/*global nuts*/


nuts.SymmetricMatrix = function (size, data) {
    'use strict';
    this.size = size;

    if (data === undefined) {
        this.data = [];
        var i;
        for (i = 0; i < (this.size * (this.size + 1)) / 2; i += 1) {
            this.data.push(0);
        }
    } else {
        this.data = data.slice();
    }
};

nuts.SymmetricMatrix.prototype =  {

    get : function (row, column) {
        'use strict';
        if (row <= column) {
            return this.data[row * this.size - (row - 1) * row / 2 + column - row];
        }
        return this.data[column * this.size - (column - 1) * column / 2 + row - column];

    },

    set : function (row, column, value) {
        'use strict';
        if (row <= column) {
            this.data[row * this.size - (row - 1) * row / 2 + column - row] = value;
        } else {
            this.data[column * this.size - (column - 1) * column / 2 + row - column] = value;
        }
    },

    addAt : function (row, column, value) {
        'use strict';
        if (row <= column) {
            this.data[row * this.size - (row - 1) * row / 2 + column - row] += value;
        } else {
            this.data[column * this.size - (column - 1) * column / 2 + row - column] += value;
        }

    },

    addValueOnDiagonal : function (value) {
        'use strict';
        var i;
        for (i = 0; i < this.size; i += 1) {
            this.addAt(i, i, value);
        }

    },

    add : function (thatSymmetricMatrix) {
        'use strict';
        var result = new nuts.SymmetricMatrix(this.size),
            i,
            j;
        for (i = 0; i < this.size; i += 1) {
            for (j = 0; j <= i; j += 1) {
                result.set(i, j, this.get(i, j) + thatSymmetricMatrix.get(i, j));
            }
        }
        return result;
    },

    clone : function () {
        'use strict';
        return new nuts.SymmetricMatrix(this.size, this.data);
    },

    squareMatrix : function () {
        'use strict';
        var i,
            j,
            result = new nuts.SquareMatrix(this.size);
        for (i = 0; i < this.size; i += 1) {
            for (j = 0; j < this.size; j += 1) {
                result.set(i, j, this.get(i, j));
            }
        }
        return result;
    },

    multiplyByVector : function (v) {
        'use strict';
        var result = [],
            i,
            j,
            temp;
        for (i = 0; i < this.size; i += 1) {
            temp = 0;
            for (j = 0; j < this.size; j += 1) {
                temp += this.get(i, j) * v[j];
            }
            result.push(temp);
        }
        return result;

    },

    multiplyByScalar : function (value) {
        'use strict';
        var result = new nuts.SymmetricMatrix(this.size),
            i,
            j;
        for (i = 0; i < this.size; i += 1) {
            for (j = 0; j <= i; j += 1) {
                result.set(i, j, this.get(i, j) * value);
            }
        }
        return result;
    },

    quadraticForm : function (v) {
        'use strict';
        var result = 0,
            i,
            j;
        for (i = 1; i < this.size; i += 1) {
            for (j = 0; j < i; j += 1) {
                result += this.get(i, j) * v[i] * v[j];
            }
        }
        result *= 2;
        for (i = 0; i < this.size; i += 1) {
            result += this.get(i, i) * Math.pow(v[i], 2);
        }
        return result;
    },

    numRows : function () {
        'use strict';
        return this.size;
    },

    numCols : function () {
        'use strict';
        return this.size;
    }

};
