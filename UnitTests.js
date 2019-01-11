// www.numericelements.com
// <nu>meric elemen<ts> = nuts
/*global nuts*/
/*global console*/

var nuts = {};

nuts.TrustRegionSubproblem_Test = function () {
    'use strict';
    var hessian = new nuts.SymmetricMatrix(2),
        gradient = [0, 1],
        trustRegionSubproblem,
        cauchyPoint,
        trustRegionResult,
        theta,
        r = new nuts.SquareMatrix(2),
        rt = new nuts.SquareMatrix(2),
        mr,
        mrs;


    hessian.set(0, 0, -2);
    hessian.set(1, 1, 2);
    trustRegionSubproblem = new nuts.TrustRegionSubproblem(gradient, hessian);
    cauchyPoint = trustRegionSubproblem.cauchyPoint(1);


    if (Math.abs(cauchyPoint[0]) > 1e-8 || Math.abs(-0.5 - cauchyPoint[1]) > 1e-8) {
        console.log("UnitTests for TrustRegionSubproblem_Test, cauchyPoint failed");
        console.log(cauchyPoint);
    }

    gradient = [1, 0];
    trustRegionSubproblem = new nuts.TrustRegionSubproblem(gradient, hessian);
    trustRegionResult = trustRegionSubproblem.solve(1);

    if (Math.abs(-1 - trustRegionResult.step[0]) > 1e-8 || Math.abs(trustRegionResult.step[1]) > 1e-8) {
        console.log("UnitTests for TrustRegionSubproblem_Test, solve failed");
        console.log(trustRegionResult);
    }

    hessian.set(0, 0, 2);
    gradient = [0, 1];

    trustRegionSubproblem = new nuts.TrustRegionSubproblem(gradient, hessian);
    trustRegionResult = trustRegionSubproblem.solve(1);

    if (Math.abs(trustRegionResult.step[0]) > 1e-8 || Math.abs(-0.5 - trustRegionResult.step[1]) > 1e-8) {
        console.log("UnitTests for TrustRegionSubproblem_Test, solve failed");
        console.log(trustRegionResult);
    }

    hessian.set(0, 0, 2);
    hessian.set(1, 1, -0.5);
    hessian.set(1, 0, 0);
    theta = 45 / 180 * Math.PI;
    r.set(0, 0, Math.cos(theta));
    rt.set(0, 0, Math.cos(theta));
    r.set(1, 1, Math.cos(theta));
    rt.set(1, 1, Math.cos(theta));
    r.set(0, 1, -Math.sin(theta));
    rt.set(1, 0, -Math.sin(theta));
    r.set(1, 0, Math.sin(theta));
    rt.set(0, 1, Math.sin(theta));

    mr = r.multiplyByMatrix(hessian).multiplyByMatrix(rt);
    mrs = new nuts.SymmetricMatrix(2);
    mrs.set(0, 0, mr.get(0, 0));
    mrs.set(1, 0, mr.get(1, 0));
    mrs.set(1, 1, mr.get(1, 1));

    gradient = [0, 0];
    trustRegionSubproblem = new nuts.TrustRegionSubproblem(gradient, mrs);
    trustRegionResult = trustRegionSubproblem.solve(1);

    if (Math.abs(trustRegionResult.step[0] - Math.sin(theta)) > 1e-8 || Math.abs(trustRegionResult.step[1] + Math.cos(theta)) > 1e-8) {
        console.log("UnitTests for TrustRegionSubproblem_Test, solve failed");
        console.log(trustRegionResult);
    }



};

nuts.CholeskyDecomposition_Test = function () {
    'use strict';
    var matrix,
        choleskyDecomposition,
        b = [2, 3],
        b_verify,
        a,
        i,
        g;
    matrix = new nuts.SymmetricMatrix(2);
    matrix.set(0, 0, 1);
    matrix.set(1, 0, 0);
    matrix.set(1, 1, 1);
    choleskyDecomposition = new nuts.CholeskyDecomposition(matrix);

    if (choleskyDecomposition.success === false) {
        console.log("UnitTests for CholeskyDecomposition_Test, constructor failed");
        console.log(choleskyDecomposition);
    }

    a = choleskyDecomposition.solve(b);

    if (a[0] !== b[0] || a[1] !== b[1]) {
        console.log("UnitTests for CholeskyDecomposition_Test, solve failed");
        console.log(a);
    }

    matrix.set(1, 1, -1);

    choleskyDecomposition = new nuts.CholeskyDecomposition(matrix);

    if (choleskyDecomposition.success === true) {
        console.log("UnitTests for CholeskyDecomposition_Test, constructor failed");
        console.log(choleskyDecomposition);
    }

    matrix = new nuts.SymmetricMatrix(3);
    matrix.set(0, 0, 4);
    matrix.set(1, 0, 12);
    matrix.set(1, 1, 37);
    matrix.set(2, 0, -16);
    matrix.set(2, 1, -43);
    matrix.set(2, 2, 98);

    choleskyDecomposition = new nuts.CholeskyDecomposition(matrix);

    if (choleskyDecomposition.success === false) {
        console.log("UnitTests for CholeskyDecomposition_Test, constructor failed");
        console.log(choleskyDecomposition);
    }

    b = [3, 4, 6];
    a = choleskyDecomposition.solve(b);
    b_verify = matrix.multiplyByVector(a);

    g = choleskyDecomposition.getG();

    if (g.get(0, 0) !== 2 || g.get(1, 0) !== 6  || g.get(2, 0) !== -8 || g.get(1, 1) !== 1 || g.get(2, 1) !== 5 || g.get(2, 2) !== 3) {
        console.log("UnitTests for CholeskyDecomposition_Test, solve failed");
        console.log(g.get(0, 0));
        console.log(g.get(1, 0));
        console.log(g.get(2, 0));
        console.log(g.get(1, 1));
        console.log(g.get(2, 1));
        console.log(g.get(2, 2));
    }



    for (i = 0; i < 3; i += 1) {
        if (Math.abs(b_verify[i] - b[i]) > 10e-8) {
            console.log("UnitTests for CholeskyDecomposition_Test, solve failed");
            //console.log(choleskyDecomposition);
            console.log(b);
            console.log(b_verify);
            console.log(a);
            console.log(matrix);
            console.log(choleskyDecomposition);
            break;
        }
    }
};




nuts.runAllUnitTests = function () {
    'use strict';
    nuts.TrustRegionSubproblem_Test();
    nuts.CholeskyDecomposition_Test();
    console.log("Run of all unit tests is complete")
};
