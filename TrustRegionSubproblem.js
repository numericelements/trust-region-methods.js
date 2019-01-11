// www.numericelements.com
// <nu>meric elemen<ts> = nuts
/*global nuts*/


nuts.TrustRegionSubproblem = function (gradient, hessian, k_easy, k_hard) {
    'use strict';
    this.gradient = gradient;
    this.hessian = hessian;
    this.k_easy = k_easy || 0.01;
    this.k_hard = k_hard || 0.02;
    this.UPDATE_COEFF = 0.01;
    this.CLOSE_TO_ZERO = 10e-5;
};

nuts.TrustRegionSubproblem.prototype = {
    cauchyPoint : function (trustRegionRadius) {
        'use strict';
        // Cauchy Point : The minimizer along the steepest descent (-gradient) direction subject to trust-region bound.
        // See Numerical Optimizatoin, second edition, Nocedal and Wright, p. 71
        var gHg = this.hessian.quadraticForm(this.gradient),
            gNorm = nuts.norm(this.gradient),
            result = nuts.multiplyVectorByScalar(this.gradient, -trustRegionRadius / gNorm),
            tau;
        if (gHg <= 0) {
            return result;
        }
        tau = Math.pow(gNorm, 3) / trustRegionRadius / gHg;
        if (tau < 1) {
            return nuts.multiplyVectorByScalar(result, tau);
        }
        return result;
    },

    solve : function (trustRegionRadius) {
        'use strict';
        // Conn A., Gould N., Toint P., Trust Region Methods p. 193
        var gNorm = nuts.norm(this.gradient),
            p,
            q,
            pSquaredNorm,
            pNorm,
            hitsBoundary = true,
            lambda_new,
            lambda = this.initialLambdas(trustRegionRadius),
            hessianPlusLambda,
            choleskyDecomposition,
            delta_lambda,
            s_min,
            intersection,
            stepLength,
            relative_error,
            quadraticTerm,
            sls,
            vNorm;
            //nit = 0;



        while (true) {

            /*
            nit += 1;
            if (nit > 1000) {
                break;
            }
            if (nit > 100) {
                console.log(100);
            }
            */



            hessianPlusLambda = this.hessian.clone();
            hessianPlusLambda.addValueOnDiagonal(lambda.current);
            choleskyDecomposition = new nuts.CholeskyDecomposition(hessianPlusLambda);

            if (choleskyDecomposition.success && gNorm > this.CLOSE_TO_ZERO) {
                p = choleskyDecomposition.solve(nuts.multiplyVectorByScalar(this.gradient, -1));
                pSquaredNorm = nuts.squaredNorm(p);
                pNorm = Math.sqrt(pSquaredNorm);

                // Check for interior convergence
                if (pNorm < trustRegionRadius && lambda.current === 0) {
                    hitsBoundary = false;
                    break;
                }

                q = choleskyDecomposition.solve_LT_result_equal_b(p);

                delta_lambda = pSquaredNorm / nuts.squaredNorm(q) * (pNorm - trustRegionRadius) / trustRegionRadius;
                lambda_new = lambda.current + delta_lambda;

                if (lambda_new < 0) {
                    lambda_new = (lambda.lb + lambda.ub) / 2;
                }

                if (pNorm < trustRegionRadius) { // Inside Boundary
                    s_min = this.estimateSmallestSingularValue(choleskyDecomposition.getG());
                    intersection = this.getBoundariesIntersections(p, s_min.vector, trustRegionRadius);
                    if (Math.abs(intersection.tmin) < Math.abs(intersection.tmax)) {
                        stepLength = intersection.tmin;
                    } else {
                        stepLength = intersection.tmax;
                    }
                    quadraticTerm = hessianPlusLambda.quadraticForm(p);

                    // Check stop criteria

                    relative_error = Math.pow(stepLength * s_min.value, 2) / (quadraticTerm + lambda.current * Math.pow(trustRegionRadius, 2));
                    if (relative_error <= this.k_hard) {
                        nuts.saxpy(stepLength, s_min.vector, p);
                        break;
                    }

                    // Update uncertainty bounds

                    lambda.ub = lambda.current;
                    lambda.lb = Math.max(lambda.lb, lambda.current - Math.pow(s_min.value, 2));

                    // Compute Cholesky factorization

                    hessianPlusLambda = this.hessian.clone();
                    hessianPlusLambda.addValueOnDiagonal(lambda_new);
                    choleskyDecomposition = new nuts.CholeskyDecomposition(hessianPlusLambda);

                    if (choleskyDecomposition.success) {
                        lambda.current = lambda_new;
                    } else {
                        lambda.lb = Math.max(lambda.lb, lambda_new);
                        lambda.current = Math.max(lambda.lb * lambda.ub, lambda.lb + this.UPDATE_COEFF * (lambda.ub - lambda.lb));
                    }

                } else { // Outside Boundary
                    //Check stop criteria
                    relative_error = Math.abs(pNorm - trustRegionRadius) / trustRegionRadius;
                    if (relative_error <= this.k_easy) {
                        break;
                    }

                    //Update uncertainty bounds
                    lambda.lb = lambda.current;
                    lambda.current = lambda_new;
                }
            } else if (choleskyDecomposition.success && nuts.norm(this.gradient) <= this.CLOSE_TO_ZERO) {
                  // Check for interior convergence
                if (lambda.current === 0) {
                    p = nuts.zeroVector(this.gradient.length);
                    hitsBoundary = false;
                    break;
                }
                s_min = this.estimateSmallestSingularValue(choleskyDecomposition.getG());
                stepLength = trustRegionRadius;

                  // Check stop criteria

                if (Math.pow(stepLength * s_min.value, 2) <= this.k_hard * lambda.current * Math.pow(trustRegionRadius, 2)) {
                    p = nuts.multiplyVectorByScalar(s_min.vector, stepLength);
                    break;
                }

                // update uncertainty bounds

                lambda.ub = lambda.current;
                lambda.lb = Math.max(lambda.lb, lambda.current - Math.pow(s_min.value, 2));

                // Update damping factor
                lambda.current = Math.max(Math.sqrt(lambda.lb * lambda.ub), lambda.lb + this.UPDATE_COEFF * (lambda.ub - lambda.lb));
            } else {// Unsuccesful factorization
                // Compute auxiliary term
                sls = this.singularLeadingSubmatrix(hessianPlusLambda, choleskyDecomposition.getG(), choleskyDecomposition.firstNonPositiveDefiniteLeadingSubmatrixSize);
                vNorm = nuts.norm(sls.vector);
                lambda.lb = Math.max(lambda.lb, lambda.current + sls.delta / Math.pow(vNorm, 2));
                lambda.current = Math.max(Math.sqrt(lambda.lb * lambda.ub), lambda.lb + this.UPDATE_COEFF * (lambda.ub - lambda.lb));
            }

        }

        return {
            step : p,
            hitsBoundary : hitsBoundary
        };
    },

    initialLambdas : function (trustRegionRadius) {
        'use strict';
        var lambda_initial,
            lambda_lb,
            lambda_ub,
            hessianFrobeniusNorm = 0,
            hessianInfiniteNorm = 0,
            minHessianDiagonal = this.hessian.get(0, 0),
            i,
            j,
            tempInfiniteNorm,
            gershgorin,
            EPSILON = 10e-4;

        for (i = 0; i < this.hessian.size; i += 1) {
            tempInfiniteNorm = 0;
            for (j = 0; j < this.hessian.size; j += 1) {
                hessianFrobeniusNorm += Math.pow(this.hessian.get(i, j), 2);
                tempInfiniteNorm += Math.abs(this.hessian.get(i, j));
            }
            hessianInfiniteNorm = Math.max(hessianInfiniteNorm, tempInfiniteNorm);
            minHessianDiagonal = Math.min(minHessianDiagonal, this.hessian.get(i, i));
        }
        hessianFrobeniusNorm = Math.sqrt(hessianFrobeniusNorm);
        gershgorin = this.gershgorin_bounds(this.hessian);

        //upper bound
        lambda_ub = Math.max(0, nuts.norm(this.gradient) / trustRegionRadius + Math.min(-gershgorin.lb, Math.min(hessianFrobeniusNorm, hessianInfiniteNorm)));

        //lower bound
        lambda_lb = Math.max(0, Math.max(-minHessianDiagonal, nuts.norm(this.gradient) / trustRegionRadius - Math.min(gershgorin.ub, Math.min(hessianFrobeniusNorm, hessianInfiniteNorm))));

        if (lambda_ub !== 0) {
            lambda_ub += EPSILON;
        }

        if (lambda_lb === 0) {
            lambda_initial = 0;
        } else {
            lambda_initial = Math.max(Math.sqrt(lambda_lb * lambda_ub), lambda_lb + this.UPDATE_COEFF * (lambda_ub - lambda_lb));
        }

        return {
            current : lambda_initial,
            lb : lambda_lb,
            ub : lambda_ub
        };

    },

    solveUpperTriangular : function (lowerTriangular, y) {
        'use strict';
        var x = y.slice(),
            n = lowerTriangular.size,
            i,
            k,
            sum;
        // LT x = y
        for (i = n - 1; i >= 0; i -= 1) {
            sum = x[i];
            for (k = i + 1; k < n; k += 1) {
                sum -= lowerTriangular.get(k, i) * x[k];
            }
            x[i] = sum / lowerTriangular.get(i, i);
        }
        return x;


    },

    solveLowerTriangular : function (lowerTriangular, b) {
        'use strict';
        var x = b.slice(),
            n = lowerTriangular.size,
            i,
            k,
            sum;
        // L x = b
        for (i = 0; i  < n; i += 1) {
            sum = b[i];
            for (k = i - 1; k >= 0; k -= 1) {
                sum -= lowerTriangular.get(i, k) * x[k];
            }
            x[i] = sum / lowerTriangular.get(i, i);
        }
        return x;

    },

    gershgorin_bounds : function (matrix) {
        'use strict';

        // Given a squared matrix, compute upper and lower bounds for its eigenvalues

        var lowerBound = 0,
            upperBound = 0,
            matrixDiagonal = [],
            matrixDiagonalAbsolute = [],
            matrixRowSums = [],
            lb = [],
            ub = [],
            i,
            j,
            rowSum;

        for (i = 0; i < matrix.numRows(); i += 1) {
            rowSum = 0;
            for (j = 0; j < matrix.numCols(); j += 1) {
                rowSum += Math.abs(matrix.get(i, j));
            }
            matrixRowSums.push(rowSum);
        }

        for (i = 0; i < matrix.numRows(); i += 1) {
            matrixDiagonal.push(matrix.get(i, i));
            matrixDiagonalAbsolute.push(Math.abs(matrix.get(i, i)));
        }

        for (i = 0; i < matrix.numRows(); i += 1) {
            lb.push(matrixDiagonal[i] + matrixDiagonalAbsolute[i] - matrixRowSums[i]);
            ub.push(matrixDiagonal[i] - matrixDiagonalAbsolute[i] + matrixRowSums[i]);
        }

        lowerBound = Math.min.apply(null, lb);
        upperBound = Math.max.apply(null, ub);
        return {
            lb : lowerBound,
            ub : upperBound
        };

    },

    estimateSmallestSingularValue : function (lowerTriangular) {
        'use strict';
        // see :  Golub, G. H., Van Loan, C. F. (2013), "Matrix computations". Forth Edition. JHU press. pp. 140-142.
        // https://github.com/scipy/scipy/blob/master/scipy/optimize/_trustregion_exact.py

        var n = lowerTriangular.size,
            p = nuts.zeroVector(n),
            y = nuts.zeroVector(n),
            y_plus,
            y_minus,
            i,
            k,
            p_plus = [],
            p_minus = [],
            v,
            vNorm,
            yNorm;

        for (k = 0; k < n; k += 1) {
            y_plus = (1 - p[k]) / lowerTriangular.get(k, k);
            y_minus = (-1 - p[k]) / lowerTriangular.get(k, k);
            for (i = k + 1; i < n; i += 1) {
                p_plus.push(p[i] + lowerTriangular.get(i, k) * y_plus);
                p_minus.push(p[i] + lowerTriangular.get(i, k) * y_minus);
            }

            if (Math.abs(y_plus) + this.norm1(p_plus) >= Math.abs(y_minus) + this.norm1(p_minus)) {
                y[k] = y_plus;
                for (i = k + 1; i < n; i += 1) {
                    p[i] = p_plus[i - k - 1];
                }
            } else {
                y[k] = y_minus;
                for (i = k + 1; i < n; i += 1) {
                    p[i] = p_minus[i - k - 1];
                }
            }

        }

        v = this.solveUpperTriangular(lowerTriangular, y);
        vNorm = nuts.norm(v);
        yNorm = nuts.norm(y);

        return {
            value : yNorm / vNorm,
            vector :  nuts.divideVectorByScalar(v, vNorm)
        };

    },

    norm1 : function (v) {
        'use strict';
        var result = 0,
            i;
        for (i = 0; i < v.length; i += 1) {
            result += Math.abs(v[i]);
        }
        return result;
    },

    getBoundariesIntersections : function (z, d, trustRegionRadius) {
        'use strict';
        //Solve the scalar quadratic equation ||z + t d|| == trust_radius.
        //This is like a line-sphere intersection.
        //Return the two values of t, sorted from low to high.

        var a = nuts.squaredNorm(d),
            b = 2 * nuts.dotProduct(z, d),
            c = nuts.squaredNorm(z) - trustRegionRadius * trustRegionRadius,
            sqrtDiscriminant = Math.sqrt(b * b - 4 * a * c),
            aux = b + sqrtDiscriminant * Math.sign(b),
            ta = -aux / (2 * a),
            tb = -2 * c / aux;
        return {
            tmin : Math.min(ta, tb),
            tmax : Math.max(ta, tb)
        };
    },

    singularLeadingSubmatrix : function (A, L, k) {
        'use strict';
        var delta = 0,
            j,
            l = new nuts.SquareMatrix(k),
            v = [],
            vtemp,
            i,
            u = nuts.zeroVector(k);
        for (j = 0; j < k - 1; j += 1) {
            delta += Math.pow(L.get(k - 1, j), 2);
        }
        delta -= A.get(k - 1, k - 1);

        for (i = 0; i < k - 1; i += 1) {
            for (j = 0; j <= i; j += 1) {
                l.set(i, j, L.get(i, j));
            }
            u[i] = L.get(k - 1, i);
        }

        v = nuts.zeroVector(A.size);
        v[k - 1] = 1;

        if (k !== 1) {
            vtemp = this.solveLowerTriangular(l, u);
            for (i = 0; i < k - 1; i += 1) {
                v[i] = vtemp[i];
            }
        }

        return {
            delta : delta,
            vector : v
        };

    }



};
