using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Gauss_Seidel_Serial
{
    class Gauss_Seidel
    {
        // return true if it converges. Output: solution matrix, errors, loops it took
       public static Boolean solve(Matrix A, Matrix b, out Matrix x, out Matrix err, out int loops)
        {
            // check sanity
            if (!A.isSquare() || !b.isColumn() || (A.Height != b.Height))
            {
                Exception e = new Exception("Matrix A must be square! Matrix b must be a column matrix with the same height as matrix A!");
                throw e;
            }

            // follow samples in Wikipedia step by step https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method

            // decompose A into the sum of a lower triangular component L* and a strict upper triangular component U
            Matrix L, U;
            Matrix.Decompose(A, out L, out U);
            Matrix L_1 = ~L; // inverse of L*

            // x (at step k+1) = T * x (at step k) + C
            // where T = - (inverse of L*) * U and C = (inverse of L*) * b
            x = Matrix.zeroLike(b); // at step k
            Matrix new_x = Matrix.zeroLike(x); // at step k + 1
            Matrix T = -L_1 * U;
            Matrix C = L_1 * b;

            loops = 0;
            Boolean converge = false;
            int ITERATION_LIMIT = 1000; // if it still doesn't converge after this many loops, assume it won't converge and give up
            for (; loops < ITERATION_LIMIT; loops++)
            {
                new_x = T * x + C;
                if (converge = Matrix.AllClose(new_x, x, 1e-16)) // converge
                    break;
                x = new_x;
            }

            err = A * x - b;

            return converge;
        }

        // if u don't care about error
        public static Boolean solve(Matrix A, Matrix b, out Matrix x, out int loops)
        {
            Matrix err;
            return solve(A, b, out x, out err, out loops);
        }

        // just use this when u don't care about convergence and loop
        public static Matrix solve(Matrix A, Matrix b)
        {
            Matrix x;
            int loops;
            solve(A, b, out x, out loops);
            return x;
        }

        public static Boolean convergence(Matrix A)
        {
            //The convergence properties of the Gauss-Seidel method are dependent on the matrix A. Namely, the procedure is known to converge if either:
            //  A is symmetric positive-definite,[4] or
            //  A is strictly or irreducibly diagonally dominant.
            //The Gauss–Seidel method sometimes converges even if these conditions are not satisfied.
            return A.isDiagonallyDominant() || A.isPositiveDefinite();
        }
    }
}
