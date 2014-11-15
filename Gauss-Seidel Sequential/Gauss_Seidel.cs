using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Gauss_Seidel_Sequential
{
    class Gauss_Seidel
    {
        public static bool showBenchmark = false;

        // return true if it converges. Output: solution matrix, errors, loops it took
        public static Boolean solve(Matrix A, Matrix b, out Matrix x, out Matrix err, out int loops)
        {
            // check sanity
            if (!A.isSquare || !b.isColumn || (A.Height != b.Height))
            {
                Exception e = new Exception("Matrix A must be square! Matrix b must be a column matrix with the same height as matrix A!");
                throw e;
            }

            // follow samples in Wikipedia step by step https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method

            benchmark bm = new benchmark(), bm2 = new benchmark(), bm3 = new benchmark();
            double sequential = 0, parallel = 0;
            bm.start();

            bm2.start();
            // decompose A into the sum of a lower triangular component L* and a strict upper triangular component U
            int size = A.Height;
            Matrix L, U;
            Matrix.Decompose(A, out L, out U);
            sequential += bm2.getElapsedSeconds();

            // Inverse matrix L*
            Matrix L_1 = Matrix.InverseAlt(L, ref sequential, ref parallel);

            // Main iteration: x (at step k+1) = T * x (at step k) + C
            // where T = - (inverse of L*) * U, and C = (inverse of L*) * b

            bm2.start();
            // init necessary variables
            x = Matrix.zeroLike(b); // at step k
            Matrix new_x; // at step k + 1
            sequential += bm2.getElapsedSeconds();
            bm2.start();
            Matrix T = -L_1 * U;
            Matrix C = L_1 * b;
            parallel += bm2.getElapsedSeconds();

            // the actual iteration
            // if it still doesn't converge after this many loops, assume it won't converge and give up
            bm2.start();
            loops = 0;
            Boolean converge = false;
            int loopLimit = 100;
            sequential += bm2.getElapsedSeconds();
            bm2.start();
            for (; loops < loopLimit; loops++)
            {
                new_x = T * x + C; // yup, only one line

                // consider it's converged if it changes less than threshold (1e-15)
                if (converge = Matrix.AllClose(new_x, x, 1e-15))
                {
                    x = new_x;
                    loops++;
                    break;
                }

                // save result
                x = new_x;
            }
            parallel += bm2.getElapsedSeconds();

            bm2.start();
            // round the result slightly
            x.Round(1e-14);
            err = A * x - b;
            err.Round(1e-14);
            sequential += bm2.getElapsedSeconds();

            bm.pause();
            if (showBenchmark)
            {
                Console.WriteLine("Sequential part took " + sequential + " secs.");
                Console.WriteLine("Parallel part took " + parallel + " secs.");
                Console.WriteLine("Total: " + bm.getResult() + " (" + bm.getElapsedSeconds() + " secs). Seq + Parallel: " + (sequential + parallel));
            }

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

        // Check if equations stored in matrix A should converge. Might still do sometimes even when this return false
        public static Boolean convergence(Matrix A)
        {
            //The convergence properties of the Gauss-Seidel method are dependent on the matrix A. Namely, the procedure is known to converge if either:
            //  A is symmetric positive-definite,[4] or
            //  A is strictly or irreducibly diagonally dominant.
            //The Gauss–Seidel method sometimes converges even if these conditions are not satisfied.
            return A.isDiagonallyDominant || A.isPositiveDefinite;
        }
    }
}
