using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Gauss_Seidel_Sequential;
using MPI;

namespace Gauss_Seidel_Parallel
{
    class Gauss_Seidel_Parallel
    {
        public static bool showBenchmark = false;

        // return true if it converges. Output: solution matrix, errors, loops it took
        public static Boolean solve(Matrix A, Matrix b, out Matrix x, out Matrix err, out int loops, Intracommunicator comm)
        {
            // check sanity. rank 0 only
            if (comm.Rank == 0 && (!A.isSquare || !b.isColumn || (A.Height != b.Height)))
            {
                Exception e = new Exception("Matrix A must be square! Matrix b must be a column matrix with the same height as matrix A!");
                throw e;
            }

            // follow samples in Wikipedia step by step https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method

            benchmark bm = new benchmark(), bm2 = new benchmark(), bm3 = new benchmark();
            double sequential = 0, parallel = 0, communication = 0;
            bm.start();

            bm2.start();
            // decompose A into the sum of a lower triangular component L* and a strict upper triangular component U
            int size = 0; Matrix L = null, U = null, L_1;
            if (comm.Rank == 0)
            {
                size = A.Height;
                Matrix.Decompose(A, out L, out U);
            }
            bm2.pause();
            sequential += bm2.getElapsedSeconds();

            bm2.start();
            comm.Broadcast(ref size, 0);
            comm.Broadcast(ref U, 0);
            comm.Broadcast(ref b, 0);
            bm2.pause();
            communication += bm2.getElapsedSeconds();

            // Inverse matrix L*
            comm.Barrier();
            L_1 = MatrixParallel.Inverse(L, comm, ref sequential, ref parallel, ref communication);

            // Main iteration: x (at step k+1) = T * x (at step k) + C
            // where T = - (inverse of L*) * U, and C = (inverse of L*) * b

            // split T & C into groups of rows, each for one slave, according to the nature of this algorithm
            // each slave will have one piece of T & one piece of C stored locally. the rest of T & C is not needed
            // there might be cases where jobs > slaves, so some might get no job at all
            // Changes: only split L_1. Slaves will calculate T & C (pieces) themselves
            bm2.start();
            Matrix jobDistro = Utils.splitJob(size, comm.Size);
            int startRow = 0, endRow = 0, myJobSize = (int)jobDistro[0, comm.Rank];
            for (int p = 0; p < comm.Size; p++)
            {
                if (p != comm.Rank)
                {
                    startRow += (int)jobDistro[0, p];
                }
                else
                {
                    endRow = startRow + (int)jobDistro[0, p] - 1;
                    break;
                }
            }
            Matrix[] L_1Ps = new Matrix[comm.Size];
            if (comm.Rank == 0)
            {
                int slaveStart = 0;
                for (int p = 0; p < comm.Size; p++)
                {
                    L_1Ps[p] = Matrix.extractRows(L_1, slaveStart, slaveStart + (int)jobDistro[0, p] - 1);
                    slaveStart += (int)jobDistro[0, p];
                }
            }
            bm2.pause();
            sequential += bm2.getElapsedSeconds();

            bm2.start();
            Matrix L_1P = comm.Scatter(L_1Ps, 0);
            bm2.pause();
            communication += bm2.getElapsedSeconds();
            bm2.start();
            Matrix T = -L_1P * U; Matrix C = L_1P * b;
            bm2.pause();
            parallel += bm2.getElapsedSeconds();

            // the actual iteration
            // if it still doesn't converge after this many loops, assume it won't converge and give up
            Boolean converge = false;
            int loopLimit = 100;
            x = Matrix.zeroLike(b); // at step k
            for (loops = 0; loops < loopLimit; loops++)
            {
                bm3.start();
                // (re-)distributing x vector. Must be done every single loop
                // this loop needs x from the previous loop
                comm.Broadcast(ref x, 0);
                bm3.pause();
                communication += bm3.getElapsedSeconds();

                // calculation step
                bm3.start();
                comm.Barrier();
                Matrix new_x = T * x + C;

                // check convergence
                converge = Matrix.SomeClose(new_x, x, 1e-15, startRow);

                // collect result x
                comm.Barrier();
                x = comm.Reduce(new_x, Matrix.Concatenate, 0);

                // collect convergence. consider converged if ALL slaves claim so
                converge = comm.Reduce(converge, bothTrue, 0);
                comm.Broadcast(ref converge, 0); // make sure EVERYONE breaks/coninues
                bm3.pause();
                parallel += bm3.getElapsedSeconds();
                if (converge)
                {
                    loops++;
                    break;
                }
            }

            bm2.start();
            // round the result slightly
            err = null;
            if (comm.Rank == 0)
            {
                x.Round(1e-14);
                err = A * x - b;
                err.Round(1e-14);
            }
            bm2.pause();
            sequential += bm2.getElapsedSeconds();

            bm.pause();
            if (showBenchmark)
            {
                Console.WriteLine("Sequential part took " + sequential + " secs.");
                Console.WriteLine("Parallel part took " + parallel + " secs.");
                Console.WriteLine("Communication took " + communication + " secs.");
                Console.WriteLine("Total: " + bm.getResult() + " (" + bm.getElapsedSeconds() + " secs). Seq + Parallel: " + (sequential + parallel));
            }

            return converge;
        }

        private static Boolean bothTrue(Boolean v1, Boolean v2)
        {
            return v1 && v2;
        }

        // if u don't care about error
        public static Boolean solve(Matrix A, Matrix b, out Matrix x, out int loops, Intracommunicator comm)
        {
            Matrix err;
            return solve(A, b, out x, out err, out loops, comm);
        }

        // just use this when u don't care about convergence and loop
        public static Matrix solve(Matrix A, Matrix b, Intracommunicator comm)
        {
            Matrix x;
            int loops;
            solve(A, b, out x, out loops, comm);
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
