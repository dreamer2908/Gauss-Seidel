using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Gauss_Seidel_Serial;
using MPI;

namespace Gauss_Seidel_Parallel
{
    class Gauss_Seidel_Parallel
    {
        public static bool showBenchmark = false;

        // return true if it converges. Output: solution matrix, errors, loops it took
        // rank 0 (master) calls this, while ranks 1+ call solveSub
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

            comm.Broadcast(ref size, 0);
            comm.Broadcast(ref U, 0);
            comm.Broadcast(ref b, 0);
            sequential += bm2.getElapsedSeconds();

            // Inverse matrix L*
            comm.Barrier();
            L_1 = MatrixParallel.Inverse(L, comm, ref sequential, ref parallel, ref communication);
            comm.Broadcast(ref L_1, 0);

            // Main iteration: x (at step k+1) = T * x (at step k) + C
            // where T = - (inverse of L*) * U, and C = (inverse of L*) * b

            bm2.start();
            // init necessary variables
            x = Matrix.zeroLike(b); // at step k

            // split T & C into groups of rows, each for one slave, according to the nature of this algorithm
            // each slave will have one piece of T & one piece of C stored locally. the rest of T & C is not needed
            // there might be cases where jobs > slaves, so some might get no job at all
            // Changes: only split L_1. Slaves will calculate T & C (pieces) themselves
            // Changes: slaves will split L_1 themselves
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
            Matrix L_1P = Matrix.extractRows(L_1, startRow, endRow);
            Matrix T = -L_1P * U; Matrix C = L_1P * b;
            bm2.pause();
            parallel += bm2.getElapsedSeconds();

            // the actual iteration
            // if it still doesn't converge after this many loops, assume it won't converge and give up
            Boolean converge = false;
            int loopLimit = 100;
            for (loops = 0; loops < loopLimit; loops++)
            {
                bm3.start();
                // (re-)distributing x vector. Must be done every single loop
                // this loop needs x from the previous loop
                comm.Broadcast(ref x, 0);
                comm.Barrier();
                bm3.pause();
                communication += bm3.getElapsedSeconds();

                // calculation step
                bm3.start();
                Matrix new_x = T * x + C;

                // check convergence
                converge = Matrix.SomeClose(new_x, x, 1e-15, startRow);

                // collect result x
                comm.Barrier();
                x = comm.Reduce(new_x, Matrix.Concatenate, 0);

                // collect convergence. consider converged if ALL slaves claim so
                converge = comm.Reduce(converge, bothTrue, 0);
                comm.Broadcast(ref converge, 0);
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
            sequential += bm2.getElapsedSeconds();

            bm.pause();
            if (showBenchmark)
            {
                Console.WriteLine("Sequential part took " + sequential + " secs.");
                Console.WriteLine("Parallel part took " + parallel + " secs.");
                Console.WriteLine("Communication took " + communication + " secs.");
                Console.WriteLine("Total: " + bm.getResult() + " (" + bm.getElapsedSeconds() + " secs). Sum: " + (sequential + parallel));
            }

            return converge;
        }

        private static Boolean bothTrue(Boolean v1, Boolean v2)
        {
            return v1 && v2;
        }

        private static void sendMsgToAllSlaves(Intracommunicator comm, string msg)
        {
            for (int i = 1; i < comm.Size; i++)
            {
                comm.Send(msg, i, 10);
            }
        }

        // method ranks 1+ call. does the actual iterating calculation
        public static void solveSub(Intracommunicator comm)
        {
            // receive x, C, and rows of T from main
            // loop new_x = T * x + C
            // only do something when master (rank 0) tells it to
            // command: continue, sendx, receivex, receivet, receivec, exit

            Matrix T = null, C = null, x = null, new_x = null, L_1 = null, U = null, b = null;
            int startRow = 0;
            bool converge = false;
            double convergeThreshold = 1e-15, timeC = 0;
            benchmark bm = new benchmark();

            string command;
            do
            {
                command = comm.Receive<string>(0, 10);
                switch (command)
                {
                    case "calc_x": new_x = T * x + C; converge = Matrix.SomeClose(new_x, x, convergeThreshold, startRow); comm.Send("done", 0, 10); break;
                    case "calc_tc": T = -L_1 * U; C = L_1 * b; break;
                    case "send_x": comm.Send(new_x, 0, 10); break;
                    case "converge": comm.Send(converge, 0, 16); break;
                    case "rec_t": bm.start(); T = comm.Receive<Matrix>(0, 11); timeC += bm.getElapsedSeconds(); break;
                    case "rec_c": bm.start(); C = comm.Receive<Matrix>(0, 12); timeC += bm.getElapsedSeconds(); break;
                    case "rec_x": bm.start(); x = comm.Receive<Matrix>(0, 13); timeC += bm.getElapsedSeconds(); break;
                    case "set_startrow": startRow = comm.Receive<int>(0, 14); break;
                    case "set_threshold": convergeThreshold = comm.Receive<Double>(0, 15); break;
                    case "rec_l": bm.start(); L_1 = comm.Receive<Matrix>(0, 16); timeC += bm.getElapsedSeconds(); break;
                    case "rec_u": bm.start(); U = comm.Receive<Matrix>(0, 17); timeC += bm.getElapsedSeconds(); break;
                    case "rec_b": bm.start(); b = comm.Receive<Matrix>(0, 18); timeC += bm.getElapsedSeconds(); break;
                    case "inverse": MatrixParallel.Inverse(comm, ref timeC); break;
                }
            } while (command != "exit");
            if (showBenchmark)
                Console.WriteLine("Communication time @ rank " + comm.Rank + ": " + timeC.ToString());
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
