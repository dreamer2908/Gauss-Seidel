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
        public static Boolean solve(Matrix A, Matrix b, out Matrix x, out Matrix err, out int loops, ref Intracommunicator comm)
        {
            // check sanity
            if (!A.isSquare || !b.isColumn || (A.Height != b.Height))
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
            int size = A.Height;
            Matrix L, U, L_1;
            Matrix.Decompose(A, out L, out U);
            sequential += bm2.getElapsedSeconds();

            bm2.start();
            // Inverse matrix L*
            sendMsgToAllSlaves(comm, "inverse");
            bm2.pause();
            parallel += bm2.getElapsedSeconds();
            communication += bm2.getElapsedSeconds();
            L_1 = MatrixParallel.Inverse(L, comm, ref sequential, ref parallel, ref communication);

            // Main iteration: x (at step k+1) = T * x (at step k) + C
            // where T = - (inverse of L*) * U, and C = (inverse of L*) * b

            bm2.start();
            // init necessary variables
            x = Matrix.zeroLike(b); // at step k

            // split T & C into groups of rows, each for one slave, according to the nature of this algorithm
            // each slave will have one piece of T & one piece of C stored locally. the rest of T & C is not needed
            // there might be cases where jobs > slaves, so some might get no job at all
            // Changes: only split L_1. Slaves will calculate T & C (pieces) themselves
            int slaves = comm.Size - 1;
            Matrix jobDistro = Utils.splitJob(size, slaves);
            List<Matrix> L_1Pieces = new List<Matrix>();
            List<int> startRows = new List<int>();
            int start = 0;
            for (int p = 0; p < slaves; p++)
            {
                Matrix L_1P = null;
                startRows.Add(start);
                if (jobDistro[0, p] > 0) // only attempt to give job(s) if it has been assigned at least one
                {
                    int end = start + (int)jobDistro[0, p] - 1;
                    L_1P = Matrix.extractRows(L_1, start, end);
                    start = end + 1;
                }
                L_1Pieces.Add(L_1P);
            }
            bm2.pause();
            sequential += bm2.getElapsedSeconds();
            communication += bm2.getElapsedSeconds();

            bm2.start();
            // distribute pieces to slaves
            for (int p = 0; p < slaves; p++)
            {
                if (jobDistro[0, p] > 0)
                {
                    // it's not really necessary to tag data differently as we notice slaves before sending data
                    // but whatever :v
                    comm.Send("rec_l", p + 1, 10);
                    comm.Send(L_1Pieces[p], p + 1, 16);
                    comm.Send("rec_u", p + 1, 10);
                    comm.Send(U, p + 1, 17);
                    comm.Send("rec_b", p + 1, 10);
                    comm.Send(b, p + 1, 18);
                    comm.Send("set_startrow", p + 1, 10);
                    comm.Send(startRows[p], p + 1, 14);
                    comm.Send("set_threshold", p + 1, 10);
                    comm.Send(1e-15, p + 1, 15);
                    comm.Send("calc_tc", p + 1, 10);
                }
            }
            bm2.pause();
            parallel += bm2.getElapsedSeconds();
            communication += bm2.getElapsedSeconds();

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
                bm3.start();
                // (re-)distributing x vector. Must be done every single loop
                // this loop needs x from the previous loop
                for (int p = 0; p < slaves; p++)
                {
                    if (jobDistro[0, p] > 0)
                    {
                        comm.Send("rec_x", p + 1, 10);
                        comm.Send(x, p + 1, 13);
                    }
                }

                // "ask" slaves to execute the calculation step
                for (int p = 0; p < slaves; p++)
                {
                    if (jobDistro[0, p] > 0)
                        comm.Send("calc_x", p + 1, 10);
                }
                communication += bm3.getElapsedSeconds();

                // waiting for done message
                for (int p = 0; p < slaves; p++)
                {
                    if (jobDistro[0, p] > 0)
                    {
                        comm.Receive<string>(p + 1, 10);
                    }
                }
                bm3.start();
                // collect result x
                int offset = 0;
                for (int p = 0; p < slaves; p++)
                {
                    if (jobDistro[0, p] > 0)
                    {
                        // collect piece of new_x from this slave
                        comm.Send("send_x", p + 1, 10);
                        Matrix xp = comm.Receive<Matrix>(p + 1, 10);
                        // copy to its correct place in x
                        for (int r = 0; r < xp.Height; r++)
                            for (int c = 0; c < xp.Width; c++)
                                x[offset + r, c] = xp[r, c];
                        offset += xp.Height;
                    }
                }

                // collect convergence. consider converged if ALL slaves claim so
                for (int p = 0; p < slaves; p++)
                {
                    if (jobDistro[0, p] > 0)
                    {
                        comm.Send("converge", p + 1, 10);
                        converge = comm.Receive<bool>(p + 1, 16);
                        if (!converge)
                            break;
                    }
                }
                bm3.pause();
                communication += bm3.getElapsedSeconds();
                if (converge)
                {
                    loops++;
                    break;
                }
            }

            bm3.start();
            // command them to exit
            sendMsgToAllSlaves(comm, "exit");
            communication += bm3.getElapsedSeconds();
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
                Console.WriteLine("Communication took " + communication + " secs.");
                Console.WriteLine("Total: " + bm.getResult() + " (" + bm.getElapsedSeconds() + " secs). Sum: " + (sequential + parallel));
            }

            return converge;
        }

        private static void sendMsgToAllSlaves(Intracommunicator comm, string msg)
        {
            for (int i = 1; i < comm.Size; i++)
            {
                comm.Send(msg, i, 10);
            }
        }

        // method ranks 1+ call. does the actual iterating calculation
        public static void solveSub(ref Intracommunicator comm)
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
        public static Boolean solve(Matrix A, Matrix b, out Matrix x, out int loops, ref Intracommunicator comm)
        {
            Matrix err;
            return solve(A, b, out x, out err, out loops, ref comm);
        }

        // just use this when u don't care about convergence and loop
        public static Matrix solve(Matrix A, Matrix b, ref Intracommunicator comm)
        {
            Matrix x;
            int loops;
            solve(A, b, out x, out loops, ref comm);
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
