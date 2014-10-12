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

            benchmark bm = new benchmark();

            // decompose A into the sum of a lower triangular component L* and a strict upper triangular component U
            int size = A.Height;
            Matrix L, U;
            Matrix.Decompose(A, out L, out U);

            bm.start();
            // Inverse matrix L*
            Matrix L_1 = ~L;
            if (showBenchmark)
                Console.WriteLine("Matrix inversion took " + bm.getResult());

            // Main iteration: x (at step k+1) = T * x (at step k) + C
            // where T = - (inverse of L*) * U, and C = (inverse of L*) * b

            // init necessary variables
            x = Matrix.zeroLike(b); // at step k
            Matrix T = -L_1 * U;
            Matrix C = L_1 * b;

            // split T & C into groups of rows, each for one slave, according to the nature of this algorithm
            // each slave will have one piece of T & one piece of C stored locally. the rest of T & C is not needed
            // there might be cases where jobs > slaves, so some might get no job at all
            int slaves = comm.Size - 1;
            Matrix jobDistro = Utils.splitJob(size, slaves);
            List<Matrix> tPieces = new List<Matrix>(), cPieces = new List<Matrix>();
            List<int> startRows = new List<int>();
            int start = 0;
            for (int p = 0; p < slaves; p++)
            {
                Matrix tP = null, cP = null;
                startRows.Add(start);
                if (jobDistro[0, p] > 0) // only attempt to give job(s) if it has been assigned at least one
                {
                    int end = start + (int)jobDistro[0, p] - 1;
                    tP = Matrix.extractRows(T, start, end);
                    cP = Matrix.extractRows(C, start, end);
                    start = end + 1;
                }
                tPieces.Add(tP);
                cPieces.Add(cP);
            }

            // distribute pieces to slaves
            for (int p = 0; p < slaves; p++)
            {
                if (jobDistro[0, p] > 0)
                {
                    comm.Send("receivet", p + 1, 10); // it's not really necessary to tag data differently
                    comm.Send(tPieces[p], p + 1, 11); // as we notice slaves before sending data
                    comm.Send("receivec", p + 1, 10); // but whatever :v
                    comm.Send(cPieces[p], p + 1, 12);
                    comm.Send("setstartrow", p + 1, 10);
                    comm.Send(startRows[p], p + 1, 14);
                    comm.Send("setthreshold", p + 1, 10);
                    comm.Send(1e-15, p + 1, 15);
                }
            }

            // the actual iteration
            // if it still doesn't converge after this many loops, assume it won't converge and give up
            loops = 0;
            Boolean converge = false;
            int loopLimit = 100;
            bm.start();
            for (; loops < loopLimit; loops++)
            {
                // (re-)distributing x vector. Must be done every single loop
                // this loop needs x from the previous loop
                for (int p = 0; p < slaves; p++)
                {
                    if (jobDistro[0, p] > 0)
                    {
                        comm.Send("receivex", p + 1, 10);
                        comm.Send(x, p + 1, 13);
                    }
                }

                // "ask" slaves to execute the calculation step
                for (int p = 0; p < slaves; p++)
                {
                    if (jobDistro[0, p] > 0)
                        comm.Send("continue", p + 1, 10);
                }

                // collect result x
                int offset = 0;
                for (int p = 0; p < slaves; p++)
                {
                    if (jobDistro[0, p] > 0)
                    {
                        // collect piece of new_x from this slave
                        comm.Send("sendx", p + 1, 10);
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
                if (converge)
                {
                    loops++;
                    break;
                }
            }

            // command them to exit
            for (int p = 0; p < slaves; p++)
            {
                comm.Send("exit", p + 1, 10);
            }

            if (showBenchmark)
                Console.WriteLine("Iteration took " + bm.getResult());

            // round the result slightly
            x.Round(1e-14);
            err = A * x - b;
            err.Round(1e-14);

            return converge;
        }

        // method ranks 1+ call. does the actual iterating calculation
        public static void solveSub(ref Intracommunicator comm)
        {
            // receive x, C, and rows of T from main
            // loop new_x = T * x + C
            // only do something when master (rank 0) tells it to
            // command: continue, sendx, receivex, receivet, receivec, exit

            Matrix T = null, C = null, x = null, new_x = null;
            int startRow = 0;
            bool converge = false;
            double convergeThreshold = 1e-15;

            string command;
            do
            {
                command = comm.Receive<string>(0, 10);
                switch (command)
                {
                    case "continue": new_x = T * x + C; converge = Matrix.SomeClose(new_x, x, convergeThreshold, startRow); break;
                    case "sendx": comm.Send(new_x, 0, 10); break;
                    case "converge": comm.Send(converge, 0, 16); break;
                    case "receivet": T = comm.Receive<Matrix>(0, 11); break;
                    case "receivec": C = comm.Receive<Matrix>(0, 12); break;
                    case "receivex": x = comm.Receive<Matrix>(0, 13); break;
                    case "setstartrow": startRow = comm.Receive<int>(0, 14); break;
                    case "setthreshold": convergeThreshold = comm.Receive<Double>(0, 15); break;
                }
            } while (command != "exit");
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
