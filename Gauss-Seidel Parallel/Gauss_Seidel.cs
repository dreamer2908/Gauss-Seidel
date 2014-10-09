using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Gauss_Seidel_Serial;
using MPI;

namespace Gauss_Seidel_Parallel
{
    class Gauss_Seidel
    {
        // return true if it converges. Output: solution matrix, errors, loops it took
        public static Boolean solve(Matrix A, Matrix b, out Matrix x, out Matrix err, out int loops, ref Intracommunicator comm)
        {
            // check sanity
            if (!A.isSquare || !b.isColumn || (A.Height != b.Height))
            {
                Exception e = new Exception("Matrix A must be square! Matrix b must be a column matrix with the same height as matrix A!");
                throw e;
            }

            // follow samples in Wikipedia step by step https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method

            // decompose A into the sum of a lower triangular component L* and a strict upper triangular component U
            int size = A.Height;
            Matrix L, U;
            Matrix.Decompose(A, out L, out U);
            //benchmark bm = new benchmark();
            //bm.start();
            Matrix L_1 = ~L; // inverse of L*
            //Console.WriteLine("Matrix inversion took " + bm.getResult());

            // x (at step k+1) = T * x (at step k) + C
            // where T = - (inverse of L*) * U and C = (inverse of L*) * b
            x = Matrix.zeroLike(b); // at step k
            Matrix new_x = Matrix.zeroLike(b); // at step k + 1
            Matrix T = -L_1 * U;
            Matrix C = L_1 * b;

            // split T into pieces (groups of rows)
            int slaves = comm.Size - 1;
            Matrix jobDistro = Utils.splitJob(size, slaves);
            List<Matrix> tPieces = new List<Matrix>(), cPieces = new List<Matrix>();
            int start = 0;
            for (int p = 0; p < slaves; p++)
            {
                Matrix tP = null, cP = null;
                if (jobDistro[0, p] > 0)
                {
                    int end = start + (int)jobDistro[0, p] - 1;
                    // Console.WriteLine("start = " + start.ToString() + ", end = " + end.ToString());
                    tP = Matrix.extractRows(T, start, end);
                    cP = Matrix.extractRows(C, start, end);
                    start = end + 1;
                }
                tPieces.Add(tP);
                cPieces.Add(cP);
            }

            // distribute pieces
            for (int p = 0; p < slaves; p++)
            {
                //Console.WriteLine("Distributing pieces to slave #" + (p + 1).ToString());
                if (jobDistro[0, p] > 0)
                {
                    comm.Send("receivet", p + 1, 10);
                    comm.Send(tPieces[p], p + 1, 11);
                    comm.Send("receivec", p + 1, 10);
                    comm.Send(cPieces[p], p + 1, 12);
                }
            }

            loops = 0;
            Boolean converge = false;
            int ITERATION_LIMIT = 100; // if it still doesn't converge after this many loops, assume it won't converge and give up
            //bm.start();
            for (; loops < ITERATION_LIMIT; loops++)
            {
                //Console.WriteLine("Running loop #" + loops.ToString());
                //Console.WriteLine(x.ToString());
                // distributing x
                for (int p = 0; p < slaves; p++)
                {
                    if (jobDistro[0, p] > 0)
                    {
                        comm.Send("receivex", p + 1, 10);
                        comm.Send(x, p + 1, 13);
                    }
                }

                // "ask" slaves to run this step
                for (int p = 0; p < slaves; p++)
                {
                    if (jobDistro[0, p] > 0)
                        comm.Send("continue", p + 1, 10);
                }

                // gather result
                int offset = 0;
                new_x = Matrix.zeroLike(x);
                for (int p = 0; p < slaves; p++)
                {
                    if (jobDistro[0, p] > 0)
                    {
                        comm.Send("sendx", p + 1, 10);
                        Matrix _xp = comm.Receive<Matrix>(p + 1, 10);
                        // copy _xp to x
                        for (int r = 0; r < _xp.Height; r++)
                            for (int c = 0; c < _xp.Width; c++)
                                new_x[offset + r, c] = _xp[r, c];
                        offset += _xp.Height;
                    }
                }

                //Console.WriteLine(new_x.ToString());
                // check converge
                if (converge = Matrix.AllClose(new_x, x, 1e-16))
                {
                    //Console.WriteLine("Converged.");
                    break;
                }

                x = new_x;
            }

            // ask them to exit
            for (int p = 0; p < slaves; p++)
            {
                comm.Send("exit", p + 1, 10);
            }

            //Console.WriteLine("Iteration took " + bm.getResult());
            err = A * x - b;

            return converge;
        }

        public static void solveSub(ref Intracommunicator comm)
        {
            // receive x, C, and rows of T from main
            // loop new_x = T * x + C
            // at the end of each loop, wait for main's commands
            // (continue, sendx, receivex, receivet, receivec, exit)

            Matrix T = null;
            Matrix C = null;
            Matrix x = null;
            Matrix new_x = null;

            string command;
            do
            {
                command = comm.Receive<string>(0, 10);
                // Console.WriteLine("Rank " + comm.Rank.ToString() + " received command \"" + command + "\"");
                if (command == "continue")
                {
                    new_x = T * x + C;
                }
                else if (command == "sendx")
                {
                    comm.Send(new_x, 0, 10);
                }
                else if (command == "receivex")
                {
                    x = comm.Receive<Matrix>(0, 13);
                }
                else if (command == "receivet")
                {
                    T = comm.Receive<Matrix>(0, 11);
                }
                else if (command == "receivec")
                {
                    C = comm.Receive<Matrix>(0, 12);
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
