using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Gauss_Seidel_Serial;
using MPI;

namespace Gauss_Seidel_Parallel
{
    [Serializable]
    class MatrixParallel : Matrix
    {
        public MatrixParallel(int dim1, int dim2)
            : base(dim1, dim2)
        {
        }

        public static Matrix Inverse(Matrix matrix, Intracommunicator comm)
        {
            double timeS = 0, timeP = 0, timeC = 0;
            return Inverse(matrix, comm, ref timeS, ref timeP, ref timeC);
        }

        public static Matrix Inverse(Matrix matrix, Intracommunicator comm, ref double timeS, ref double timeP, ref double timeC)
        {
            if (!matrix.isSquare)
            {
                gtfo(comm);
                Exception e = new Exception("Matrix must be square!");
                throw e;
            }

            benchmark bm = new benchmark(), bm2 = new benchmark();
            bm.start();

            int n = matrix.dim1;
            Matrix result = zeroLike(matrix);
            int[] perm;
            int toggle;
            Matrix lum = LUPDecompose(matrix, out perm, out toggle);
            if (lum == null)
            {
                gtfo(comm);
                return zeroLike(matrix);
            }

            Double det = Determinant(lum, perm, toggle);
            if (det == 0) // not invertible
            {
                // still return for the sake of simplicity
                // Zero matrix * any matrix = zero matrix
                // so it's never a valid answer
                gtfo(comm);
                return zeroLike(matrix);
            }

            int slaves = comm.Size - 1;
            Matrix jobDistro = Utils.splitJob(n, slaves);

            timeS += bm.getElapsedSeconds();
            bm.start();
            bm2.start();
            int start = 0;
            for (int p = 0; p < slaves; p++)
            {
                if (jobDistro[0, p] > 0) // only attempt to give job(s) if it has been assigned at least one
                {
                    comm.Send("recv_lum", p + 1, 10);
                    comm.Send(lum, p + 1, 11);
                    comm.Send("recv_perm", p + 1, 10);
                    comm.Send(perm, p + 1, 12);
                    comm.Send("set_offset", p + 1, 10);
                    comm.Send(start, p + 1, 13);
                    comm.Send("set_size", p + 1, 10);
                    comm.Send((int)jobDistro[0, p], p + 1, 14);
                    comm.Send("start", p + 1, 10);
                    start += (int)jobDistro[0, p];
                }
            }
            bm2.pause();
            timeC += bm2.getElapsedSeconds();
            for (int p = 0; p < slaves; p++)
            {
                if (jobDistro[0, p] > 0) // waiting for done message
                {
                    comm.Receive<string>(p + 1, 10);
                }
            }
            bm2.start();
            int offset = 0;
            for (int p = 0; p < slaves; p++)
            {
                if (jobDistro[0, p] > 0)
                {
                    // collect piece of result from this slave
                    comm.Send("send_re", p + 1, 10);
                    Matrix xp = comm.Receive<Matrix>(p + 1, 10);
                    // copy to its correct place in result
                    for (int c = 0; c < xp.Width; c++)
                        for (int r = 0; r < xp.Height; r++)
                            result[r, offset + c] = xp[r, c];
                    offset += xp.Width;
                }
            }

            gtfo(comm);
            bm.pause();
            timeP += bm.getElapsedSeconds();
            timeC += bm2.getElapsedSeconds();
            return result;
        }

        private static void gtfo(Intracommunicator comm)
        {
            for (int i = 1; i < comm.Size; i++)
            {
                comm.Send("exit", i, 10);
            }
        }

        public static void Inverse(Intracommunicator comm, ref double timeC)
        {
            Matrix lum = null, result = null;
            int[] perm = null;
            int offset = 0, size = 0, n = 0;
            benchmark bm = new benchmark();
            string command;
            do
            {
                command = comm.Receive<string>(0, 10);
                switch (command)
                {
                    case "start":
                        {
                            result = new Matrix(n, size);
                            for (int i = offset; i < offset + size; ++i)
                            {
                                double[] b = new double[n];
                                for (int j = 0; j < n; ++j)
                                {
                                    if (i == perm[j])
                                        b[j] = 1.0;
                                    else
                                        b[j] = 0.0;
                                }
                                double[] x = HelperSolve(lum, b);
                                for (int j = 0; j < n; ++j)
                                {
                                    result[j, i - offset] = x[j];
                                }
                            }
                            bm.pause();
                            comm.Send("done", 0, 10);
                            break;
                        }
                    case "send_re": bm.start(); comm.Send(result, 0, 10); timeC += bm.getElapsedSeconds(); break;
                    case "recv_lum": bm.start(); lum = comm.Receive<Matrix>(0, 11); timeC += bm.getElapsedSeconds(); n = lum.dim1; perm = new int[n]; break;
                    case "recv_perm": bm.start(); comm.Receive(0, 12, ref perm); timeC += bm.getElapsedSeconds(); break;
                    case "set_offset": offset = comm.Receive<int>(0, 13); break;
                    case "set_size": size = comm.Receive<int>(0, 14); break;
                }
            } while (command != "exit");
        }
    }
}
