using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Gauss_Seidel_Sequential;
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
            if (comm.Rank == 0 && !matrix.isSquare)
            {
                Exception e = new Exception("Matrix must be square!");
                throw e;
            }

            benchmark bm = new benchmark(), bm2 = new benchmark();
            bm.start();

            int n = 0;
            int[] perm = new int[10]; int toggle = 0; Matrix lum = null;
            if (comm.Rank == 0)
            {
                n = matrix.dim1;
                lum = LUPDecompose(matrix, out perm, out toggle);
            }
            bm.pause();
            timeS += bm.getElapsedSeconds();

            bm.start();
            comm.Broadcast(ref n, 0);
            comm.Broadcast(ref lum, 0);
            if (comm.Rank != 0)
            {
                perm = new int[n];
            }
            comm.Broadcast(ref perm, 0);
            comm.Broadcast(ref toggle, 0);
            comm.Barrier();
            bm.pause();
            timeC += bm.getElapsedSeconds();

            if (lum == null)
            {
                return zeroLike(matrix);
            }

            bm.start();
            Double det = 0;
            if (comm.Rank == 0)
            {
                det = Determinant(lum, perm, toggle);
            }
            bm.pause();
            timeS += bm.getElapsedSeconds();
            bm.start();
            comm.Broadcast(ref det, 0);
            comm.Barrier();
            bm.pause();
            timeC += bm.getElapsedSeconds();
            if (det == 0) // not invertible
            {
                // still return for the sake of simplicity
                // Zero matrix * any matrix = zero matrix
                // so it's never a valid answer
                return zeroLike(matrix);
            }

            bm.pause();
            int slaves = comm.Size;
            Matrix jobDistro = Utils.splitJob(n, slaves);
            int startCol = 0, endCol = 0, size = (int)jobDistro[0, comm.Rank];
            for (int p = 0; p < slaves; p++)
            {
                if (p != comm.Rank)
                {
                    startCol += (int)jobDistro[0, p];
                }
                else
                {
                    endCol = startCol + (int)jobDistro[0, p] - 1;
                    break;
                }
            }
            bm.pause();
            timeP += bm.getElapsedSeconds();

            bm.start();
            Matrix result = new Matrix(n, size);
            for (int i = startCol; i < startCol + size; ++i)
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
                    result[j, i - startCol] = x[j];
                }
            }
            bm.pause();
            timeP += bm.getElapsedSeconds();

            bm.start();
            // collect result
            result = comm.Reduce(result, ConcatenateColumn, 0);
            bm.pause();
            timeP += bm.getElapsedSeconds();

            return result;
        }
    }
}
