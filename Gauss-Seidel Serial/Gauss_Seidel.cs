using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Gauss_Seidel_Serial
{
    class Gauss_Seidel
    {
        public static Matrix solve(Matrix A, Matrix b)
        {
            // check sanity
            if (!A.isSquare() || !b.isColumn() || (A.Height != b.Height))
            {
                Exception e = new Exception("Matrix A must be square! Matrix b must be a column matrix with the same height as matrix A!");
                throw e;
            }

            int ITERATION_LIMIT = 1000;

            Matrix x = Matrix.zeroLike(b);
            Matrix L, U, L_1, T, C, new_x;
            Matrix.Decompose(A, out L, out U);
            L_1 = ~L;
            T = -L_1 * U;
            C = L_1 * b;

            int loops = 0;
            for (; loops < ITERATION_LIMIT; loops++)
            {
                new_x = T * x + C;
                if (Matrix.AllClose(new_x, x, 1e-8))
                    break;
                x = new_x;
            }

            return x;
        }
    }
}
