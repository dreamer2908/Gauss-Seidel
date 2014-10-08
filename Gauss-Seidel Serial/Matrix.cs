using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Gauss_Seidel_Serial
{
    public class Matrix
    {
        private readonly Double[,] _matrix;

        public Matrix(int dim1, int dim2)
        {
            _matrix = new Double[dim1, dim2];
        }

        // Height (dim1) = number of rows, width (dim2) = number of columns
        public int Height { get { return _matrix.GetLength(0); } }
        public int Width { get { return _matrix.GetLength(1); } }
        public int dim1 { get { return _matrix.GetLength(0); } }
        public int dim2 { get { return _matrix.GetLength(1); } }

        private static Random _r = new Random(); // random number generator had better be static

        public Double this[int x, int y]
        {
            get { return _matrix[x, y]; }
            set { _matrix[x, y] = value; }
        }

        #region unit, zero, random matrix
        static public Matrix unit(int size)
        {
            Matrix re = new Matrix(size, size);
            for (int x = 0; x < size; x++)
                re[x, x] = 1;
            return re;
        }

        static public Matrix zero(int size)
        {
            Matrix re = new Matrix(size, size);
            for (int x = 0; x < size; x++)
                for (int y = 0; y < size; y++)
                    re[x, y] = 0;
            return re;
        }

        static public Matrix zero(int dim1, int dim2)
        {
            Matrix re = new Matrix(dim1, dim2);
            for (int x = 0; x < dim1; x++)
                for (int y = 0; y < dim2; y++)
                    re[x, y] = 0;
            return re;
        }

        static public Matrix zeroLike(Matrix m)
        {
            int dim1 = m.dim1, dim2 = m.dim2;
            Matrix re = new Matrix(dim1, dim2);
            for (int x = 0; x < dim1; x++)
                for (int y = 0; y < dim2; y++)
                    re[x, y] = 0;
            return re;
        }

        public void zeroFill()
        {
            for (int x = 0; x < dim1; x++)
                for (int y = 0; y < dim2; y++)
                    this[x, y] = 0;
        }

        // fill with random double number in range [0, 1]
        public void randomFill()
        {
            for (int x = 0; x < dim1; x++)
                for (int y = 0; y < dim2; y++)
                    this[x, y] = _r.NextDouble();
        }

        static public Matrix random(int dim1, int dim2, Double _min, Double _max)
        {
            return random(dim1, dim2, _min, _max, false);
        }

        static public Matrix random(int dim1, int dim2, Double _min, Double _max, bool round)
        {
            Matrix re = new Matrix(dim1, dim2);
            re.randomFill();
            // check sanity
            Double min = (_min < _max) ? _min : _max;
            Double max = (_max > _min) ? _max : _min;
            if (min == max)
                max += 1;
            // scale to max - min range
            re *= Math.Abs(max - min);
            re = Matrix.OffsetValue(re, min);
            // round
            if (round)
                re.Round(1);
            return re;
        }
        #endregion

        #region Matrix properties
        public Boolean isSquare()
        {
            return (this.Height == this.Width);
        }

        public Boolean isColumn()
        {
            return (this.dim2 == 1);
        }

        public Boolean isRow()
        {
            return (this.dim1 == 1);
        }

        public Boolean isZero()
        {
            for (int i = 0; i < this.Height; i++)
                for (int j = 0; j < this.Width; j++)
                    if (this[i, j] != 0)
                        return false;
            return true;
        }

        public Boolean isSymmetric()
        {
            return (this.ToString() != Matrix.Transpose(this).ToString());
        }

        // calculate the determinant of this matrix. Supports any size
        public Double determinant()
        {
            if (!this.isSquare())
            {
                Exception e = new Exception("Matrix must be square!");
                throw e;
            }

            // see https://en.wikipedia.org/wiki/Determinant
            // and http://ctec.tvu.edu.vn/ttkhai/TCC/63_Dinh_thuc.htm
            // and http://mathworld.wolfram.com/Determinant.html
            // Using recursive implemention
            Double det = 0;
            int n = this.Height;

            if (n == 1)
            {
                return this[0, 0];
            }

            for (int j = 0; j < n; j++)
            {
                // copy everything but row i column j to create minorMatrix
                Matrix minorMatrix = new Matrix(n - 1, n - 1);
                for (int x = 0; x < n - 1; x++)
                {
                    int r = (x < 1) ? x : x + 1;
                    for (int y = 0; y < n - 1; y++)
                    {
                        int c = (y < j) ? y : y + 1;
                        minorMatrix[x, y] = this[r, c];
                    }
                }
                Double minorDet = minorMatrix.determinant();
                Double cofactor = (int)Math.Pow(-1, 1 + j) * minorDet;
                det = det + cofactor * this[1, j];
            }

            return det;
        }

        // check if this matrix can be inversed
        public Boolean inversible()
        {
            return (this.determinant() != 0);
        }
        #endregion

        #region operators

        // removed because m == null is a pain
        //public static Boolean operator ==(Matrix m1, Matrix m2)
        //{
        //    return (m1.ToString() == m2.ToString());
        //}

        //public static Boolean operator !=(Matrix m1, Matrix m2)
        //{
        //    return (m1.ToString() != m2.ToString());
        //}

        public static Matrix operator +(Matrix m1, Matrix m2)
        {
            return Add(m1, m2);
        }

        public static Matrix operator -(Matrix m1, Matrix m2)
        {
            return m1 + (m2 * -1);
        }

        public static Matrix operator -(Matrix m1)
        {
            return (m1 * -1);
        }

        public static Matrix operator *(Matrix m1, Matrix m2)
        {
            return Multiply(m1, m2);
        }

        public static Matrix operator *(Matrix m1, Double scalar)
        {
            return Multiply(m1, scalar);
        }

        public static Matrix operator *(Double scalar, Matrix m1)
        {
            return Multiply(m1, scalar);
        }

        public static Matrix operator /(Matrix m1, Matrix m2)
        {
            return Divide(m1, m2);
        }

        public static Matrix operator /(Matrix m1, Double scalar)
        {
            return Divide(m1, scalar);
        }

        public static Matrix operator /(Double scalar, Matrix m1)
        {
            return Divide(m1, scalar);
        }

        // Get the inverse matrix of m
        public static Matrix operator ~(Matrix m)
        {
            return Matrix.Inverse(m);
        }
        #endregion

        #region Basic matrix functions
        // Return the sum of m1 and m2. They must be in the same size
        public static Matrix Add(Matrix m1, Matrix m2)
        {
            if (m1.Width != m2.Width || m1.Height != m2.Height)
            {
                Exception e = new Exception("Two matrixes must be the same in size!");
                throw e;
            }

            Matrix re = new Matrix(m1.Height, m2.Width);
            for (int i = 0; i < re.Height; i++)
            {
                for (int j = 0; j < re.Width; j++)
                {
                    re[i, j] = m1[i, j] + m2[i, j];
                }
            }
            return re;
        }

        // Return the product of matrix m1 and m2.
        public static Matrix Multiply(Matrix m1, Matrix m2)
        {
            if (m1.Width != m2.Height) // wrong size
            {
                Exception e = new Exception("First matrix's width must be equal with second matrix's height!");
                throw e;
            }

            Matrix re = new Matrix(m1.Height, m2.Width);
            for (int i = 0; i < re.Height; i++)
            {
                for (int j = 0; j < re.Width; j++)
                {
                    re[i, j] = 0;
                    for (int k = 0; k < m1.Width; k++)
                    {
                        re[i, j] += m1[i, k] * m2[k, j];
                    }
                }
            }
            return re;
        }

        // Scalar multiply matrix m1 by scalar
        public static Matrix Multiply(Matrix m1, Double scalar)
        {
            Matrix re = new Matrix(m1.Height, m1.Width);
            for (int i = 0; i < re.Height; i++)
                for (int j = 0; j < re.Width; j++)
                    re[i, j] = m1[i, j] * scalar;
            return re;
        }

        // Scalar multiply matrix m1 by 1/scalar
        public static Matrix Divide(Matrix m1, Double scalar)
        {
            return Multiply(m1, 1/scalar);
        }

        // Calculate product of matrix m1 and the inverse matrix of m2.
        public static Matrix Divide(Matrix m1, Matrix m2)
        {
            if (m1.Width != m2.Height || !m2.isSquare()) // wrong size
            {
                Exception e = new Exception("First matrix's width must be equal with second matrix's height AND the second matrix must be square!");
                throw e;
            }

            return m1 * Matrix.Inverse(m2);
        }

        // Inverse matrix m. By definition, (inverse of m) * m = unit matrix
        public static Matrix Inverse(Matrix m)
        {
            int size = m.Height;
            if (!m.isSquare())
            {
                Exception e = new Exception("Matrix must be square!");
                throw e;
            }

            Double det = Determinant(m);// m.determinant();
            if (det == 0) // not inversible 
            {
                // still return for the sake of simplicity
                // Zero matrix * any matrix = zero matrix
                // so it's never a valid answer
                return Matrix.zero(size);
            }

            Matrix re = new Matrix(size, size);
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                {
                    // copy everything but row *j* column *i* to create minorMatrix 
                    // It's NOT the same as in determinant()
                    Matrix minorMatrix = new Matrix(size - 1, size - 1);
                    for (int x = 0; x < size - 1; x++)
                    {
                        int r = (x < j) ? x : x + 1;
                        for (int y = 0; y < size - 1; y++)
                        {
                            int c = (y < i) ? y : y + 1;
                            minorMatrix[x, y] = m[r, c];
                        }
                    }
                    Double minorDet = Determinant(minorMatrix);
                    re[i, j] = (int)Math.Pow(-1, i + j) * minorDet;
                }
            return re / det;
        }

        // add offset to every entries of m 
        public static Matrix OffsetValue(Matrix m, Double offset)
        {
            Matrix re = new Matrix(m.Height, m.Width);
            for (int i = 0; i < re.Height; i++)
                for (int j = 0; j < re.Width; j++)
                    re[i, j] = m[i, j] + offset;
            return re;
        }

        // reflect A over its main diagonal (which runs from top-left to bottom-right) to obtain A transpose
        public static Matrix Transpose(Matrix m)
        {
            Matrix re = new Matrix(m.Width, m.Height);
            for (int i = 0; i < re.Height; i++)
                for (int j = 0; j < re.Width; j++)
                    re[i, j] = m[j, i];
            return re;
        }

        // Return a separated copy of matrix m
        public static Matrix Duplicate(Matrix m)
        {
            Matrix re = new Matrix(m.Height, m.Width);
            for (int i = 0; i < re.Height; i++)
                for (int j = 0; j < re.Width; j++)
                    re[i, j] = m[i, j];
            return re;
        }

        // swap row #r1 and #r2
        public void swapRows(int r1, int r2)
        {
            for (int j = 0; j < this.dim2; j++)
            {
                Double tmp = this[r1, j];
                this[r1, j] = this[r2, j];
                this[r2, j] = tmp;
            }
        }

        #endregion

        #region Advanced matrix functions
        // Decompose matrix m into lower triangular matrix L and strict upper triangular matrix U
        public static void Decompose(Matrix m, out Matrix L, out Matrix U)
        {
            int size = m.Height;
            if (!m.isSquare())
            {
                Exception e = new Exception("Matrix must be square!");
                throw e;
            }
            if (size < 2)
            {
                Exception e = new Exception("Matrix must be at least 2x2!");
                throw e;
            }

            // lower triangular, so i >= j, or 0 <= i <= size, 0 <= j <= i
            L = new Matrix(size, size);

            for (int i = 0; i < size; i++)
                for (int j = 0; j <= i; j++)
                    L[i, j] = m[i, j];

            // strict upper triangular, so i < j
            U = new Matrix(size, size);
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    if (i < j)
                        U[i, j] = m[i, j];
        }

        // Check if the difference between each pair of values between two matrixes is less than maxDiff
        public static Boolean AllClose(Matrix m1, Matrix m2, Double maxDiff)
        {
            if (m1.Width != m2.Width || m1.Height != m2.Height)
            {
                Exception e = new Exception("Two matrixes must be the same in size!");
                throw e;
            }

            for (int i = 0; i < m1.Height; i++)
                for (int j = 0; j < m1.Width; j++)
                    if (Math.Abs(m1[i, j] - m2[i, j]) > Math.Abs(maxDiff))
                        return false;
            return true;
        }

        public static Boolean CholeskyDecompose(Matrix A, out Matrix L)
        {
            // see http://en.wikipedia.org/wiki/Cholesky_decomposition

            L = Matrix.zeroLike(A);

            if (!A.isSquare()) // only square matrix can be symmetric
            {
                // Console.WriteLine("only square matrix can be symmetric");
                return false;
            }
            if (A.ToString() != Matrix.Transpose(A).ToString()) // matrix A is symmetric <=> A = AT
            {
                // Console.WriteLine("matrix A is symmetric <=> A = AT");
                return false;
            }

            int size = A.Width;
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                {
                    if (i == j)
                    {
                        // calculate a sum or whatever
                        Double sum = 0;
                        for (int k = 0; k < j; k++)
                        {
                            sum += L[j, k] * L[j, k];
                        }
                        L[j, j] = Math.Sqrt(A[j, j] - sum);
                    }
                    else if (i > j)
                    {
                        // calculate a sum or whatever
                        Double sum = 0;
                        for (int k = 0; k < j; k++)
                        {
                            sum += L[i, k] * L[j, k];
                        }
                        L[i, j] = (1 / L[j, j]) * (A[i, j] - sum);
                    }
                }

            return true;
        }

        public Boolean isPositiveDefinite()
        {
            // In linear algebra, a symmetric n × n real matrix M is said to be positive definite if zTMz is positive for every non-zero column vector z of n real numbers. Here zT denotes the transpose of z.

            if (!this.isSquare()) // only square matrix can be symmetric
                return false;
            if (this.ToString() != Matrix.Transpose(this).ToString()) // matrix A is symmetric <=> A = AT
                return false;

            // now check if it's positive definite

            // generate vector z with random number from -1 to 1. make sure z is non-zero
            Matrix z = Matrix.random(this.dim1, 1, -1, 1);
            while (z.isZero())
                z = Matrix.random(this.dim1, 1, -1, 1);
            Matrix zT = Matrix.Transpose(z);

            // check zTMz
            Matrix zTMz = zT * this * z; // should be 1x1

            Boolean re = zTMz[0, 0] > 0;

            if (re)
            {
                // recheck with CholeskyDecompose
                Matrix L;
                Boolean worked = Matrix.CholeskyDecompose(this, out L);
                re = re && worked;

                Matrix mul = L * Matrix.Transpose(L);
                mul.Round(0.00000001);
                Matrix mul2 = Matrix.Duplicate(this);
                mul2.Round(0.00000001);
                worked = (mul.ToString() == mul2.ToString());
                re = re && worked;
            }

            return re;
        }

        // check if this matrix is diagonally dominant
        public Boolean isDiagonallyDominant()
        {
            // see https://en.wikipedia.org/wiki/Diagonally_dominant_matrix
            // In mathematics, a matrix is said to be diagonally dominant if for every row of the matrix,
            // the magnitude of the diagonal entry in a row is larger than or equal to the sum of the
            // magnitudes of all the other (non-diagonal) entries in that row.
            // More precisely, the matrix A is diagonally dominant if
            // | A[i,i] | >= sum(j!=i) of | A[i,j] |
            for (int i = 0; i < this.dim1; i++)
            {
                Double sumOfNonDiagonalEntriesInRow = 0;
                for (int j = 0; j < this.dim2; j++)
                    if (j != i)
                        sumOfNonDiagonalEntriesInRow += Math.Abs(this[i, j]);
                if (Math.Abs(this[i, i]) < sumOfNonDiagonalEntriesInRow)
                    return false;
            }
            return true;
        }

        // generate a symmetric positive-definite matrix
        public static Matrix generateSymmetricPositiveDefiniteMatrix(int size)
        {
            // Follow Daryl's answer here
            // https://math.stackexchange.com/questions/357980/matlab-code-for-generating-random-symmetric-positive-definite-matrix

            // generate a random n x n matrix
            Matrix A = Matrix.random(size, size, -10, 10);
            // construct a symmetric matrix from it
            A = A * A;
            // since A(i,j) < 1 by construction and a symmetric diagonally dominant matrix
            // is symmetric positive definite, which can be ensured by adding nI
            A = A + size * Matrix.unit(size);
            // well done
            return A;
        }

        // generate a diagonally dominant matrix
        public static Matrix generateDiagonallyDominantMatrix(int size, bool round, double min, double max)
        {
            // generate a random n x n matrix
            Matrix A = Matrix.random(size, size, min, max);
            if (round)
                A.Round(1);

            // modify diagonal line to make it dominant
            for (int i = 0; i < A.dim1; i++)
            {
                Double sumOfNonDiagonalEntriesInRow = 0;
                for (int j = 0; j < A.dim2; j++)
                    if (j != i)
                        sumOfNonDiagonalEntriesInRow += Math.Abs(A[i, j]);
                if (Math.Abs(A[i, i]) < sumOfNonDiagonalEntriesInRow)
                    A[i, i] = sumOfNonDiagonalEntriesInRow + _r.Next(1, 10);
            }

            return A;
        }

        // Doolittle LUP decomposition.
        public static Matrix LUPDecompose(Matrix matrix, out int[] perm, out int toggle)
        {
            // assumes matrix is square.
            int n = matrix.Width; // convenience
            Matrix result = Matrix.Duplicate(matrix);
            perm = new int[n];
            for (int i = 0; i < n; ++i) { perm[i] = i; }
            toggle = 1;
            for (int j = 0; j < n - 1; ++j) // each column
            {
                double colMax = Math.Abs(result[j, j]); // largest val in col j
                int pRow = j;
                for (int i = j + 1; i < n; ++i)
                {
                    if (result[i, j] > colMax)
                    {
                        colMax = result[i, j];
                        pRow = i;
                    }
                }
                if (pRow != j) // swap rows
                {
                    result.swapRows(pRow, j);
                    int tmp = perm[pRow]; // and swap perm info
                    perm[pRow] = perm[j];
                    perm[j] = tmp;
                    toggle = -toggle; // row-swap toggle
                }
                if (Math.Abs(result[j, j]) < 1.0E-20)
                    return null; // consider a throw
                for (int i = j + 1; i < n; ++i)
                {
                    result[i, j] /= result[j, j];
                    for (int k = j + 1; k < n; ++k)
                        result[i, k] -= result[i, j] * result[j, k];
                }
            } // main j column loop
            return result;
        }

        public static double Determinant(Matrix matrix)
        {
            int[] perm;
            int toggle;
            Matrix lum = LUPDecompose(matrix, out perm, out toggle);
            if (lum == null)
                return 0; // throw new Exception("Unable to compute MatrixDeterminant");
            double result = toggle;
            for (int i = 0; i < lum.Width; ++i)
                result *= lum[i, i];
            return result;
        }

        #endregion

        #region rouding
        // It's to fix stuff like 0.000000000001 or 0.999999999999999
        // >>>>>>>float
        public static Matrix Round(Matrix m, Double rouding)
        {
            Matrix re = new Matrix(m.Height, m.Width);
            for (int i = 0; i < re.Height; i++)
                for (int j = 0; j < re.Width; j++)
                    re[i, j] = RoundNum(m[i, j], rouding);
            return re;
        }

        // result will be rounded to multiple of rounding arg
        // for example: rounding = 0.05, result will be like 0.05, 0.1, 0.15, 0.2, 0.25, 0.3
        public void Round(Double rouding)
        {
            for (int i = 0; i < this.Height; i++)
                for (int j = 0; j < this.Width; j++)
                    this[i, j] = RoundNum(this[i, j], rouding);
        }

        private static Double RoundNum(Double num, Double rounding)
        {
            return Math.Floor(num / rounding + 0.5) * rounding;
        }
        #endregion

        #region overriding methods
        public override bool Equals(object obj)
        {
            return (this.GetHashCode() == obj.GetHashCode());
        }

        public override int GetHashCode()
        {
            return this.ToString().GetHashCode();
        }

        public override string ToString()
        {
            String re = "";
            for (int i = 0; i < this.Height; i++)
            {
                for (int j = 0; j < this.Width; j++)
                {
                    re += Format(this[i, j]) + " ";
                }
                if (i < this.Height - 1)
                    re += "\n";
            }
            return re;
        }

        public string ToString(double rouding)
        {
            String re = "";
            for (int i = 0; i < this.Height; i++)
            {
                for (int j = 0; j < this.Width; j++)
                {
                    re += Format(RoundNum(this[i, j], rouding)) + " ";
                }
                if (i < this.Height - 1)
                    re += "\n";
            }
            return re;
        }

        private string Format(Double n)
        {
            return String.Format("{0:0.###############}", n);
        }
        #endregion
    }
}
