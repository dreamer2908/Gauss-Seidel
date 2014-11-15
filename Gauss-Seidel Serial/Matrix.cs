using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Gauss_Seidel_Serial
{
    [Serializable]
    public class Matrix
    {
        private readonly Double[,] _matrix;

        public Matrix(int dim1, int dim2)
        {
            _matrix = new Double[dim1, dim2];
        }

        public Double this[int x, int y]
        {
            get { return _matrix[x, y]; }
            set { _matrix[x, y] = value; }
        }

        // Height (dim1) = number of rows, width (dim2) = number of columns
        public int Height { get { return _matrix.GetLength(0); } }
        public int Width { get { return _matrix.GetLength(1); } }
        public int dim1 { get { return _matrix.GetLength(0); } }
        public int dim2 { get { return _matrix.GetLength(1); } }
        public int row { get { return _matrix.GetLength(0); } }
        public int column { get { return _matrix.GetLength(1); } }

        private static Random _r = new Random(); // random number generator had better be static

        #region Unit, zero, random matrix
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

        public void numFill(int n)
        {
            for (int x = 0; x < dim1; x++)
                for (int y = 0; y < dim2; y++)
                    this[x, y] = n;
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
        public Boolean isSquare
        {
            get { return (this.Height == this.Width); }
        }

        public Boolean isColumn
        {
            get { return (this.dim2 == 1); }
        }

        public Boolean isRow
        {
            get { return (this.dim1 == 1); }
        }

        public Boolean isZero
        {
            get
            {
                for (int i = 0; i < this.Height; i++)
                    for (int j = 0; j < this.Width; j++)
                        if (this[i, j] != 0)
                            return false;
                return true;
            }
        }

        public Boolean isSymmetric
        {
            get { return (this.ToString() != Matrix.Transpose(this).ToString()); }
        }

        // check if this matrix can be inversed
        public Boolean invertible
        {
            get { return (Matrix.Determinant(this) != 0); }
        }

        // get average value
        public double avgValue
        {
            get
            {
                double val = 0;
                for (int i = 0; i < this.Height; i++)
                    for (int j = 0; j < this.Width; j++)
                        val += this[i, j];
                return val / (this.Height * this.Width);
            }
        }

        public double totalValue
        {
            get
            {
                double val = 0;
                for (int i = 0; i < this.Height; i++)
                    for (int j = 0; j < this.Width; j++)
                        val += this[i, j];
                return val;
            }
        }

        public double minValue
        {
            get
            {
                double min = this[0,0];
                for (int i = 0; i < this.Height; i++)
                    for (int j = 0; j < this.Width; j++)
                        if (this[i, j] < min)
                            min = this[i, j];
                return min;
            }
        }

        public double maxValue
        {
            get
            {
                double max = this[0, 0];
                for (int i = 0; i < this.Height; i++)
                    for (int j = 0; j < this.Width; j++)
                        if (this[i, j] > max)
                            max = this[i, j];
                return max;
            }
        }
        #endregion

        #region Operators

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
            return UnsafeMultiplication(m1, m2);
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
            return Matrix.InverseAlt(m);
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

        // see http://www.bratched.com/en/home/dotnet/48-fun-with-matrix-multiplication-and-unsafe-code.html
        public unsafe static Matrix UnsafeMultiplication(Matrix m1, Matrix m2)
        {
            int h = m1.Height;
            int w = m2.Width;
            int l = m1.Width;
            Matrix resultMatrix = new Matrix(h, w);
            unsafe
            {
                fixed (double* pm = resultMatrix._matrix, pm1 = m1._matrix, pm2 = m2._matrix)
                {
                    int i1, i2;
                    for (int i = 0; i < h; i++)
                    {
                        i1 = i * l;
                        for (int j = 0; j < w; j++)
                        {
                            i2 = j;
                            double res = 0;
                            for (int k = 0; k < l; k++, i2 += w)
                            {
                                res += pm1[i1 + k] * pm2[i2];
                            }
                            pm[i * w + j] = res;
                        }
                    }
                }
            }
            return resultMatrix;
        }

        // Scalar multiply matrix m1 by 1/scalar
        public static Matrix Divide(Matrix m1, Double scalar)
        {
            return Multiply(m1, 1/scalar);
        }

        // Calculate product of matrix m1 and the inverse matrix of m2.
        public static Matrix Divide(Matrix m1, Matrix m2)
        {
            if (m1.Width != m2.Height || !m2.isSquare) // wrong size
            {
                Exception e = new Exception("First matrix's width must be equal with second matrix's height AND the second matrix must be square!");
                throw e;
            }

            return m1 * Matrix.Inverse(m2);
        }

        // calculate the determinant of this matrix. Supports any size
        // but slow with large matrix (like 9+)
        public Double determinant()
        {
            if (!this.isSquare)
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

        // Inverse matrix m. By definition, (inverse of m) * m = unit matrix
        // slow with large matrix
        public static Matrix Inverse(Matrix m)
        {
            int size = m.Height;
            if (!m.isSquare)
            {
                Exception e = new Exception("Matrix must be square!");
                throw e;
            }

            Double det = m.determinant();
            if (det == 0) // not invertible
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
                    Double minorDet = minorMatrix.determinant();
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
            // Changed from traditional element-by-element copying to memory copying
            Buffer.BlockCopy(m._matrix, 0, re._matrix, 0, m.Width * m.Height * sizeof(Double));
            //for (int i = 0; i < re.Height; i++)
            //    for (int j = 0; j < re.Width; j++)
            //        re[i, j] = m[i, j];
            return re;
        }

        // swap row #r1 and #r2
        public void swapRows(int r1, int r2)
        {
            int length = this.dim2;
            Double[] tmp = new Double[length];
            // Changed from traditional element-by-element copying to memory copying
            Buffer.BlockCopy(_matrix, (r1 * length) * sizeof(Double), tmp, 0, length * sizeof(Double));
            Buffer.BlockCopy(_matrix, (r2 * length) * sizeof(Double), _matrix, (r1 * length) * sizeof(Double), length * sizeof(Double));
            Buffer.BlockCopy(tmp, 0, _matrix, (r2 * length) * sizeof(Double), length * sizeof(Double));
            //for (int j = 0; j < this.dim2; j++)
            //{
            //    Double tmp = this[r1, j];
            //    this[r1, j] = this[r2, j];
            //    this[r2, j] = tmp;
            //}
        }

        // absolute value
        public void Abs()
        {
            for (int i = 0; i < this.Height; i++)
                for (int j = 0; j < this.Width; j++)
                    this[i, j] = Math.Abs(this[i, j]);
        }

        public static Matrix Abs(Matrix m)
        {
            Matrix re = zeroLike(m);
            for (int i = 0; i < m.Height; i++)
                for (int j = 0; j < m.Width; j++)
                    re[i, j] = Math.Abs(m[i, j]);
            return re;
        }

        public static Matrix extractRows(Matrix m, int _start, int _end)
        {
            int start = _start, end = _end;
            if (start > end)
            {
                int tmp = start;
                start = end;
                end = tmp;
            }
            if (start < 0 || end > m.Height)
            {
                Exception e = new Exception("Start can't be negative, and end can't exceed the height of matrix m!");
                throw e;
            }
            int rows = end - start + 1;
            if (rows < 1)
            {
                Exception e = new Exception("Number of rows to extract must be at least 1!");
                throw e;
            }
            Matrix re = new Matrix(rows, m.Width);
            // Changed from traditional element-by-element copying to memory copying
            Buffer.BlockCopy(m._matrix, (start * m.Width) * sizeof(Double), re._matrix, 0, rows * m.Width * sizeof(Double));
            //for (int i = 0; i < rows; i++)
            //    for (int j = 0; j < m.Width; j++)
            //        re[i, j] = m[i + start, j];
            return re;
        }

        public static Matrix extractColumns(Matrix m, int _start, int _end)
        {
            int start = _start, end = _end;
            if (start > end)
            {
                int tmp = start;
                start = end;
                end = tmp;
            }
            if (start < 0 || end > m.column)
            {
                Exception e = new Exception("Start can't be negative, and end can't exceed the width of matrix m!");
                throw e;
            }
            int columns = end - start + 1;
            if (columns < 1)
            {
                Exception e = new Exception("Number of columns to extract must be at least 1!");
                throw e;
            }
            Matrix re = new Matrix(m.row, columns);
            for (int i = 0; i < m.row; i++)
                for (int j = 0; j < columns; j++)
                    re[i, j] = m[i, j + start];
            return re;
        }

        public static Matrix Concatenate(Matrix m1, Matrix m2)
        {
            Matrix re = new Matrix(m1.row + m2.row, m1.column);
            // use memory copying
            Buffer.BlockCopy(m1._matrix, 0, re._matrix, 0, m1.row * m1.column * sizeof(Double));
            Buffer.BlockCopy(m2._matrix, 0, re._matrix, m1.row * m1.column * sizeof(Double), m2.row * m2.column * sizeof(Double));
            return re;
        }

        public static Matrix ConcatenateColumn(Matrix m1, Matrix m2)
        {
            Matrix re = new Matrix(m1.row, m1.column + m2.column);
            // can't use memory copying
            // actually, can
            for (int i = 0; i < m1.row; i++)
            {
                Buffer.BlockCopy(m1._matrix, i * m1.column * sizeof(Double), re._matrix, i * re.column * sizeof(Double), m1.column * sizeof(Double));
                Buffer.BlockCopy(m2._matrix, i * m2.column * sizeof(Double), re._matrix, (i * re.column + m1.column) * sizeof(Double), m2.column * sizeof(Double));
            }
            return re;
        }

        public void ImportRows(Matrix m2, int offset)
        {
            // use memory copying
            Buffer.BlockCopy(m2._matrix, 0, this._matrix, offset * this.column * sizeof(Double), m2.row * m2.column * sizeof(Double));
        }

        #endregion

        #region Advanced matrix functions
        // Decompose matrix m into lower triangular matrix L and strict upper triangular matrix U
        public static void Decompose(Matrix m, out Matrix L, out Matrix U)
        {
            int size = m.Height;
            if (!m.isSquare)
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

        public static Boolean SomeClose(Matrix m1, Matrix m2, Double maxDiff, int startRow)
        {
            if (m1.Width != m2.Width)
            {
                Exception e = new Exception("Two matrixes must have the same width!");
                throw e;
            }
            if (m1.Height > m2.Height)
            {
                Exception e = new Exception("Matrix m2's height must be equal to or greater than m1's!");
                throw e;
            }

            for (int i = 0; i < m1.Height; i++)
                for (int j = 0; j < m1.Width; j++)
                    if (Math.Abs(m1[i, j] - m2[startRow + i, j]) > Math.Abs(maxDiff))
                        return false;
            return true;
        }

        // Decompose matrix into L (and L*) using Cholesky algorithm
        public static Boolean CholeskyDecompose(Matrix A, out Matrix L)
        {
            // see http://en.wikipedia.org/wiki/Cholesky_decomposition

            L = Matrix.zeroLike(A);

            if (!A.isSquare) // only square matrix can be symmetric
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

        // check if this matrix is positive definite using definition and CholeskyDecompose
        public Boolean isPositiveDefinite
        {
            // In linear algebra, a symmetric n × n real matrix M is said to be positive definite if zTMz is positive for every non-zero column vector z of n real numbers. Here zT denotes the transpose of z.
            get
            {
                if (!this.isSquare) // only square matrix can be symmetric
                    return false;
                if (this.ToString() != Matrix.Transpose(this).ToString()) // matrix A is symmetric <=> A = AT
                    return false;

                // now check if it's positive definite

                // generate vector z with random number from -1 to 1. make sure z is non-zero
                Matrix z = Matrix.random(this.dim1, 1, -1, 1);
                while (z.isZero)
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
        }

        // check if this matrix is diagonally dominant
        public Boolean isDiagonallyDominant
        {
            // see https://en.wikipedia.org/wiki/Diagonally_dominant_matrix
            // In mathematics, a matrix is said to be diagonally dominant if for every row of the matrix,
            // the magnitude of the diagonal entry in a row is larger than or equal to the sum of the
            // magnitudes of all the other (non-diagonal) entries in that row.
            // More precisely, the matrix A is diagonally dominant if
            // | A[i,i] | >= sum(j!=i) of | A[i,j] |
            get
            {
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
        // see http://msdn.microsoft.com/en-us/magazine/jj863137.aspx
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
            if (!matrix.isSquare)
            {
                Exception e = new Exception("Matrix must be square!");
                throw e;
            }
            int size = matrix.Height;
            if (size == 1)
                return matrix[0, 0];
            else if (size == 2)
                return matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0];
            int[] perm;
            int toggle;
            Matrix lum = LUPDecompose(matrix, out perm, out toggle);
            if (lum == null)
                return 0; // throw new Exception("Unable to compute MatrixDeterminant");
            return Determinant(lum, perm, toggle);
        }

        public static double Determinant(Matrix lum, int[] perm, int toggle)
        {
            if (lum == null)
                return 0;
            double result = toggle;
            for (int i = 0; i < lum.Width; ++i)
                result *= lum[i, i];
            return result;
        }

        public static double[] HelperSolve(Matrix luMatrix, double[] b)
        {
            // solve luMatrix * x = b
            int n = luMatrix.dim1;
            double[] x = new double[n];
            b.CopyTo(x, 0);
            for (int i = 1; i < n; ++i)
            {
                double sum = x[i];
                for (int j = 0; j < i; ++j)
                    sum -= luMatrix[i, j] * x[j];
                x[i] = sum;
            }
            x[n - 1] /= luMatrix[n - 1, n - 1];
            for (int i = n - 2; i >= 0; --i)
            {
                double sum = x[i];
                for (int j = i + 1; j < n; ++j)
                    sum -= luMatrix[i, j] * x[j];
                x[i] = sum / luMatrix[i, i];
            }
            return x;
        }

        public static Matrix InverseAlt(Matrix matrix)
        {
            if (!matrix.isSquare)
            {
                Exception e = new Exception("Matrix must be square!");
                throw e;
            }

            int n = matrix.dim1;
            Matrix result = Matrix.Duplicate(matrix);
            int[] perm;
            int toggle;
            Matrix lum = LUPDecompose(matrix, out perm, out toggle);
            if (lum == null)
                return Matrix.zeroLike(matrix); //throw new Exception("Unable to compute inverse");

            Double det = Determinant(lum, perm, toggle);
            if (det == 0) // not invertible
            {
                // still return for the sake of simplicity
                // Zero matrix * any matrix = zero matrix
                // so it's never a valid answer
                return Matrix.zeroLike(matrix);
            }

            double[] b = new double[n];
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    if (i == perm[j])
                        b[j] = 1.0;
                    else
                        b[j] = 0.0;
                }
                double[] x = HelperSolve(lum, b);
                for (int j = 0; j < n; ++j)
                    result[j, i] = x[j];
            }
            return result;
        }

        #endregion

        #region Rouding
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

        #region Overriding methods
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
            StringBuilder re = new StringBuilder();
            for (int i = 0; i < this.Height; i++)
            {
                for (int j = 0; j < this.Width; j++)
                {
                    re.Append(Format(this[i, j]) + " ");
                }
                if (i < this.Height - 1)
                    re.Append("\n");
            }
            return re.ToString();
        }

        public string ToString(double rouding)
        {
            StringBuilder re = new StringBuilder();
            for (int i = 0; i < this.Height; i++)
            {
                for (int j = 0; j < this.Width; j++)
                {
                    re.Append(Format(RoundNum(this[i, j], rouding)) + " ");
                }
                if (i < this.Height - 1)
                    re.Append("\n");
            }
            return re.ToString();
        }

        private string Format(Double n)
        {
            return String.Format("{0:0.###############}", n);
        }
        #endregion
    }
}
