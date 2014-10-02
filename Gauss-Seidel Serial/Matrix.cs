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

        public void zeroFill()
        {
            for (int x = 0; x < dim1; x++)
                for (int y = 0; y < dim2; y++)
                    this[x, y] = 0;
        }

        public void randomFill()
        {
            for (int x = 0; x < dim1; x++)
                for (int y = 0; y < dim2; y++)
                    this[x, y] = _r.NextDouble();
        }

        public Boolean isSquare()
        {
            return (this.Height == this.Width);
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

        public static Boolean operator ==(Matrix m1, Matrix m2)
        {
            return (m1.ToString() == m2.ToString());
        }

        public static Boolean operator !=(Matrix m1, Matrix m2)
        {
            return (m1.ToString() != m2.ToString());
        }

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

        public static Matrix operator /(Matrix m1, Matrix m2)
        {
            return Divide(m1, m2);
        }

        public static Matrix operator /(Matrix m1, Double scalar)
        {
            return Divide(m1, scalar);
        }

        // Get the inverse matrix of m
        public static Matrix operator ~(Matrix m)
        {
            return Matrix.Inverse(m);
        }

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
            if (m.determinant() == 0) // not inversible 
            {
                // still return for the sake of simplicity
                // Zero matrix * any matrix = zero matrix
                // so it's never a valid answer
                return Matrix.zero(size);
            }

            Double det = m.determinant();

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

        // Return a separated copy of matrix m
        public static Matrix Duplicate(Matrix m)
        {
            Matrix re = new Matrix(m.Height, m.Width);
            for (int i = 0; i < re.Height; i++)
                for (int j = 0; j < re.Width; j++)
                    re[i, j] = m[i, j];
            return re;
        }

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
                re += "\n";
            }
            return re;
        }

        private string Format(Double n)
        {
            return String.Format("{0:0.###############}", n);
        }
    }
}
