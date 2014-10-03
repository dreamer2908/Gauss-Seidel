using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;

namespace Gauss_Seidel_Serial.Test
{
    [TestFixture]
    public class MatrixTest
    {
        #region Tests for Multiply

        [TestCase(3, 3, 2, 5, 0, 0, 0, 0)]
        [TestCase(3, 3, 2, 5, 7, 4, 33, 34)]
        [TestCase(25, 22, 1, 23, 7, 17, 549, 398)]
        [TestCase(-5, 7, 9, 8, 2, -3, -31, -6)]
        public void Multiply_VariousInputs2x2x2x1_ChecksThem(int m1_00, int m1_01, int m1_10, int m1_11, int m2_00, int m2_10, int m3_00, int m3_10)
        {
            Matrix m1 = new Matrix(2, 2);
            m1[0, 0] = m1_00;
            m1[0, 1] = m1_01;
            m1[1, 0] = m1_10;
            m1[1, 1] = m1_11;
            Matrix m2 = new Matrix(2, 1);
            m2[0, 0] = m2_00;
            m2[1, 0] = m2_10;
            Matrix m3 = new Matrix(2, 1);
            m3[0, 0] = m3_00;
            m3[1, 0] = m3_10;
            Matrix re = m1 * m2;
            Assert.AreEqual(m3.ToString(), re.ToString());
        }

        [TestCase(-5, 8, 2, 7, 0, -3, -10, -59)]
        [TestCase(2, 566, 24, -7.5, 0.1, -3.7, 104.6, -2109.2)]
        [TestCase(2.6, -0.1, 6, -1.5, -0.1, -1.7, 15.61, -3.73)]
        [TestCase(0, -2.1, 13.1, 1.1, -0.01, 1.001, 0.021, -2.1021)]
        public void Multiply_VariousInputs1x2x2x2_ChecksThem(Double m1_00, Double m1_01, Double m2_00, Double m2_01, Double m2_10, Double m2_11, Double m3_00, Double m3_01)
        {
            Matrix m1 = new Matrix(1, 2);
            m1[0, 0] = m1_00;
            m1[0, 1] = m1_01;
            Matrix m2 = new Matrix(2, 2);
            m2[0, 0] = m2_00;
            m2[0, 1] = m2_01;
            m2[1, 0] = m2_10;
            m2[1, 1] = m2_11;
            Matrix m3 = new Matrix(1, 2);
            m3[0, 0] = m3_00;
            m3[0, 1] = m3_01;
            Matrix re = m1 * m2;
            Assert.AreEqual(m3.ToString(), re.ToString());
        }
        #endregion

        #region Tests for Add/Subtract

        [TestCase(12, 23, -33, 1, -2, 5, 4, 5, 10, 28, -29, 6)]
        [TestCase(3.3, 2.3, -3.3, 1, -1.2, -0.5, 0.4, -0.5, 2.1, 1.8, -2.9, 0.5)]
        public void Add_2x2Inputs_ChecksThem(Double m1_00, Double m1_01, Double m1_10, Double m1_11, Double m2_00, Double m2_01, Double m2_10, Double m2_11, Double m3_00, Double m3_01, Double m3_10, Double m3_11)
        {
            Matrix m1 = new Matrix(2, 2);
            m1[0, 0] = m1_00;
            m1[0, 1] = m1_01;
            m1[1, 0] = m1_10;
            m1[1, 1] = m1_11;
            Matrix m2 = new Matrix(2, 2);
            m2[0, 0] = m2_00;
            m2[0, 1] = m2_01;
            m2[1, 0] = m2_10;
            m2[1, 1] = m2_11;
            Matrix m3 = new Matrix(2, 2);
            m3[0, 0] = m3_00;
            m3[0, 1] = m3_01;
            m3[1, 0] = m3_10;
            m3[1, 1] = m3_11;
            Matrix re = m1 + m2;
            Assert.AreEqual(m3.ToString(), re.ToString());
        }
        
        [TestCase(12, 23, -33, 1, 2, -5, -4, -5, 10, 28, -29, 6)]
        [TestCase(3.3, 2.3, -3.3, 1, 1.2, 0.5, -0.4, 0.5, 2.1, 1.8, -2.9, 0.5)]
        public void Subtract_2x2Inputs_ChecksThem(Double m1_00, Double m1_01, Double m1_10, Double m1_11, Double m2_00, Double m2_01, Double m2_10, Double m2_11, Double m3_00, Double m3_01, Double m3_10, Double m3_11)
        {
            Matrix m1 = new Matrix(2, 2);
            m1[0, 0] = m1_00;
            m1[0, 1] = m1_01;
            m1[1, 0] = m1_10;
            m1[1, 1] = m1_11;
            Matrix m2 = new Matrix(2, 2);
            m2[0, 0] = m2_00;
            m2[0, 1] = m2_01;
            m2[1, 0] = m2_10;
            m2[1, 1] = m2_11;
            Matrix m3 = new Matrix(2, 2);
            m3[0, 0] = m3_00;
            m3[0, 1] = m3_01;
            m3[1, 0] = m3_10;
            m3[1, 1] = m3_11;
            Matrix re = m1 - m2;
            Assert.AreEqual(m3.ToString(), re.ToString());
        }

        #endregion

        #region Tests for Inverse

        [TestCase(5, 2, -7, -3)]
        [TestCase(5, 8, 17, 3)]
        [TestCase(5, 3, 5, 3)]
        public void Inverse_2x2Inputs_ChecksThem(int k1, int k2, int k3, int k4)
        {
            Matrix m1 = new Matrix(2, 2);
            m1[0, 0] = k1;
            m1[0, 1] = k2;
            m1[1, 0] = k3;
            m1[1, 1] = k4;
            Matrix unitM = Matrix.unit(2);
            Matrix zeroM = Matrix.zero(2);
            Matrix re = ~m1;
            Matrix product = m1 * re;
            product.Round(0.001);
            if (m1.inversible())
            {
                Assert.AreEqual(product.ToString(), unitM.ToString());
            }
            else // return zero matrix if it's not inversible
            {
                Assert.AreEqual(re.ToString(), zeroM.ToString());
            }
        }

        [TestCase(-2, 3, -1, 5, -1, 4, 4, -8, 2)]
        [TestCase(10, 0, -3, -2, -4, 1, 3, 0, 2)]
        [TestCase(2, -3, -2, -6, 3, 3, -2, -3, -2)]
        [TestCase(-4, 5, 2, -3, 4, 2, -1, 2, 5)]
        [TestCase(1, -3, -6, -1, 5, 5, -1, 6, 5)]
        [TestCase(17, 17, 5, 21, 18, 21, 2, 2, 19)]
        public void Inverse_3x3Inputs_ChecksThem(int m1_00, int m1_01, int m1_02, int m1_10, int m1_11, int m1_12, int m1_20, int m1_21, int m1_22)
        {
            Matrix m1 = new Matrix(3, 3);
            m1[0, 0] = m1_00;
            m1[0, 1] = m1_01;
            m1[0, 2] = m1_02;
            m1[1, 0] = m1_10;
            m1[1, 1] = m1_11;
            m1[1, 2] = m1_12;
            m1[2, 0] = m1_20;
            m1[2, 1] = m1_21;
            m1[2, 2] = m1_22;
            Matrix re = ~m1;
            // System.Windows.Forms.MessageBox.Show(re.ToString());
            Matrix unitM = Matrix.unit(3);
            Matrix zeroM = Matrix.zero(3);
            Matrix product = m1 * re;
            product.Round(0.001);
            if (m1.inversible())
            {
                Assert.AreEqual(product.ToString(), unitM.ToString());
            }
            else // return zero matrix if it's not inversible
            {
                Assert.AreEqual(re.ToString(), zeroM.ToString());
            }
        }

        [TestCase(17, 17, 5, 1, 21, 18, 21, 0, 2, 2, 19, 3, 4, 6, -1, 5)]
        [TestCase(3, -5, -15, 1, 0, 2, -9, 10, 2, 3, -19, 23, -4, 0, -3, 5)]
        [TestCase(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)]
        public void Inverse_4x4Inputs_ChecksThem(int m1_00, int m1_01, int m1_02, int m1_03, int m1_10, int m1_11, int m1_12, int m1_13, int m1_20, int m1_21, int m1_22, int m1_23, int m1_30, int m1_31, int m1_32, int m1_33)
        {
            Matrix m1 = new Matrix(4, 4);
            m1[0, 0] = m1_00;
            m1[0, 1] = m1_01;
            m1[0, 2] = m1_02;
            m1[0, 3] = m1_03;
            m1[1, 0] = m1_10;
            m1[1, 1] = m1_11;
            m1[1, 2] = m1_12;
            m1[1, 3] = m1_13;
            m1[2, 0] = m1_20;
            m1[2, 1] = m1_21;
            m1[2, 2] = m1_22;
            m1[2, 3] = m1_23;
            m1[3, 0] = m1_30;
            m1[3, 1] = m1_31;
            m1[3, 2] = m1_32;
            m1[3, 3] = m1_33;
            Matrix re = ~m1;
            // System.Windows.Forms.MessageBox.Show(re.ToString());
            Matrix unitM = Matrix.unit(4);
            Matrix zeroM = Matrix.zero(4);
            Matrix product = m1 * re;
            product.Round(0.001);
            if (m1.inversible())
            {
                Assert.AreEqual(product.ToString(), unitM.ToString());
            }
            else // return zero matrix if it's not inversible
            {
                Assert.AreEqual(re.ToString(), zeroM.ToString());
            }
        }

        [Test]
        public void InverseMatrix_WrongSize__ThrowsFluent()
        {
            Matrix m = new Matrix(3, 2);
            var ex = Assert.Catch<Exception>(() => Matrix.Inverse(m));
            Assert.That(ex.Message, Is.StringContaining("Matrix must be square"));
        }
        #endregion

        # region Tests for determinant
        [TestCase(1, 1, 1, 1)]
        [TestCase(5, 2, -7, -3)]
        [TestCase(5, 8, 17, 3)]
        [TestCase(5, 3, 5, 3)]
        public void determinant_2x2Inputs_ChecksThem(int m2_00, int m2_01, int m2_10, int m2_11)
        {
            Matrix m2 = new Matrix(2, 2);
            m2[0, 0] = m2_00;
            m2[0, 1] = m2_01;
            m2[1, 0] = m2_10;
            m2[1, 1] = m2_11;
            Double expected = m2[0, 0] * m2[1, 1] - m2[0, 1] * m2[1, 0];
            Assert.AreEqual(expected, m2.determinant());
        }

        [TestCase(-2, 2, -3, -1, 1, 3, 2, 0, -1, 18)]
        [TestCase(1, 2, 3, 2, 3, 1, 3, 1, 2, -18)]
        [TestCase(-2, 3, -1, 5, -1, 4, 4, -8, 2, -6)]
        [TestCase(10, 0, -3, -2, -4, 1, 3, 0, 2, -116)]
        [TestCase(2, -3, -2, -6, 3, 3, -2, -3, -2, 12)]
        [TestCase(-4, 5, 2, -3, 4, 2, -1, 2, 5, -3)]
        [TestCase(1, -3, -6, -1, 5, 5, -1, 6, 5, 1)]
        public void determinant_3x3Inputs_ChecksThem(int m1_00, int m1_01, int m1_02, int m1_10, int m1_11, int m1_12, int m1_20, int m1_21, int m1_22, int expected)
        {
            Matrix m1 = new Matrix(3, 3);
            m1[0, 0] = m1_00;
            m1[0, 1] = m1_01;
            m1[0, 2] = m1_02;
            m1[1, 0] = m1_10;
            m1[1, 1] = m1_11;
            m1[1, 2] = m1_12;
            m1[2, 0] = m1_20;
            m1[2, 1] = m1_21;
            m1[2, 2] = m1_22;
            Assert.AreEqual(expected, m1.determinant());
        }

        [TestCase(9, 3, 5, 1, -6, -9, 7, 2, -1, -8, 1, 3, 9, 3, 5, 0, -615)]
        [TestCase(0, -3, 25, 8, -1, -2, 9, -2, 11, -4, 1, 0, -5, 3, -5, 1, 0)]
        [TestCase(20, -3, 0, 18, -1, -12, 0, 2, 4, 4, -1, 12, -25, 13, 0, 21, 11107)]
        public void determinant_4x4Inputs_ChecksThem(int m1_00, int m1_01, int m1_02, int m1_03, int m1_10, int m1_11, int m1_12, int m1_13, int m1_20, int m1_21, int m1_22, int m1_23, int m1_30, int m1_31, int m1_32, int m1_33, int expected)
        {
            Matrix m1 = new Matrix(4, 4);
            m1[0, 0] = m1_00;
            m1[0, 1] = m1_01;
            m1[0, 2] = m1_02;
            m1[0, 3] = m1_03;
            m1[1, 0] = m1_10;
            m1[1, 1] = m1_11;
            m1[1, 2] = m1_12;
            m1[1, 3] = m1_13;
            m1[2, 0] = m1_20;
            m1[2, 1] = m1_21;
            m1[2, 2] = m1_22;
            m1[2, 3] = m1_23;
            m1[3, 0] = m1_30;
            m1[3, 1] = m1_31;
            m1[3, 2] = m1_32;
            m1[3, 3] = m1_33;
            Assert.AreEqual(expected, m1.determinant());
        }
        #endregion

        #region Tests for Round

        [TestCase(0.05, 1.135, -3.04, -6.12, -1.4, 5.0001, 5.89, -1.21, 6.096, 5.01, 1.15, -3.05, -6.10, -1.4, 5.00, 5.90, -1.20, 6.10, 5.00)]
        public void Round_3x3Inputs_ChecksThem(Double rounding, Double m1_00, Double m1_01, Double m1_02, Double m1_10, Double m1_11, Double m1_12, Double m1_20, Double m1_21, Double m1_22, Double m2_00, Double m2_01, Double m2_02, Double m2_10, Double m2_11, Double m2_12, Double m2_20, Double m2_21, Double m2_22)
        {
            Matrix m1 = new Matrix(3, 3);
            m1[0, 0] = m1_00;
            m1[0, 1] = m1_01;
            m1[0, 2] = m1_02;
            m1[1, 0] = m1_10;
            m1[1, 1] = m1_11;
            m1[1, 2] = m1_12;
            m1[2, 0] = m1_20;
            m1[2, 1] = m1_21;
            m1[2, 2] = m1_22;
            Matrix m2 = new Matrix(3, 3);
            m2[0, 0] = m2_00;
            m2[0, 1] = m2_01;
            m2[0, 2] = m2_02;
            m2[1, 0] = m2_10;
            m2[1, 1] = m2_11;
            m2[1, 2] = m2_12;
            m2[2, 0] = m2_20;
            m2[2, 1] = m2_21;
            m2[2, 2] = m2_22;
            m1.Round(rounding);
            Assert.AreEqual(m1.ToString(), m2.ToString());
        }

        #endregion

        #region Tests for Transpose

        [Test]
        public void Transpose_Sample1InWikipedia_ChecksIt()
        {
            Matrix A = new Matrix(1, 2);
            A[0, 0] = 1;
            A[0, 1] = 2;
            Matrix expected = new Matrix(2, 1);
            expected[0, 0] = 1;
            expected[1, 0] = 2;
            Matrix re = Matrix.Transpose(A);
            Assert.AreEqual(expected.ToString(), re.ToString());
        }

        [Test]
        public void Transpose_Sample2InWikipedia_ChecksIt()
        {
            Matrix A = new Matrix(2, 2);
            A[0, 0] = 1;
            A[0, 1] = 2;
            A[1, 0] = 3;
            A[1, 1] = 4;
            Matrix expected = new Matrix(2, 2);
            expected[0, 0] = 1;
            expected[0, 1] = 3;
            expected[1, 0] = 2;
            expected[1, 1] = 4;
            Matrix re = Matrix.Transpose(A);
            Assert.AreEqual(expected.ToString(), re.ToString());
        }

        [Test]
        public void Transpose_Sample3InWikipedia_ChecksIt()
        {
            Matrix A = new Matrix(3, 2);
            A[0, 0] = 1;
            A[0, 1] = 2;
            A[1, 0] = 3;
            A[1, 1] = 4;
            A[2, 0] = 5;
            A[2, 1] = 6;
            Matrix expected = new Matrix(2, 3);
            expected[0, 0] = 1;
            expected[0, 1] = 3;
            expected[0, 2] = 5;
            expected[1, 0] = 2;
            expected[1, 1] = 4;
            expected[1, 2] = 6;
            Matrix re = Matrix.Transpose(A);
            Assert.AreEqual(expected.ToString(), re.ToString());
        }

        #endregion
    }
}
