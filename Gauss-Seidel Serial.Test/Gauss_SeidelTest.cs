using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;

namespace Gauss_Seidel_Serial.Test
{
    [TestFixture]
    class Gauss_SeidelTest
    {
        [Test]
        public void solve_sample1InWiki_converges() 
        {
            Matrix A = new Matrix(2, 2);
            A[0, 0] = 16;
            A[0, 1] = 3;
            A[1, 0] = 7;
            A[1, 1] = -11;
            Matrix b = new Matrix(2, 1);
            b[0, 0] = 11;
            b[1, 0] = 13;
            Matrix re = Gauss_Seidel.solve(A, b);
            re.Round(0.0001);
            Console.WriteLine("What it returns:");
            Console.WriteLine(re.ToString());
            Matrix expected = new Matrix(2, 1);
            expected[0, 0] = 0.8122;
            expected[1, 0] = -0.6650;
            expected.Round(0.0001);
            Console.WriteLine("Correct solution:");
            Console.WriteLine(expected.ToString());
            Assert.AreEqual(expected.ToString(), re.ToString());
        }

        [Test]
        public void solve_sample2InWiki_notConverge()
        {
            Matrix A = new Matrix(2, 2);
            A[0, 0] = 2;
            A[0, 1] = 3;
            A[1, 0] = 5;
            A[1, 1] = 7;
            Matrix b = new Matrix(2, 1);
            b[0, 0] = 11;
            b[1, 0] = 13;
            Matrix re = Gauss_Seidel.solve(A, b);
            re.Round(0.0001);
            Console.WriteLine("What it returns:");
            Console.WriteLine(re.ToString());
            Matrix expected = new Matrix(2, 1);
            expected[0, 0] = -38;
            expected[1, 0] = 29;
            expected.Round(0.0001);
            Console.WriteLine("Correct solution:");
            Console.WriteLine(expected.ToString());
            Assert.AreNotEqual(expected.ToString(), re.ToString());
        }

        [Test]
        public void solve_sample3InWiki_converges()
        {
            Matrix A = new Matrix(4, 4);
            A[0, 0] = 10;
            A[0, 1] = -1;
            A[0, 2] = 2;
            A[0, 3] = 0;
            A[1, 0] = -1;
            A[1, 1] = 11;
            A[1, 2] = -1;
            A[1, 3] = 3;
            A[2, 0] = 2;
            A[2, 1] = -1;
            A[2, 2] = 10;
            A[2, 3] = -1;
            A[3, 0] = 0;
            A[3, 1] = 3;
            A[3, 2] = -1;
            A[3, 3] = 8;
            Matrix b = new Matrix(4, 1);
            b[0, 0] = 6;
            b[1, 0] = 25;
            b[2, 0] = -11;
            b[3, 0] = 15;
            Matrix re = Gauss_Seidel.solve(A, b);
            re.Round(0.0001);
            Console.WriteLine("What it returns:");
            Console.WriteLine(re.ToString());
            Matrix expected = new Matrix(4, 1);
            expected[0, 0] = 2;
            expected[1, 0] = 1;
            expected[2, 0] = -1;
            expected[3, 0] = 1;
            expected.Round(0.0001);
            Console.WriteLine("Correct solution:");
            Console.WriteLine(expected.ToString());
            Assert.AreNotEqual(expected.ToString(), re.ToString());
        }
    }
}
