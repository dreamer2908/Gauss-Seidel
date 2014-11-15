using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;

namespace Gauss_Seidel_Sequential.Test
{
    [TestFixture]
    class UtilsTest
    {
        [Test]
        public void parseInput12_Sample1_ChecksThem()
        {
            string sample = "4\n10 -1 2 0\n-1 11 -1 3\n2 -1 10 -1\n0 3 -1 8\n6 25 -11 15\n1 2 -1 1\n-------------";
            Matrix A, b, sol;
            Boolean re = Utils.parseInput(sample, out A, out b, out sol);
            try
            {
                Console.WriteLine("Matrix A:");
                Console.WriteLine(A.ToString());
            }
            catch (Exception)
            {
            }
            try
            {
                Console.WriteLine("Matrix b:");
                Console.WriteLine(b.ToString());
            }
            catch (Exception)
            {
            }
            try
            {
                Console.WriteLine("Matrix sol:");
                Console.WriteLine(sol.ToString());
            }
            catch (Exception)
            {
            }

            Matrix _A = new Matrix(4, 4);
            _A[0, 0] = 10;
            _A[0, 1] = -1;
            _A[0, 2] = 2;
            _A[0, 3] = 0;
            _A[1, 0] = -1;
            _A[1, 1] = 11;
            _A[1, 2] = -1;
            _A[1, 3] = 3;
            _A[2, 0] = 2;
            _A[2, 1] = -1;
            _A[2, 2] = 10;
            _A[2, 3] = -1;
            _A[3, 0] = 0;
            _A[3, 1] = 3;
            _A[3, 2] = -1;
            _A[3, 3] = 8;
            Matrix _b = new Matrix(4, 1);
            _b[0, 0] = 6;
            _b[1, 0] = 25;
            _b[2, 0] = -11;
            _b[3, 0] = 15;
            Matrix _sol = new Matrix(4, 1);
            _sol[0, 0] = 1;
            _sol[1, 0] = 2;
            _sol[2, 0] = -1;
            _sol[3, 0] = 1;

            Assert.AreEqual(true, re);
            Assert.AreEqual(_A.ToString(), A.ToString());
            Assert.AreEqual(_b.ToString(), b.ToString());
            Assert.AreEqual(_sol.ToString(), sol.ToString());
        }

        [Test]
        public void parseInput34_Sample1_ChecksThem()
        {
            string sample = "4\n10 -1 2 0\n-1 11 -1 3\n2 -1 10 -1\n0 3 -1 8\n6 25 -11 15\n1 2 -1 1\n-------------\n";
            sample += "\n3\n1 3 5\n3 -5 1\n0 5 2\n1 9 3\n-------";
            sample += "\n2\n2 1\n3 4\n1 6\n4 4\n";
            List<Matrix> _A = new List<Matrix>(), _b = new List<Matrix>(), _sol = new List<Matrix>();
            Utils.parseInput(sample, out _A, out _b, out _sol);
            for (int i = 0; i < _A.Count; i++)
            {
                Console.WriteLine("System #" + i.ToString());
                Matrix A = _A[i], b = _b[i], sol = _sol[i];

                Console.WriteLine("Matrix A:");
                Console.WriteLine(A.ToString());
                Console.WriteLine("Matrix b:");
                Console.WriteLine(b.ToString());

                Console.WriteLine("Matrix sol:");
                try { Console.WriteLine(sol.ToString()); }
                catch (Exception) {Console.WriteLine("None.\n"); }
            }
        }

        [TestCase(100, 3)]
        [TestCase(999, 10)]
        [TestCase(999, 3)]
        [TestCase(991, 13)]
        [TestCase(5, 7)]
        public void splitJob_VariousInputs_ChecksThem(int jobs, int machines)
        {
            Matrix re = Utils.splitJob(jobs, machines);
            Console.WriteLine(re.ToString());
            Assert.LessOrEqual(re.maxValue, jobs);
            Assert.GreaterOrEqual(re.maxValue, 0);
            Assert.LessOrEqual(Math.Abs(re.maxValue - re.minValue), 1);
            Assert.AreEqual(jobs, re.totalValue);
        }
    }
}
