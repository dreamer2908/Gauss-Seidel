using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Gauss_Seidel_Serial
{
    class Utils
    {
        public static Boolean parseInput(string input, out Matrix A, out Matrix b, out Matrix sol)
        {
            // split input into lines. All empty lines are removed
            string[] inputArray = input.Split(new char[] {'\n', '\r'},  StringSplitOptions.RemoveEmptyEntries);
            return parseInput(inputArray, out A, out b, out sol);
        }

        public static bool parseInput(string[] inputArray, out Matrix A, out Matrix b, out Matrix sol)
        {
            // assume its format is correct FOR NOW
            int size = 0;
            if (inputArray.Length > 3 && int.TryParse(inputArray[0], out size))
            {
                A = new Matrix(size, size); // size of A = number of entries on each line
                for (int i = 0; i < size; i++)
                {
                    string line = inputArray[i + 1];
                    string[] entries = line.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                    for (int j = 0; j < size; j++)
                    {
                        Double number;
                        if (Double.TryParse(entries[j], out number))
                            A[i, j] = number;
                        else
                            A[i, j] = 0;
                    }
                }
                // now RHS vector
                string RHSLine = inputArray[size + 1];
                string[] RHSEntries = RHSLine.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                b = new Matrix(size, 1);
                for (int i = 0; i < size; i++)
                {
                    Double number;
                    if (Double.TryParse(RHSEntries[i], out number))
                        b[i, 0] = number;
                    else
                        b[i, 0] = 0;
                }
                // now solution vector
                if (inputArray.Length > size + 2)
                {
                    string solLine = inputArray[size + 2];
                    if (!solLine.StartsWith("--"))
                    {
                        string[] solEntries = solLine.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        sol = new Matrix(size, 1);
                        for (int i = 0; i < size; i++)
                        {
                            Double number;
                            if (Double.TryParse(solEntries[i], out number))
                                sol[i, 0] = number;
                            else
                                sol[i, 0] = 0;
                        }
                    }
                    else
                    {
                        sol = null;
                    }
                }
                else
                {
                    sol = null;
                }
                
                return true;
            }
            else
            {
                A = null;
                b = null;
                sol = null;
                return false;
            }
        }

        public static void parseInput(string input, out List<Matrix> A, out List<Matrix> b, out List<Matrix> sol)
        {
            // split input into lines. All empty lines are removed
            string[] inputArray = input.Split(new char[] { '\n', '\r' }, StringSplitOptions.RemoveEmptyEntries);
            parseInput(inputArray, out A, out b, out sol);
        }

        public static void parseInput(string[] inputArray, out List<Matrix> A, out List<Matrix> b, out List<Matrix> sol)
        {
            A = new List<Matrix>();
            b = new List<Matrix>();
            sol = new List<Matrix>();
            // split input into different groups (each is supposed to be a system of equations) and parse them separatedly
            int separatorPos = -1;
            for (int i = 0; i < inputArray.Length; i++)
            {
                if (inputArray[i].StartsWith("--") || i == inputArray.Length - 1)
                {
                    if (i - separatorPos - 1 > 2) // at least 2 lines between them
                    {
                        string[] _inputArray = new string[i - separatorPos - 1];
                        Array.Copy(inputArray, separatorPos + 1, _inputArray, 0, i - separatorPos - 1);
                        Matrix _A, _b, _sol;
                        //Console.WriteLine();
                        //Console.WriteLine(string.Join("\n", _inputArray));
                        parseInput(_inputArray, out _A, out _b, out _sol);
                        A.Add(_A);
                        b.Add(_b);
                        sol.Add(_sol);
                    }
                    separatorPos = i;
                }
            }
        }

        public static Matrix splitJob(int jobs, int machines)
        {
            Matrix re = new Matrix(1, machines);
            if (jobs >= machines)
            {
                re.numFill((int)Math.Floor(1.0 * jobs / machines));
                if (jobs % machines != 0) // there're still some jobs unassigned
                {
                    double remainJobsForEach = 1.0 * (jobs % machines) / machines;
                    int groupSize = (int)(1.0 / remainJobsForEach); // one machine from each group will be assigned one more job
                    for (int i = 0; i < machines; i++)
                    {
                        if (i % groupSize == 0)
                            re[0, i] += 1;
                        if (i == machines - 1) //last one takes the rest
                            re[0, i] += jobs - re.totalValue;
                    }
                }
            }
            else
            {
                for (int i = 0; i < jobs; i++)
                    re[0, i] = 1;
            }
            return re;
        }
    }
}
