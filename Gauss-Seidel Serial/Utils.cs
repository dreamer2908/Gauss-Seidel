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
            // assume its format is correct FOR NOW
            int offset = 0; // loop later
            int size = 0;
            if (int.TryParse(inputArray[offset], out size))
            {
                A = new Matrix(size, size); // size of A = number of entries on each line
                for (int i = 0; i < size; i++)
                {
                    string line = inputArray[offset + i + 1];
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
                string RHSLine = inputArray[offset + size + 1];
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
                string solLine = inputArray[offset + size + 2];
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
    }
}
