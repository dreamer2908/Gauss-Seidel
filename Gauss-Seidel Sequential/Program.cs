﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace Gauss_Seidel_Sequential
{
    class Program
    {
        static void Main(string[] _args)
        {
            string programName = "Gauss-Seidel Linear System of Equations Solver";
            string programVer = "1.0 (sequential)";
            string programAuthor = "Quy N.H.";
            Console.WriteLine(programName + " v" + programVer + " by " + programAuthor + "\n");

            bool testing = false;
            string[] args = _args;
            if (testing)
                args = "-o output.txt -b 200 -t 10".Split(new char[] { ' ' });
            // parse args
            string inputFile = "", outputFile = "";
            bool benchmarkMode = false, showEquation = false, generateInput = false, showBenchmark = false;
            int benchmarkSize = 3;
            int benchmarkTime = 1;
            int i = 0;
            while (i < args.Length)
            {
                string arg = args[i];
                if (arg.StartsWith("--"))
                {
                    arg = arg.Substring(2);
                    switch (arg)
                    {
                        case "input": if (i + 1 < args.Length) inputFile = args[i + 1]; break;
                        case "output": if (i + 1 < args.Length) outputFile = args[i + 1]; break;
                        case "show-equation": showEquation = true; break;
                        case "show-benchmark": showBenchmark = true; break;
                        case "generate-input": generateInput = true; break;
                        case "benchmark": if (i + 1 < args.Length && int.TryParse(args[i + 1], out benchmarkSize)) { benchmarkMode = true; i++; }; break;
                        case "times": if (i + 1 < args.Length && int.TryParse(args[i + 1], out benchmarkTime)) { benchmarkMode = true; i++; }; break;
                    }
                }
                else if (arg.StartsWith("-"))
                {
                    arg = arg.Substring(1);
                    switch (arg)
                    {
                        case "i": if (i + 1 < args.Length) inputFile = args[i + 1]; break;
                        case "o": if (i + 1 < args.Length) outputFile = args[i + 1]; break;
                        case "e": showEquation = true; break;
                        case "m": showBenchmark = true; break;
                        case "g": generateInput = true; break;
                        case "b": if (i + 1 < args.Length && int.TryParse(args[i + 1], out benchmarkSize)) { benchmarkMode = true; i++; }; break;
                        case "t": if (i + 1 < args.Length && int.TryParse(args[i + 1], out benchmarkTime)) { benchmarkMode = true; i++; }; break;
                    }
                }
                i++;
            }

            // get input(s)
            List<Matrix> As = new List<Matrix>(), bs = new List<Matrix>(), sols = new List<Matrix>(), xs = new List<Matrix>();

            if (benchmarkMode)
            {
                // generate input
                for (int j = 0; j < benchmarkTime; j++)
                {
                    As.Add(Matrix.generateDiagonallyDominantMatrix(benchmarkSize, true, -100, 100));
                    bs.Add(Matrix.random(benchmarkSize, 1, -100, 100, true));
                }
                Console.WriteLine("Generated " + benchmarkTime.ToString() + " random system(s) to solve.");
            }
            else if (inputFile.Length > 0 && File.Exists(inputFile))
            {
                // parse input
                string inputArray = File.ReadAllText(inputFile);
                Utils.parseInput(inputArray, out As, out bs, out sols);
                Console.WriteLine("Got " + As.Count.ToString() + " system(s) from input file.");
            }
            else
            {
                // yell at user
                Console.WriteLine("Give me some inputs!");
                Console.ReadKey();
                Environment.Exit(1);
            }

            // do the calculation
            List<bool> converges = new List<bool>();
            List<int> loopses = new List<int>();
            List<Matrix> errs = new List<Matrix>();

            int equCounts = As.Count;
            benchmark bm = new benchmark();
            string bmResult = "";

            Gauss_Seidel.showBenchmark = showBenchmark;

            bm.start();
            for (int j = 0; j < equCounts; j++)
            {
                Console.Write("Solving system #" + (j + 1).ToString() + "... ");
                Matrix x , err;
                int loops = 0;
                bool converge = Gauss_Seidel.solve(As[j], bs[j], out x, out err, out loops);
                xs.Add(x);
                loopses.Add(loops);
                converges.Add(converge);
                errs.Add(err);
                Console.WriteLine("Done.");
            }
            bmResult = bm.getResult();

            // write output
            if (!generateInput)
            {
                // show the result as usual
                writeOutput(outputFile, "\n");
                for (int j = 0; j < equCounts; j++)
                {
                    Matrix x = xs[j], err = errs[j];
                    int loops = loopses[j];
                    bool converge = converges[j];
                    string strResult = "";
                    if (showEquation)
                        strResult += "\nEquation:\n" + Utils.writeEquation(As[j], bs[j]);
                    strResult += "\nNo. equations: " + x.Height.ToString();
                    strResult += "\nSolution: " + Matrix.Transpose(x).ToString(1e-14);
                    strResult += "\nErrors: " + Matrix.Transpose(err).ToString(1e-14);
                    strResult += "\nMean absolute error: " + string.Format("{0:0.##############}", Matrix.Abs(err).avgValue);
                    strResult += "\nConverged: " + converge.ToString();
                    strResult += "\nLoops: " + loops.ToString();
                    writeOutput(outputFile, strResult);
                }
                writeOutput(outputFile, "\nElapsed time: " + bmResult + " (" + string.Format("{0:0.###}", bm.getElapsedSeconds() / equCounts) + " sec / equation).");
                writeOutput(outputFile, "");
            }
            else
            {
                // create a valid input file
                for (int j = 0; j < equCounts; j++)
                {
                    Matrix x = xs[j], err = errs[j];
                    int loops = loopses[j];
                    bool converge = converges[j];
                    string strResult = "\n-----------\n";
                    strResult += x.Height.ToString();
                    strResult += "\n";
                    strResult += As[j].ToString();
                    strResult += "\n";
                    strResult += Matrix.Transpose(bs[j]).ToString();
                    strResult += "\n";
                    strResult += Matrix.Transpose(x).ToString(1e-14);
                    strResult += "\n";
                    writeOutput(outputFile, strResult);
                }
                writeOutput("", "\nElapsed time: " + bmResult + " (" + string.Format("{0:0.###}", bm.getElapsedSeconds() / equCounts) + " sec / equation).");
                writeOutput("", "");
            }

            Console.WriteLine("Done. Press a key to exit...");
            Console.ReadKey();
        }

        private static void writeOutput(string outputFile, string strResult)
        {
            if (outputFile != "")
            {
                File.AppendAllText(outputFile, strResult.Replace("\n", "\r\n"));
            }
            else
            {
                Console.WriteLine(strResult);
            }
        }
    }
}
