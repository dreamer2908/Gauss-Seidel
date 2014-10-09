using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using Gauss_Seidel_Serial;
using MPI;

namespace Gauss_Seidel_Parallel
{
    class Program
    {
        static void Main(string[] _args)
        {
            using (new MPI.Environment(ref _args))
            {
                Intracommunicator comm = Communicator.world;
                if (comm.Rank == 0)
                {
                    // program for rank 0
                    string programName = "Gauss-Seidel Linear System of Equations Solver";
                    string programVer = "1.0 (parallel)";
                    string programAuthor = "Quy N.H.";
                    Console.WriteLine(programName + " v" + programVer + " by " + programAuthor + "\n");

                    // check number of processes in this communicator
                    if (comm.Size < 2)
                    {
                        Console.WriteLine("Please run at least 2 processes of me.");
                        Console.WriteLine("Exiting...");
                        MPI.Environment.Abort(1);
                    }

                    bool testing = true;
                    string[] args = _args;
                    if (testing)
                        args = "-o output.txt -b 200 -t 10".Split(new char[] { ' ' });
                    // parse args
                    string inputFile = "", outputFile = "";
                    bool benchmarkMode = false, showEquation = false;
                    int benchmarkSize = 3;
                    int benchmarkTime = 1;
                    int i = 0;
                    while (i < args.Length)
                    {
                        string arg = args[i];
                        if (arg.StartsWith("-"))
                        {
                            arg = arg.Substring(1);
                            switch (arg)
                            {
                                case "i": if (i + 1 < args.Length) inputFile = args[i + 1]; break;
                                case "o": if (i + 1 < args.Length) outputFile = args[i + 1]; break;
                                case "s": showEquation = true; break;
                                case "b": if (i + 1 < args.Length && int.TryParse(args[i + 1], out benchmarkSize)) { benchmarkMode = true; i++; }; break;
                                case "t": if (i + 1 < args.Length && int.TryParse(args[i + 1], out benchmarkTime)) { benchmarkMode = true; i++; }; break;
                            }
                        }
                        else if (arg.StartsWith("--"))
                        {
                            arg = arg.Substring(2);
                            switch (arg)
                            {
                                case "input": if (i + 1 < args.Length) inputFile = args[i + 1]; break;
                                case "output": if (i + 1 < args.Length) outputFile = args[i + 1]; break;
                                case "show-equation": showEquation = true; break;
                                case "benchmark": if (i + 1 < args.Length && int.TryParse(args[i + 1], out benchmarkSize)) { benchmarkMode = true; i++; }; break;
                                case "times": if (i + 1 < args.Length && int.TryParse(args[i + 1], out benchmarkTime)) { benchmarkMode = true; i++; }; break;
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
                    }
                    else if (inputFile.Length > 0)
                    {
                        // parse input
                        string[] inputArray = File.ReadAllLines(inputFile);
                        Utils.parseInput(inputArray, out As, out bs, out sols);
                    }
                    else
                    {
                        // yell at user
                        Console.WriteLine("Give me some inputs!");
                        Console.WriteLine("Exiting...");
                        MPI.Environment.Abort(1);
                    }

                    // do the calculation
                    List<bool> converges = new List<bool>();
                    List<int> loopses = new List<int>();
                    List<Matrix> errs = new List<Matrix>();

                    int equCounts = As.Count;
                    benchmark bm = new benchmark();
                    string bmResult = "";

                    bm.start();
                    for (int j = 0; j < equCounts; j++)
                    {
                        Console.WriteLine("Solving system equation #" + (j + 1).ToString());
                        Matrix x, err;
                        int loops = 0;
                        for (int r = 1; r < comm.Size; r++)
                        {
                            comm.Send("start", r, 0);
                        }
                        bool converge = Gauss_Seidel.solve(As[j], bs[j], out x, out err, out loops, ref comm);
                        xs.Add(x);
                        loopses.Add(loops);
                        converges.Add(converge);
                        errs.Add(err);
                    }
                    bmResult = bm.getResult();

                    // write output
                    for (int j = 0; j < equCounts; j++)
                    {
                        Matrix x = xs[j], err = errs[j];
                        int loops = loopses[j];
                        bool converge = converges[j];
                        string strResult = "";
                        if (showEquation)
                            strResult += "\nEquation:\n" + writeEquation(As[j], bs[j]);
                        strResult += "\nNo. equations: " + x.Height.ToString();
                        strResult += "\nSolution: " + Matrix.Transpose(x).ToString(1e-16);
                        strResult += "\nErrors: " + Matrix.Transpose(err).ToString(1e-16);
                        strResult += "\nMean absolute error: " + string.Format("{0:0.################}", Matrix.Abs(err).avgValue);
                        strResult += "\nConverged: " + converge.ToString();
                        strResult += "\nLoops: " + loops.ToString();
                        writeOutput(outputFile, strResult);
                    }
                    writeOutput(outputFile, "\nTotal time: " + bmResult);
                    writeOutput(outputFile, "\nAvg time: " + string.Format("{0:0.###}", bm.getElapsedSeconds() / equCounts) + " second per equation.");
                    writeOutput(outputFile, "\n");

                    Console.WriteLine("Done. Exiting...");

                    // tell other ranks to exit
                    for (int r = 1; r < comm.Size; r++)
                    {
                        comm.Send("exit", r, 0);
                    }
                }
                else
                {
                    // program for all other ranks
                    // wait for command (start (solveSub), exit)
                    string command = null;
                    do
                    {
                        command = comm.Receive<string>(0, 0); // receive command from rank 0
                        if (command == "start")
                        {
                            Gauss_Seidel.solveSub(ref comm);
                        }
                    } while (command != "exit");
                }
            }
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

        private static string writeEquation(Matrix A, Matrix b)
        {
            // assume valid inputs
            StringBuilder re = new StringBuilder();
            for (int i = 0; i < A.Height; i++)
            {
                if (i > 0)
                    re.Append("\n");
                re.Append(string.Format("{0:0.###}", A[i, 0]) + " * x" + 1.ToString());
                for (int j = 1; j < A.Width; j++)
                {
                    re.Append(" + " + string.Format("{0:0.###}", A[i, j]) + " * x" + (j + 1).ToString());
                }
                re.Append(" = " + string.Format("{0:0.###}", b[i, 0]));
            }
            return re.ToString();
        }
    }
}
