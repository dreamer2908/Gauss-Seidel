using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using Gauss_Seidel_Serial;
using MPI;

namespace Gauss_Seidel_Parallel.Test
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
                    string programVer = "1.0 (test)";
                    string programAuthor = "Quy N.H.";
                    Console.WriteLine(programName + " v" + programVer + " by " + programAuthor + "\n");

                    // check number of processes in this communicator
                    if (comm.Size < 2)
                    {
                        Console.WriteLine("Please run at least 2 processes of me.");
                        Console.WriteLine("Exiting...");
                        MPI.Environment.Abort(1);
                    }

                    string inputFile = "", outputFile = "";
                    bool benchmarkMode = true, showEquation = false, generateInput = false;
                    int benchmarkSize = 100;
                    int benchmarkTime = 100;

                    // get input(s)
                    List<Matrix> As = new List<Matrix>(), bs = new List<Matrix>(), sols = new List<Matrix>(), xs_p = new List<Matrix>(), xs_s = new List<Matrix>();

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
                        string inputArray = File.ReadAllText(inputFile);
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
                    List<bool> converges_p = new List<bool>(), converges_s = new List<bool>();
                    List<int> loopses_p = new List<int>(), loopses_s = new List<int>();
                    List<Matrix> errs_p = new List<Matrix>(), errs_s = new List<Matrix>();

                    int equCounts = As.Count;
                    benchmark bm = new benchmark();
                    string bmResult = "";

                    bm.start();
                    for (int j = 0; j < equCounts; j++)
                    {
                        Console.Write("Solving system #" + (j + 1).ToString() + "... ");
                        Matrix x, err;
                        int loops = 0;
                        for (int r = 1; r < comm.Size; r++)
                        {
                            comm.Send("start", r, 0);
                        }
                        Console.Write("Parallel... ");
                        bool converge = Gauss_Seidel_Parallel.solve(As[j], bs[j], out x, out err, out loops, comm);
                        xs_p.Add(x);
                        loopses_p.Add(loops);
                        converges_p.Add(converge);
                        errs_p.Add(err);
                        Console.Write("Serial... ");
                        converge = Gauss_Seidel.solve(As[j], bs[j], out x, out err, out loops);
                        xs_s.Add(x);
                        loopses_s.Add(loops);
                        converges_s.Add(converge);
                        errs_s.Add(err);
                        Console.WriteLine("Done.");
                    }
                    bmResult = bm.getResult();

                    // write output
                    Console.WriteLine("\nVerifying results:\n");
                    int total = 0, passed = 0, failed = 0;
                    for (int j = 0; j < equCounts; j++)
                    {
                        Matrix x_p = xs_p[j], err_p = errs_p[j], x_s = xs_s[j], err_s = errs_s[j];
                        int loops_p = loopses_p[j], loops_s = loopses_s[j];
                        bool converge_p = converges_p[j], converge_s = converges_s[j];
                        bool c = false, l = false, s = false;
                        Console.Write("System #" + (j + 1).ToString() + ": ");
                        if (s = x_p.ToString() == x_s.ToString())
                            Console.Write("solutions match, ");
                        else
                            Console.Write("solutions DON'T match, ");
                        if (l = loops_p == loops_s)
                            Console.Write("loop count matches, ");
                        else
                            Console.Write("loop count DOESN'T match, ");
                        if (c = converge_p == converge_s)
                            Console.Write("convergence matches. ");
                        else
                            Console.Write("convergence DOESN'T match. ");
                        if (s && c && l)
                        {
                            Console.WriteLine("Passed!");
                            passed += 1;
                        }
                        else
                        {
                            Console.WriteLine("Failed!");
                            failed += 1;
                        }
                        total += 1;
                    }

                    Console.WriteLine("\nTotal: " + total.ToString() + ". Passed: " + passed.ToString() + ". Failed: " + failed.ToString());

                    Console.WriteLine("\nElapsed time: " + bmResult + " (" + string.Format("{0:0.###}", bm.getElapsedSeconds() / equCounts) + " sec / equation).");

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
                            Matrix A = null, b = null, x, err; int loops;
                            Gauss_Seidel_Parallel.solve(A, b, out x, out err, out loops, comm);
                        }
                    } while (command != "exit");
                }
            }
        }
    }
}
