Gauss-Seidel
=====================

Yet another random university project. 

### Requirement ###

 - Microsoft Visual C# 2008 Express Edition or newer
 - MPI.NET Runtime, MPI.NET SDK
 - Microsoft Compute Cluster Pack

### CLI options ###

This applies to both parallel and sequential versions.
All options have a short and a long switch.

 - `--input inputfile.txt` (or `-i inputfile.txt`): import system(s) from inputfile.txt.
 - `--output outputfile.txt` (or `-o outputfile.txt`): write result(s) to outputfile.txt.
 - `--show-equation` (or `-e`): include the system in result, either print it to console or write to output file.
 - `--show-benchmark` (or `-m`): show how much time it spent on sequential part, parallel part, and communication.
 - `--generate-input` (or `-g`): generate a valid input for `--input` from system(s) and write it (instead of results) to `outputfile.txt` specified by `--output`.
 - `--benchmark 1000` (or `-b 1000`): generate a (bunch of) random system(s) with 1000 equations, 1000 unknowns. Change 1000 to whatever size you want. `--input` will be ignored if this is set.
 - `--times 3` (or `-t 3`): specify number of system to generate. Replace 3 with the number you want.

### Input file format ###

Input files should be in this format.

 - Size
 - A matrix (each row in a separated line, numbers in a row are separated by a space)
 - RHS vector (all in one line)
 - Expected result (all in one line, if needed)
 - A line of `-------------` between each input linear system

Example:

`4
10 -1 2 0
-1 11 -1 3
2 -1 10 -1
0 3 -1 8
6 25 -11 15
1 2 −1 1
-------------`

for

`10.0*x1 + -1.0*x2 + 2.0*x3 + 0.0*x4 = 6.0
-1.0*x1 + 11.0*x2 + -1.0*x3 + 3.0*x4 = 25.0
2.0*x1 + -1.0*x2 + 10.0*x3 + -1.0*x4 = -11.0
0.0*x1 + 3.0*x2 + -1.0*x3 + 8.0*x4 = 15.0`
The exact solution of the system is (1, 2, −1, 1).

Note:
 - Expected result line is optional and is actually not used.
 - There's no empty line between lines. See the raw file.