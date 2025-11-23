# Three-dimensional exact-cover problem solvers and generator
The repository contains two separate solvers and one generator for solving and generating three-dimensional exact-cover problem puzzles.\
solver::Solver and generator::Generator were produced out of my tiny little brain, and are thus quite unperformant compared to the bleeding edge methods for solving exact cover problems.\
\
I wrote these on my own intuition on how such problems could be solved, as a result of which they are quite brute-forcey.\
Solver has gone through four iterations now, the first three of which significantly improved it's performance. The fourth iteration actually tanked it from the third, likely due to changes to the Block struct making it heavier to clone, 
although the readability was significantly improved.\
\
The latest thing I wrote is the first thing in this repository that I have based on actual research, and it utilizes methods known as "dancing links" and "algorithm x" known together as "DLX".\
The performance of this based on very small initial testing is about 200 times faster than my own solver, as it does not require the copying of any data during the solve nor precomputing massive test tables of every possible combination.\

## NOTE! solver::Solver and generator::Generators functionality may be currently broken as the common struct implementations were modified during the writing of dlx::DLX
