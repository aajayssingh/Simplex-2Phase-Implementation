# Simplex-2Phase-Implementation
This project is a C++ implementation of simplex two phase Algorithm. For the standard linear programs of maximization type.

Problem Statement : 
Implemented Two Phase Simplex algorithm for maximization LP problem , the input is assumed
to be in standard form Ax ≤ b.
1
Input description :
Assume following maximization LP problem
Maximize z = - x1 - x2 + 0s1 + 0s2 - A1 - A2
subject to the constraints :
2x1 + x2 s1 + A1 = 4,
x1 + 7x2 s2 + A2 = 7,
and x1 , x2, s1, s2 , A1, A2 0.
converted to give matrix input :
Maximize z + x1 + x2 + 0s1 + 0s2 + A1 + A2
subject to the constraints :
2x1 + x2 s1 +0s2 + A1 + 0A2 = 4,
x1 + 7x2 +0s1 s2 + 0A1 + A2 = 7,
and x1 , x2, s1, s2 , A1, A2 0.
corresponding matrix in std input is :
1100110
2 1 -1 0 1 0 4
1 7 0 -1 0 1 7
we input above nxm matrix in inp-params.txt file where n = 7, m =3, unknowns = 2,
artificial vars = 2.
are inputs asked .
2
Output Description :
output is the tableau at each simplex iteration in
x1 x2 s1 s2 a1 a2
z
1
1
0 0 1
1
following form :
a1 2
1
-1 0 1
0
a2 1
7
0 -1 0
1
phase 1 as well as in phase 2. Tableau is in
rhs
0
4
7
• row Z has rhs as the optimal value of objective function at particular iteration.
• column 1 has current basis variables.
• in phase 1 if Z ≺ 0 , implies that no feasible solution for original LP.
1• if phase 1 has Z=0 and no artificial variable in Basis implies that optimal solution exist
and we proceed to phase 2 , removing all artificial variables.
• if phase 1 has Z=0 and atleast one artificial variable at zero level , implies optimal solution
may exist. So, we proceed to phase 2 , removing non basic artificial variable (setting
corresponding column 0)and ensuring that for artificial variables which are in basis min
ratio is 0.
3
Compilation steps :
• set the input tableau in standard form in inp-params.txt.
• g++ -std=C++11 main.cpp -o exec.
4
Run
./exec
