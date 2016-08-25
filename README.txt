Directory Contents
1) gov_prob.cpp  - solves the governement problem.  Returns bond pricing and bond choice schedule. Default decision is an input. 
2) gov_prob_guess.cpp - solves the government problem where default decision is by the governement.
3) kss_dist_d.cpp - creates a distribution of wealth under default.
4) kss_dist_nd.cpp - creaets a distribution of wealth under no default.
5) kss_policy.cpp - Solves the consumer problem using a linear parametrization of aggregate capital. 
6) kss_policy2.cpp - Solves the consumer problem using a lookup table for aggregate capital.
7) findSS.m  - first script to run.  Creates mexPar.mat binaries, which are needed to run run_debt.m.  The script finds the steady state of aggregate capital for the economy.
8) plot2.m
9) run_debt.m - With compile flag turned on, will compile the .cpp files to mex.  With guess flag turned on, will use gov_prob_guess to get a guess for the pricing schedule and kss_policy, which uses linear parametrization of aggregate capital.  The full problem uses gov_poicy2 which uses a lookup table for aggregate capital.  
10) markovappr.m - script returns a markov transition matrix and the corresponding indices for an AR(1) method using Tauchen method.
11) prodpop.mat - binaries holding population and wage/age profile used in other scripts.
12) betahat.mat - binaries which hold guesses for parametrization of aggregate capital.


