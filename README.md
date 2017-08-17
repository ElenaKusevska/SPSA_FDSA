# SPSA
Fortran implementation of the SPSA and FDSA algorithms

This is a Fortran 90/95 program that implements the SPSA (Simultaneous Perturbation Stochastic Approximation) and FDSA (Finite 
Difference Stochastic Approximation algorithm). 

At the moment I have two functions to choose from, but the code can be applied to any function. I have so far only used it for 
functions with 2 variables, because I want to look at the output on a surface plot. But, the code should be able to handle 
functions of any number of variables without a problem. 

The program can deal with filenames ending in *00-*99, so that 101 different runs, with different choice of coefficients or 
initial guess can be stored in one folder for comparison. 

It takes in two input files:
  guessXX (where XX = 00-99, is the number of the run) - that contains the initial guess, and
  coeffXX (where XX = 00-99, is the number of the run) - that contains the gain coefficients, and other details for the particulaR 
    run of the program. For example, here, the procedure (spsa or fdsa) can be specified, the random number distribution to be used 
    (currently I have only implemented the Bernoulli +/-1 and the split uniform distribution), whether to define “a” for the spsa 
    algorithm in terms of the maximum number of initial steps or to just specify it a number, etc. 

It also creates three output files:
  run_detailsXX (where XX = 00-99, is the number of the run) - contains the gradient, the final number of iterations, and the 
    energy difference between the two final steps
  fn_XX (where XX = 00-99, is the number of the run, and n is the function being optimized – f1, f2…) - contains the Cartesian 
    coordinates of every step in the optimization procedure, and can be used for plotting, and
  fn_XX_last (where XX = 00-99, is the number of the run, and n is the function being optimized – f1, f2…) - contains the 
    Cartesian coordinates of the last step in the optimization procedure, and can be used for plotting

I use gnuplot to plot the results of the optimization. There is a folder included with the input and output for a couple of runs 
of the program, as well as the resulting plots and the gnuplot scripts that I have used. For the plots.

---------------------------------------------------------------------------------------------
The details that can be specified in the coeffXX input file:
---------------------------------------------------------------------------------------------

Line 1: a
a is… 

line 2: c...

