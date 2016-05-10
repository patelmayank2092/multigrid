# Simulation and Scientific Computing 2 
#Assignment 1
#Group: Sachin Sanil, Shivang Vyas, Mayank Patel, Vinayak Gholap

To run the code follow the steps. 
1. Open the terminal in the code directory.. 
2. First type "make clean" (without the inverted commas) to remove all the existing object files (.o) and the text files (.txt).
3. Then to compile all the files type "make all" (without the inverted commas). This will compile all the .h and .cpp files. This is done by the Makefile. Appropriate flags are given in Makefile as required.
4. We now will have an executable file "mgsolve". 
5. Here to execute this file you will need to pass the number of levels (l) and number of V-cycles (n) for the Multigrid, type "./mgsolve l n".
6. Now you will have the Residual Norm's and Convergence rate's for each level of V-cycle printed on the terminal window. Below you will find the time taken by the Solver in seconds.
7. A text file with name "solution.txt" will be generated in the code directory / folder. This file contains the number of grid points in X & Y axes and the third column indicates U approximate after the final V-cycle as calculated by the Solver.
8. To get the Error Norm for a fixed number of V-cycle, do the step 5 with number of level varying from 3 to 8 and for a fixed V-cycle for all levels. This will generate a text file named "E-Norm.txt" with the number of grid points and Error norm in 2nd column.
9. A Error Norm graph is present named "E-Norm_plot.pdf".
10. Finally to run the code with other "l & n" inputs, type "make clean" to remove existing generated files and repeat the steps from step 3.

There are two extra plots for u approximate and the plot for u approximate and u exact for a multigrid.
