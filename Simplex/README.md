# Simplex-LP-Solver
This is a python implementation of the simplex algorithm.
The program would take a valid linear programming problem, in form of:

max(cx+z) subjected to:
Ax = b

In the future, I will improve the program so that it will accept more variations of LP.

The program would produce the optimal solution, as well as the steps.



## input:
The program can be run with no or one text file as an argument, if no text tile is provided, the program will run with the provided default example.

To run it, try something like the following:

```
Python simplex.py (name of input file)
Python simplex.py 
```

The given file should consist of lines containing only integers or decimal numbers or fractions, separated by spaces or tab. Although the program provides a basic check of the numerical input, invalid input may still cause the program to crash. (So please be nice to my program 'v'.)

The input file should be like this:

<pre>
m n
A<sub>11</sub> ... A<sub>1n</sub> 
. .     .
.   .   .
.     . .
A<sub>m1</sub> ... A<sub>mn</sub> 
b<sub>1</sub> ... b<sub>m</sub>  
c<sub>1</sub> ... c<sub>n</sub>  
z
</pre>

where m is the number of rows of A (i.e. constraints) and n is the number of columns (i.e. unknown);

A<sub>ij</sub> is the ij-entry of A;

b<sub>i</sub> is the i<sup>th</sup> element of b;

c<sub>i</sub> is the i<sup>th</sup> element of c;

z is a constant.

## functions:

#
```
readIn() -> LP
```

This function would convert the input to an LP, which consists of [m, n], A, b, c, z

#
```
simplex(LP, init_B = [], mode = 'a')
```
This function, by default, would produce the result of the LP (optimal cx and x if they exist). Set the mode = 's' to see detailed step or 'r' to see the basic result.
