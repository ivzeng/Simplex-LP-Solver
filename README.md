# Simplex-LP-Solver
This is a python implementation of the simplex algorithm.
The program would take a valid linear programming problem, in form of:

max(cx+z) subjected to:
Ax = b

In the future, I will improve the program so that it will accept more variations of LP.

## input:
The program can be run with no or one text file as an argument, if no text tile is provided, the program will run with the provided example.

The given file should consist of lines containing only integers (or maybe decimal numbers and fractions in the future version), separated by spaces. Although the program provides a basic check of the numerical input,
invalid input may still cause the program to crash. (So please be nice to my program 'v'.)

The text file should be like this:

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
