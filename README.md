# Simplex-LP-Solver
This is a python implementation of the simplex algorithm.
The program would take a valid linear programming problem, in form of:

max(cx+z) subjected to:
Ax = b

In the future, I will improve the program so that it will accept more variations of LP.

The program would produce the optimal solution, as well as the steps.



## input:
The program can be run with no or one text file as an argument, if no text tile is provided, the program will run with the provided default example.

To start the program, run the any one of the following commands:

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

## process
The program would first try to convert the system into [Canonical form](https://en.wikipedia.org/wiki/Canonical_form#:~:text=In%20mathematics%20and%20computer%20science,identified%20in%20a%20unique%20way). This requires user to choose columns from A to get a basis. (In the future I will make the program select a feasible basis itself if the user input nothing.) If the select columns can form a basis that produces a feasible basic solution for x, the program would ask for user to confirm.

After a feasible basis is selected, the program would start the iteration, one step at a time. The user can press enter to see the result of next step, until an optimal solution is found.
