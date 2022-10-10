#!/usr/bin/env python

from tabnanny import check
import numpy as np
from numpy.linalg import matrix_rank
from fractions import Fraction
import sys

#### linear programming, simplex method ####

# helper fns #

# read input file
def readIn():
	with open('defaultIn.txt' if len(sys.argv) == 1 else sys.argv[1], 'r') as input:
		dm = getList(input.readline())
		A = np.array([getList(input.readline(), type = 1) for i in range(dm[0])])
		b = np.array(getList(input.readline(), type = 1))
		c = np.array(getList(input.readline(), type = 1))
		z = getNum(input.readline(), type = 1)
	return dm, A, b, c, z

# compare vector, return [a,b,c]
#   a: # of ==
#   b: # of <
#   c: # of >
def vcmp(v1, v2):
	if len(v1) != len(v2):
		print('invalid comparing: vectors have different length')
		return -1
	count = [0,0,0]
	v1 -= v2
	for i in v1:
		if i == 0:
			count[0] += 1
		elif i < 0:
			count[1] += 1
		else:
			count[2] += 1
	return count

# count # of elements >(or >=, <, ...) 0
def cmp0(v, op):
	return sum(1 for i in v if op(i)) 

# get index of the first positive component, -1 if not exists:
def pstCpn(v):
	for i in range(len(v)):
		if v[i] > 0:
			return i
	return -1

# check input validity
def checkIn(dm, A, b, c):
	if dm[0] == 0 or dm[0] != len(b):
		print("invalid input: A's column and d has different lengths")
		return -1
	dm += [len(A[0])]
	if dm[1] == 0 or dm[1] != len(c):
		print("invalid input: c, x and A's row have different lengths")
		return -2
	return dm
# convert a floating point number into fraction
def NumToFrac(n):
	return str(Fraction(n).limit_denominator())

# print a Row
def printRow(r, wid):
	print('[', ' '.join([('{0: <'+str(wid)+'}').format(x) for x in [NumToFrac(i) for i in r]]), ']')

# convert string to a number, type = 0 if the number is int, otherwise the number is float
def getNum(n, type = 0):
	if type == 0:
		return int(n)
	n = [float(i) for i in n.split('/')]
	if len(n) == 1:
		return n[0]
	return n[0]/n[1]

# convert to an integer list
def getList(line, type = 0):
	return [getNum(i, type) for i in line.split()]

# print system
def printls(A, b, c, z, yt = 0,x = 0):
	print('max (c*x + z)')
	print('subjected to:')
	print('Ax = b')
	print('where:')
	print('A = [') 
	for r in A:
		printRow(r, 8)
	print(']')
	print('b = ')
	printRow(b, 8)
	if type(yt) != int:
		print('yt =')
		printRow(yt, 8)
	print('c = ')
	printRow(c, 8)
	print('z = ', Fraction(z).limit_denominator())
	if type(x) != int:
		print('basic solution: bx =')
		printRow(x, 8)

# convert the system into canonical form
def canonical(A, b, c, z, dm, idxs_in = []):
	def convertSystem(A_b):
		A_b_iv = np.linalg.inv(A_b)
		nb = np.matmul(A_b_iv, b)
		bx = getBasicSolution(nb, idxs_in, dm)
		nA = np.matmul(A_b_iv, A)
		yt = np.matmul(A_b_iv.T, c[idxs_in])
		nc = c - np.matmul(yt, A)
		nz = z + np.matmul(yt, b)
		return  nA, nb, nc, nz, yt, bx, idxs_in

	if len(idxs_in) != 0:
		idxs_in.sort()
		nA, nb, nc, nz, yt, bx, idxs_in = convertSystem(getBasis(A, dm, idxs_in))
		printls(nA, nb, nc, nz, yt, bx)
		return  nA, nb, nc, nz, bx, idxs_in
	done = False
	while not done:
		print('choose', dm[0], 'columns to form a basis')
		idxs_in = getList(input())
		idxs_in.sort()
		A_b = getBasis(A, dm, idxs_in)
		# error checking
		while type(A_b) is int:
			if A_b == -1:
				msg = "invalid indexes"
			elif A_b == -2:
				msg = "invalid length"
			else:
				msg = "columns are not linearly indepent" 
			print(msg + ', please enter again: \n')
			idxs_in = getList(input())
			idxs_in.sort()
			A_b = getBasis(A, dm, idxs_in)
			print('Basis selected:\n', A_b)

		nA, nb, nc, nz, yt, bx, idxs_in = convertSystem(A_b)
		# check feasibility of the basis
		if (cmp0(nb.round(5), lambda a: a < 0) != 0):
			print("basic solution:", bx, " is not feasible, please try again")
			continue
		print('Successfully converted the system into canonical form,')
		printls(nA, nb, nc, nz, yt, bx)
		print("good? (1: yes, 0: no)")
		if input() == '1':
			done  = True
	return nA, nb, nc, nz, bx, idxs_in

# get basis from input
def getBasis(A, dm, idxs):
	if (len(idxs) != dm[0]):
		return -2
	if (min(idxs) < 0 or max(idxs) >= dm[1]):
		return -1
	B = np.array([[A[j][i] for i in idxs] for j in range(dm[0])])
	if matrix_rank(B) != dm[0]:
		return -3
	return B

# updata basis indices
def updateBI(A, b, dm, k, B):
	rt = [-1 if A[i][k] <= 0 else b[i]/A[i][k] for i in range(dm[0])]
	print('ratio test: ',rt)
	curmin, cmp = 0, 0
	for i in range(0, dm[0]):
		if rt[i] != -1 and (curmin == 0 or curmin > rt[i]):
			curmin, cmp = rt[i], i
	print('index', B[cmp], 'is removed, index', k, 'is added to from a basis')
	B[cmp] = k
	# print('B:', B)
	return B


# get basic solution
def getBasicSolution(nb, idxs_in, dm):
	bx = [0 for i in range(dm[1])]
	for i in range(dm[0]):
		bx[idxs_in[i]] = nb[i]
	return np.array(bx)

# simplex
def simplex(dm, A, b, c, z):
	dm = checkIn(dm, A, b, c)
	if (type(dm) is int):
		return -1
	print('Input system')
	printls(A, b, c, z)
	A, b, c, z, x, B = canonical(A, b, c, z, dm)
	k = pstCpn(c.round(5))
	stepCount = 1
	while k != -1:
	# print([A[i][k] for i in range(dm[0])])
		if cmp0(np.array([A[i][k] for i in range(dm[0])]).round(5), lambda a: a > 0) == 0:
			print('unbounded solution')
			return 0
		print()
		print('step', str(stepCount) + ':')
		B = updateBI(A, b, dm, k, B)
		A, b, c, z, x, B = canonical(A, b, c, z, dm, B)
		k = pstCpn(c.round(5))
		stepCount += 1
		input()
	print()
	print("The solution")
	print('z =')
	print(NumToFrac(z))
	print('with')
	print('x =')
	printRow(x, 8)
	print("is optimal")


##### solve #####
simplex(*readIn())
