import numpy as np
import matplotlib.pyplot as plt 
import sklearn as sk
import cvxpy as cp 

x = cp.Variable(5)
A = np.random.randn(5,5)
x0 = [0, 3, 0, 0.42, 0]
y = np.dot(A,x0) + np.random.randn(5)

objfun = cp.Minimize(cp.sum_squares(A@x0 - y) + 0.1 * cp.norm1(x))
problem = cp.Problem(objfun)

problem.solve(solver=cp.GUROBI,verbose=True)
print(x.value)


numbers = [0, 34, 14, 15, 16]

def fFlop(numberList):
    out = 0
    for k in numberList:
        out += k
    return out 

def fWhile(numberList):
    out, counter = 0, 0 
    while counter < len(numberList):
        out += numberList[counter]
        counter += 1
    return out

def fRec(numberList):
    if len(numberList)==0:
        return 0
    else:
        return numberList[0] + fRec(numberList[1:]) 


print(fFlop(numbers))
print(fWhile(numbers))
print(fRec(numbers))

l1, l2  =  ["a", "b", "c"],  [1, 2, 3]

def mergeAlter(list1, list2):
    out = [0]* (len(list1) + len(list2))
    out[::2] = list1
    out[1::2] = list2
    return out

print(mergeAlter(l1,l2))

ns = [50, 2, 1, 9]
from functools import cmp_to_key
def genBiggestNum(listDigits):
    out = sorted(listDigits, key = cmp_to_key(lambda k, m:-1 
                if str(m) + str(k) < str(k) + str(m) else 1)) 
    strOut = map(str,out)
    return int(''.join(strOut))

print(genBiggestNum(ns))



import itertools
digits = [1, 8, 3, 5, 6, 7, 2, 9]
operation = ["+", "-", ";"]
def give88(dig, oper):
    allPositibilties = list(itertools.combinations(oper,7))
    for k in allPositibilties:
        print(k)
    return allPositibilties

print(give88(digits,operation))