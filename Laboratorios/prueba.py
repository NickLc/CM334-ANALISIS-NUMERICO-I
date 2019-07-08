# -*- coding: utf-8 -*-

import numpy as np
import Library_metod.__metodos_Lineales as ml
import Library_metod.__metodos_Iterativos as mi

A = np.array([[10., -1., 2., 0.],
              [-1., 11., -1., 3.],
              [2., -1., 10., -1.],
              [0., 3., -1., 8.]])
# initialize the RHS vector
b = np.array([6., 25., -11., 15.])


#print(dir(mi))
#print(help(mi.Jacobi))

x1 = mi.Gradiente_Conjugado(A,b)
x2 = mi.Jacobi(A,b) 
x3 = mi.Gauss_Seidel(A,b)
print("x1 = {}\n x2 = {}\n x3 = {}".format(x1,x2,x3))
