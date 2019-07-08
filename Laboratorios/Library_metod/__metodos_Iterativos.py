#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from tools.toolNick import *
import numpy as np

def Jacobi(A, b, tau = 1e-8, n_iter=1000, v = False):
    """METODO DE RESOLUCION DE JACOBI
    Recibe A = Matriz, b = Array, tau = valor min de error, 
    n_iter = numero de iteraciones, v = verbosidad
    Se resolverá la ecuacion iterativa: X_(k+1) = J*X(k) + c 
    con G = I - inv(D)*A   &   c = inv(D)*b 
    Retorna la solucion"""

    x_0 = np.zeros_like(b)
    D, E, F = decompose(A, v)
    
    if v : print("========== Invirtiendo matriz D: =============")
    D_inv = Inverse(D, v=False)
    J = np.eye(D_inv.shape[0]) - np.matmul(D_inv, A)
    
    if v : print("J:\n{}".format(J))
    
    x = x_0

    c = np.matmul(D_inv, b)

    if v : print("c:\n{}".format(c))
    if v : print("\n========== Iterando la solucion =============\n")
    if n_iter >= 0:
        for i in range(0, n_iter):
            if v : print("Iteration {0}:x = {1}".format(i, x))
            
            if(np.linalg.norm(x - (np.matmul(J, x) + c)) < tau):
                if v : print("La solución converge con una variación de tau = {}\nSaliendo...".format(tau))
                break
            x = np.matmul(J, x) + c

    return x

def Gauss_Seidel(A, b, tau = 1e-8, n_iter=1000, v = False):
    """ METODO DE RESOLUCION DE GAUSS-SEIDEL
    Recibe A = Matriz, b = Array, tau = valor min de error, 
    n_iter = numero de iteraciones, v = verbosidad
    Se resolverá la ecuacion iterativa: X_(k+1) = G*X(k) + c 
    con G = I-(D-E)\'A y c = (D-E)\'b
    Retorna la solucion"""
        
    x_0 = np.zeros_like(b)
    D, E, F = decompose(A, v)
    
    if v : print("========== Invirtiendo matriz (D-E): =============")
    
    DE_inv = Inverse(D - E)
    G = np.matmul(DE_inv, F)
    
    if v:print("G:\n{}".format(G))
    
    x = x_0
    c = np.matmul(DE_inv, b)
    
    if v : print("c:\n{}".format(c))
    
    if n_iter >=0 :
        for i in range(0, n_iter):
            if v : print("Iteration {0}:x = {1}".format(i, x))
            lim = np.linalg.norm(x - (np.matmul(G, x) + c))
            if v : print("Lim = ",lim)
            if(lim < tau):
                if v : print("La solución converge con una variación de tau = {}\nSaliendo...".format(tau))
                break
            x = np.matmul(G, x) + c 

    return x

def SOR(A, b, w, tau = 1e-8, n_iter=1000, v=False):
    """ METODO DE RESOLUCION SOR
    Recibe A = Matriz, b = Array, w = parametro de relajacion 0<w<2, tau = valor min de error, 
    n_iter = numero de iteraciones, v = verbosidad 
    Se resolverá la ecuacion iterativa: X_(k+1) = S*X(k) + c 
    donde S = (D-wE)*[(1-w)D + wF]  y  c = (D-wE)*b 
    Retorna la solucion"""

    x_0 = np.zeros_like(b)
    D, E, F = decompose(A, v)
	
    D_wE = D - w*E

	#=========================================
    if v:
        print("(D - wE):\n{}".format(D_wE))		
        print("========== Invirtiendo matriz (D - wE): =============")
	#=========================================
    D_wE_inv = Inverse(D_wE, v)

    S = np.matmul( D_wE_inv, ( ((1-w)*D) + (w*F) ) )
	#=========================================
    if v:
        print("S:\n{}".format(S))
    #=========================================

    x = x_0

    b_hat = w*np.matmul(D_wE_inv, b)
    
    #=========================================
    if v:
        print("c:\n{}".format(b_hat))
    #=========================================

    if n_iter >= 0:
        for i in range(0, n_iter):
            if v:
                print("x_{}\t{}".format(i, x))
            if(np.linalg.norm(x - (np.matmul(S, x) + b_hat)) < tau):
                if v : print("La solución converge con una variación de tau = {}\nSaliendo...".format(tau))
                break
            x = np.matmul(S, x) + b_hat
            
    return x

def SSOR(A, b, w, tau = 1e-8, n_iter=1000, v=False):
    """ METODO DE RESOLUCION SSOR
    Recibe A = Matriz, b = Array, w = parametro de relajacion 0<w<2, tau = valor min de error, 
    n_iter = numero de iteraciones, v = verbosidad 
    Se resolverá la ecuacion iterativa: X_(k+1) = S*X(k) + c 
    donde SS = (D-wF)\'[(1-w)D + wE](D-wE)\'[(1-w)D + wF]  
    c = w(D-wF)\'{ [(1-w)D + wE](D-wE)\' + I }\'b
    Retorna la solucion"""
    
    x_0 = np.zeros_like(b)
    D, E, F = decompose(A, v)
    
    D_wE = D - w*E
    D_wF = D - w*F

    #=========================================
    if v:
        print("(D - wE):\n{}".format(D_wE))
        print("(D - wF):\n{}".format(D_wF))		
        print("========== Invirtiendo matriz (D - wE): =============")
    #=========================================
    D_wE_inv = Inverse(D_wE, v)

    #=========================================
    if v:
        print("========== Invirtiendo matriz (D - wF): =============")
    #=========================================
    
    D_wF_inv = Inverse(D_wF, v)

    SS_1 = np.matmul( D_wF_inv, ( ((1-w)*D) + (w*E) ) )
    SS_2 = np.matmul( D_wE_inv, ( ((1-w)*D) + (w*F) ) )
    SS = np.matmul(SS_1, SS_2)
    
    #=========================================
    if v:
        print("SS:\n{}".format(SS))
    #=========================================

    x = x_0

    b_hat = np.matmul( (1-w)*D + w*E , D_wE_inv) + np.eye(D.shape[0])
    b_hat = w*np.matmul(D_wF_inv, b_hat)
    b_hat = np.matmul(b_hat, b)
    
    #=========================================
    if v:
        print("c:\n{}".format(b_hat))
    #=========================================

    if n_iter >= 0:
        for i in range(0, n_iter):
            if v:
                print("x_{}\t{}".format(i, x))
            if(np.linalg.norm(x - (np.matmul(SS, x) + b_hat)) < tau):
                if v: print("La solución converge con una variación de tau = {}\nSaliendo...".format(tau))
                break
            x = np.matmul(SS, x) + b_hat
            
    return x

def Gradiente_Conjugado(A, b, tau=1e-8, n_iter=1000, delta=1e-2, v=False):
	""" METODO DE RESOLUCION POR GRADIENTE CONJUGADO
    Recibe A = Matriz, b = Array, tau = valor min de error, 
    n_iter = numero de iteraciones, delta = limite inferior de la velocidad, v = verbosidad 
    Se resolverá la ecuacion iterativa: X_(k+1) = X(k) + t_(k)*v_(k) 
    Retorna la solucion"""
	x_0 = np.zeros_like(b)
	r = b - np.matmul(A, x_0)
	vel = r[:]
	c = r.dot(r)
	x = x_0
	
	if v : print("{}\n{}\n{}".format(r,vel,c))

	for i in range(1, n_iter):
		if np.sqrt(vel.dot(vel)) < delta:
			break
		z = np.matmul(A, vel)
		t = c/vel.dot(z)
		x = x + t*vel
		r = r - t*z
		d = r.dot(r)
		if d < tau:
			if v : print("La solución converge con una variación de tau = {}\nSaliendo...".format(tau))
			break
		vel = r + (d/c)*vel
		c = d
		if v : print("x_{}\n{}".format(i, x))
	return x

#==============================================================================================================

if __name__ == '__metodos_Iterativos' :
    print("__metodos_Iterativos se ha importado correctamente.")

if __name__ == '__main__':
	A = np.array([[10., -1., 2., 0.],
              [-1., 11., -1., 3.],
              [2., -1., 10., -1.],
              [0., 3., -1., 8.]])
	# initialize the RHS vector
	b = np.array([6., 25., -11., 15.])
	mostrar_Sistema(A,b)
	print('\nMetodos Iterativos')
	print("Jacobi:\t\t{}".format(Jacobi(A,b)))
	print('Gauss_Seidel:\t{}'.format(Gauss_Seidel(A,b)))
	print('SOR:\t\t{}'.format(SOR(A,b,w=0.9)))
	print('SSOR:\t\t{}'.format(SSOR(A,b,w=0.9)))
	print('Grad Conjugado:\t{}'.format(Gradiente_Conjugado(A,b)))
