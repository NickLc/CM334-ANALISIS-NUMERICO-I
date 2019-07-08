#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
numericoUtils

Módulo para el curso Análisis Numérico
Se encuentran implementados métodos numéricos para resolver ecuaciones lineales.
@autor: gwynplaine
"""

import numpy


def impMatrix(a):
    """recibe una lista e imprime sus componentes"""

    for x in a:
        print(x)


def pivoteo(a, i):
    """recibe una matriz a y el indice a pivotear i devuelve la matriz de pivotación"""
        
    fil = a.shape[0]
    aux_i = i
    
    for k in range(i + 1, fil):
        if abs(a[k, i]) > abs(a[aux_i, i]) :
            aux_i = k
    P = numpy.identity(fil)
    if aux_i != i :
        P[[aux_i, i]] = P[[i, aux_i]]
    return P

def getR(A, b, x):
    """recibe una matriz de coeficientes a, la matriz de resultado b y x virtual solución.
    Retorna r vector error residual."""
    
    ax = np.matmul(A, x)
    r = ax - b
    return r

def Cond(a, b, x, x_i, v = False):
    """Recibe el vector de coeficientes a, el vector resultado b, el vector solución x y el vector de solución aproximada x_i.
    Para obtener detalle del proceso, pasar el v = True como parámetro.
    Retorna True si es que la solución es aceptable y False si la solución no es aceptable según la norma infinita."""
    
    norm_b = numpy.linalg.norm(b, numpy.inf)
    condA = numpy.linalg.cond(a, numpy.inf)
    condA_inv = 1 / condA
    E = x - x_i # vector error
    
    if v :
        print("Vector error: ")
        print(E)
    
    E_r = numpy.linalg.norm(E, numpy.inf) / numpy.linalg.norm(x, numpy.inf) # Error relativo
    r = getR(a, b, x)
    
    if v :
        print("Vector error residual: ")
        print(r)
    
    norm_r = numpy.linalg.norm(r, numpy.inf)
    rb = r / norm_b
    return condA_inv * rb <= E_r and E_r <= condA * rb


def resolverMTriangularSup(U, b):
    """recibe un array U triangular superior, b vector de coeficientes
    retorna x vector solución"""
    (fil, col) = U.shape
    x = numpy.zeros(col)
    for i in reversed(range(fil)):
        suma = U[i, i + 1:].dot(x[i + 1:]) # sumatoria
        x[i] = (b[i] - suma) / U[i, i]
    return x


def inversa(a, v = False):
    """recibe una matriz a
    retorna la matriz inversa de a si esta existe
    Para una solución detallada, pasar v = True como parámetro."""
    (fil, col) = a.shape
    assert fil == col, "La matriz no es cuadrada."
    assert numpy.linalg.det(a) != 0, "Matriz singular!"
    E = numpy.identity(fil)
    for i in range(fil): # se condiera dim(A)
        if a[i, i] == 0 :
            P = pivoteo(a, i)
            
            if v :
                print("P_%d:" %i)
                print(P)

            a = numpy.matmul(P, a)
        alpha_i =  a[:, i] / a[i, i]# alpha_i es declarado  como vector fila
        alpha_i[i] = -alpha_i[i] / a[i, i] + 1
        alpha_i = alpha_i.reshape(fil, 1) # transforma alpha_i en columna
        
        if v :
            print("alpha_%d: " %i)
            print(alpha_i)
        
        e_i = numpy.zeros(fil); e_i[i] = 1
        
        if v :
            print("e_%d: " %i)
            print(e_i)
            print("alpha_%d * e_%d: " %(i, i))
            print(alpha_i * e_i)
                
        
        E_i = numpy.identity(fil) - alpha_i * e_i
        
        if v :
            print("E_%d: " %i)
            print(E_i)
        E = numpy.matmul(E_i, E)
        a = numpy.matmul(E_i, a)

        if v :
            print("E_%d * a: " %i)
            print(a)

    if v :
        print("a inversa: \n", E)
    return E




def elimGauss(a, b, p = False, v = False):
    """recibe una matriz cuadrada A
        y la matriz b. Utiliza el método de eliminación de Gauss.
        Retorna x vector solución.
        Para realizar un pivoteo total, pasar p = 'total' como parámetro.
        Para una solución detallada, pasar v = True como parámetro."""

    a_b = numpy.c_[a, b] # Matriz aumentada
    (fil, col) = a_b.shape
    
    if v :
        print("Matriz aumentada: ")
        print(a_b)

    if p :
        if v :
            print("Se pivotará antes de utilizar el método")
        
        for i in range(fil):
            P = pivoteo(a_b, i)
            numpy.matmul(P, a_b)
            if v :
                print("P_%i:" %i)
                print(P)

        if v :
            print("Fin del pivoteo. Procediendo a calcular")
            print("Matriz a_b después del pivoteo inicial:")
            print(a_b)

    for i in range(fil): # se condiera dim(a)
        if a_b[i, i] == 0 :
            P = pivoteo(a_b, i)
            
            if v :
                print("P_%d:" %i)
                print(P)

            a_b = numpy.matmul(P, a_b)
        alpha_i = numpy.zeros(fil) # alpha_i es declarado fila
        alpha_i[i + 1:] = a_b[i + 1: fil, i] / a_b[i, i]
        alpha_i = alpha_i.reshape(fil, 1) # transforma alpha_i en columna
        
        if v :
            print("alpha_%d: " %i)
            print(alpha_i)
        
        e_i = numpy.zeros(fil); e_i[i] = 1
        
        if v :
            print("e_%d: " %i)
            print(e_i)
            print("alpha_%d * e_%d: " %(i, i))
            print(alpha_i * e_i)
                
        
        L_i = numpy.identity(fil) - alpha_i * e_i
        
        if v :
            print("L_%d: " %i)
            print(L_i)
        
        a_b = numpy.matmul(L_i, a_b)

        if v :
            print("L_%d * a_b: " %i)
            print(a_b)
            
    print("a_b final: \n", a_b)
    u = a_b[:, :col - 1]
    b_sol = a_b[:, col - 1]
    return resolverMTriangularSup(u, b_sol)

def gaussJordan(a, b, p = False, v = False):
    """recibe una matriz cuadrada A
    y la matriz b. Utiliza el método Gauss-Jordan.
    Retorna x vector de resultados.
    Para realizar un pivoteo total, pasar p = 'total' como parámetro.
    Para una solución detallada, pasar, v = True como parámetro."""
    
    a_b = numpy.c_[a, b] # Matriz aumentada
    (fil, col) = a_b.shape
    
    
    if v :
        print("Matriz aumentada: ")
        print(a_b)

    if p :
        if v :
            print("Se pivotará antes de aplicar el método.")
        
        for i in range(fil):
            P = pivoteo(a_b, i)

            if v :
                print("P_%i:" %i)
                print(P)

        if v :
            print("Fin del pivoteo.")
            print("Matriz a_b despupes del pivoteo inicial:")
            print(a_b)
            print()
            print("Procediendo a calcular")
    
    for i in range(fil): # se condiera dim(A)
        if a_b[i, i] == 0 :
            P = pivoteo(a_b, i)
            
            if v :
                print("P_%d:" %i)
                print(P)

            a_b = numpy.matmul(P, a_b)
        alpha_i =  a_b[:, i] / a_b[i, i]# alpha_i es declarado  como vector fila
        alpha_i[i] = -alpha_i[i] / a_b[i, i] + 1
        alpha_i = alpha_i.reshape(fil, 1) # transforma alpha_i en columna
        
        if v :
            print("alpha_%d: " %i)
            print(alpha_i)
        
        e_i = numpy.zeros(fil); e_i[i] = 1
        
        if v :
            print("e_%d: " %i)
            print(e_i)
            print("alpha_%d * e_%d: " %(i, i))
            print(alpha_i * e_i)
                
        
        T_i = numpy.identity(fil) - alpha_i * e_i
        
        if v :
            print("T_%d: " %i)
            print(T_i)
        
        a_b = numpy.matmul(T_i, a_b)

        if v :
            print("T_%d * a_b: " %i)
            print(a_b)

            
    print("a_b final: \n", a_b)
    x = a_b[:, col - 1]
    return x
    
if __name__ == '__main__' :
    print("Se ha ejecutado numericoUtils.py como archivo principal.")
    print("Esto... No se ha programado nada para una situación como esta...")
    print("Sugiero que ejecutes python y luego importes el modulo.")
    print("Ten un buen día!")

if __name__ == 'numericoUtils' :
    print("numericoUtils se ha importado correctamente.")

