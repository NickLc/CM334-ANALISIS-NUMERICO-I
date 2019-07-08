#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

#--------------------------------------------------------------------------------------------

def pivoteo(a, i):
    """Calcula la matriz de pivotación para a en un índice dado. La pivotación
    se realiza por filas.

    Recibe una matriz a y el indice a pivotear i.
    Retorna la matriz de pivotación."""
    #La pivotacion consiste en cambiar de fila con aquella que tenga el maximo en la columna
    fil = a.shape[0]
    aux_i = i
    
    for k in range(i + 1, fil):
        if abs(a[k, i]) > abs(a[aux_i, i]) :
            aux_i = k
    P = np.identity(fil)
    if aux_i != i :
        P[[aux_i, i]] = P[[i, aux_i]]
    return P

#--------------------------------------------------------------------------------------------

def pivoteo_Total(A_b, v):
    
    #La pivotacion consiste en cambiar de fila con aquella que tenga el maximo 
    #en toda la sub matriz  a_kj  fila (i<k<n), col (i<j<n)
    (fil, col) = A_b.shape
    line = "-----------------------------------------"
    if v : print("{}\nSe pivotará antes de utilizar el método".format(line))

    for i in range(fil):
        P = pivoteo(A_b, i)
        A_b = P @ A_b
        if v : print("P_{}:\n{}\n".format(i+1,P))
        
    if v : print("Fin del pivoteo\n[A|b] despues del pivoteo total\n{}\n{}\n".format(A_b,line))
    
    return A_b

#--------------------------------------------------------------------------------------------

def get_L(A, i, v=False):
    """Esta funcion nos permite encontrar la matriz de transformacion de gauss, L_i
    Para una solución detallada, pasar v = True como parámetro.
    """
    (fil, col) = A.shape #Guarda el #filas y #columnas
    
    #Hallar alfa 
    alpha_i = np.zeros(fil) # alpha_i es declarado fila
    alpha_i[i + 1:] = A[i + 1: fil, i] / A[i, i]
    alpha_i = alpha_i.reshape(fil, 1) # transforma alpha_i en columna

    # Mostrar alfa
    if v : print("alpha_{}:\n{}".format(i+1, alpha_i))

    e_i = np.zeros(fil); e_i[i] = 1

    # Mostrar e(i)*alfa(i)
    if v : print("alpha_{} * e_{}:\n{}".format(i+1,i+1,alpha_i*e_i))

    L_i = np.identity(fil) - alpha_i * e_i

    # Mostrar L(i)
    if v : print("L_{}:\n{}".format(i+1,L_i))
    
    return L_i
#--------------------------------------------------------------------------------------------

def get_T(A, i, v=False):
    """Esta funcion nos permite encontrar la matriz de transformacion T_i
    Para una solución detallada, pasar v = True como parámetro.
    """
    (fil, col) = A.shape # Guarda el #filas y #columnas
    
    # Hallar alfa 
    alpha_i =  A[:, i] / A[i, i]# alpha_i es declarado  como vector fila
    alpha_i[i] = -alpha_i[i] / A[i, i] + 1
    alpha_i = alpha_i.reshape(fil, 1) # transforma alpha_i en columna

    # Mostrar alfa
    if v : print("alpha_{}:\n{}".format(i+1, alpha_i))

    e_i = np.zeros(fil); e_i[i] = 1

    # Mostrar e(i)*alfa(i)
    if v : print("alpha_{} * e_{}:\n{}".format(i+1,i+1,alpha_i*e_i))

    T_i = np.identity(fil) - alpha_i * e_i

    # Mostrar T(i)
    if v : print("T_{}:\n{}".format(i+1,T_i))
    
    return T_i
#--------------------------------------------------------------------------------------------

def get_L_and_P_LU(L_list, P_list):
    
    # Las matrices de permutacion son iguales a su inversa P = inv(P)   
    P_inv = P_list
    
    # La inversa de la matrice L_i es igual L_i pero cambiado el signo 
    # en toda la columan i, menos en la diagonal 
    
    fil,col = L_list[0].shape
    L = np.identity(fil)
    P = np.identity(fil) #Producto de las matrices de permutacion
    
    # Hallar la inversa de L_i 
    for i in range(len(L_list)):
        L = L_list[i]
        L[:,i:i+1] *=-1
        L[i,i] *=-1
        L_list[i] = L
    
    # L = P_1^-1 * L_1^-1 * ....... * P_n^-1 * L_n^-1
    for i in reversed(range(len(L_list))):
        L = L_list[i] @ L
        L = P_list[i] @ L
    
    # P = P_n * P_(n-1) * ........ * P_1  
    for i in range(len(P_list)):
        P = P_list[i] @ P   

    L = P @ L
       
    return L,P

#--------------------------------------------------------------------------------------------

def resolverMTriangularSup(U, b):
    """recibe un array U triangular superior, b vector de coeficientes
    retorna x vector solución"""
    (fil, col) = U.shape
    x = np.zeros(col)
    for i in reversed(range(fil)):
        suma = U[i, i + 1:].dot(x[i + 1:]) # sumatoria
        x[i] = (b[i] - suma) / U[i, i]
    return x

#--------------------------------------------------------------------------------------------

def resolverMTriangularInf(L, b):
    """recibe un array L triangular inferior, b vector de
    coeficientes.
    retorna x vector solución"""
    (fil, col) = L.shape
    x = np.zeros(col)
    for i in range(fil):
        suma = L[i, :i + 1].dot(x[:i + 1])
        x[i] = (b[i] - suma) / L[i, i]
    return x

#--------------------------------------------------------------------------------------------

def resolverMDiagonal(D, b):
    """recibe un array D diagonal, b vector decoeficientes.
    retorna x vector solución"""
    (fil, col) = D.shape
    x = np.zeros(col)
    for i in range(fil):
        x[i] = b[i] / D[i, i]
    return x
    
#--------------------------------------------------------------------------------------------

def esSimetrica(A):
    fil, col = A.shape
    assert fil == col, "La matriz de coeficientes no es cuadrada."
    bol = not(False in (A == A.T))
    return bol

#--------------------------------------------------------------------------------------------

if __name__ == 'toolNick' :
    print("toolNick se ha importado correctamente.")

