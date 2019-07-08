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

def pivoteo_parlett(a, i):
    """Calcula la matriz de pivotación para a en un índice dado. La pivotación
    se realiza por filas.

    Recibe una matriz a y el indice a pivotear i.
    Retorna la matriz de pivotación."""
    #La pivotacion consiste en cambiar de fila con aquella que tenga el maximo en la columna
    fil = a.shape[0]
    aux_i = i + 1
    
    for k in range(i + 2, fil):
        if abs(a[k, i]) > abs(a[aux_i, i]) :
            aux_i = k
    P = np.identity(fil)
    if aux_i != i + 1:
        P[[aux_i, i+1]] = P[[i+1, aux_i]]
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

# decompose
    # A: Matrix

    # return:  D = Diagonal matrix from A,
    #          E: negative of l.tri. matrix from A (no diag),
    #          F: negative of u.tri. matrix from A (no diag)
def decompose(A, verbose=False):
    E = -np.tril(A, -1)
    F = -np.triu(A, 1)
    D = np.diag(A.diagonal())
    if verbose:
        print("D: \n", D)
        print("E: \n", E)
        print("F: \n", F)

    return D,E,F

#--------------------------------------------------------------------------------------------

def Inverse(A, v = False, pti = False, pp = False):
    """A matriz de coeficientes 
    v nos muestra el procedimiento detalldo de la eliminacion
    pti se utiliza si requiere pivotacion total inicio
    pp siempre se realiza el pivoteo parcial
    inv retorna la inversa de la matriz A"""
     
    (fil, col) = A.shape #Guarda el #filas y #columnas
    line = "=============================================" 
    if v : print("========== Hallamdo la inversa: =============")
    # Creamos una lista para almacenar las matrices T
    T_list = [] 
    if pti: A = pivoteo_Total(A, v)

    for i in range(fil): # se condiera dim(a)
        if A[i, i] == 0 :
            P = pivoteo(A, i)   
            if v : print("P_\n{}:".format(P))
            A = P @ A
    
        if pp:
            P = pivoteo(A, i)   
            if v : print("P_\n{}:".format(P))
            A = P @ A
        
        # Obtenemos el T_i
        T_i = get_T(A, i, v)
        T_list.append(T_i) # Agregamos para hallar la inversa
        
        A = T_i @ A
        
        # Mostrar T(i) * (A|b)
        if v : print("T_{} * [A|b]:\n{}\n{}".format(i+1,A,line))
        
    
    A_inv = np.identity(fil)
    for i in range(fil):
        A_inv = T_list[i] @ A_inv 
    if v: print("La matriz inversa es:\n{}\n{}".format(A_inv,line))
    return A_inv

#--------------------------------------------------------------------------------------------

def get_L(A, i, v = False,r = False, vr = 0):
    """Esta funcion nos permite encontrar la matriz de transformacion de gauss, L_i
    Para una solución detallada, pasar v = True como parámetro.
    Si se quiere hacer un redondeo r=Valor de redondeo
    """
    (fil, col) = A.shape #Guarda el #filas y #columnas
    
    #Hallar alfa 
    alpha_i = np.zeros(fil) # alpha_i es declarado fila
    alpha_i[i + 1:] = A[i + 1: fil, i] / A[i, i]
    alpha_i = alpha_i.reshape(fil, 1) # transforma alpha_i en columna
    if r : alpha_i = redondeo(alpha_i, vr)
    
    # Mostrar alfa
    if v : print("alpha_{}:\n{}".format(i+1, alpha_i))

    e_i = np.zeros(fil); e_i[i] = 1

    # Mostrar e(i)*alfa(i)
    if v : print("alpha_{} * e_{}:\n{}".format(i+1,i+1,alpha_i*e_i))

    L_i = np.identity(fil) - alpha_i * e_i
    if r : L_i = redondeo(L_i, vr)
    # Mostrar L(i)
    if v : print("L_{}:\n{}".format(i+1,L_i))
    
    return L_i

#--------------------------------------------------------------------------------------------

def get_L_parlet(A, i, v = False,r = False, vr = 0):
    """Esta funcion nos permite encontrar la matriz de transformacion de gauss, L_i
    Para una solución detallada, pasar v = True como parámetro.
    Si se quiere hacer un redondeo r=Valor de redondeo
    """
    (fil, col) = A.shape #Guarda el #filas y #columnas
    
    #Hallar alfa 
    alpha_i = np.zeros(fil) # alpha_i es declarado fila
    alpha_i[i + 2:] = A[i + 2: fil, i] / A[i+1, i]
    alpha_i = alpha_i.reshape(fil, 1) # transforma alpha_i en columna
    if r : alpha_i = redondeo(alpha_i, vr)
    
    # Mostrar alfa
    if v : print("alpha_{}:\n{}".format(i+1, alpha_i))

    e_i = np.zeros(fil); e_i[i+1] = 1

    # Mostrar e(i+1)*alfa(i)
    if v : print("alpha_{} * e_{}:\n{}".format(i+1,i+2,alpha_i*e_i))

    L_i = np.identity(fil) - alpha_i * e_i
    if r : L_i = redondeo(L_i, vr)
    # Mostrar L(i)
    if v : print("L_{}:\n{}".format(i+1,L_i))
    
    return L_i

#--------------------------------------------------------------------------------------------

def get_T(A, i, v = False):
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
    r si se quiere un redondeo de vr cifras
    retorna x vector solución
    """
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

def expandirQ(A, n):
    (fil, col) = A.shape
    u = np.zeros((n, col))
    l = np.vstack((np.identity(n), np.zeros((fil, n))))
    A = np.r_[u, A]
    A = np.c_[l, A]
    return A

#--------------------------------------------------------------------------------------------

def back_solve_triangular(A, b, v=False):
    x = np.empty(A.shape[1])
    rows = A.shape[0]
    line = '==========================================='
    for i in range(rows-1, -1, -1):
        pre_sum = A[i, i+1:].dot(x[i+1:])
        x[i] = (b[i] - pre_sum)/A[i,i]
        if v : print("x_{}\n{}\n{}".format(i+1, x[i],line))
    if v : print("x: {}\n{}\n".format(x,line))
    return x
    
#--------------------------------------------------------------------------------------------

def esSimetrica(A):
    fil, col = A.shape
    assert fil == col, "La matriz de coeficientes no es cuadrada."
    bol = not(False in (A == A.T))
    return bol

#--------------------------------------------------------------------------------------------

def max_norm(A):    
    sum_row = 0
    nrow = np.shape(A)[0]
    for i in range(nrow):
        s = sum(np.abs(A[i,:]))
        if s > sum_row: sum_row = s
    return sum_row

#--------------------------------------------------------------

# calcula la condicional de la matriz A con la norma del maximo
def cond(A, v = False):    
    invA = np.linalg.inv(A) 
    condicion = max_norm(A)*max_norm(invA)
    if v: 
        if condicion < 1: print("A esta BIEN CONDICIONADA")
        if condicion < 100 : print("Es posible que A tenga un BUEN COMPORTAMIENTO") 
        else :print("Es posible que A tenga un MAL COMPORTAMIENTO")
    return condicion

#--------------------------------------------------------------

# Error relativo devido a la aplicacion del metodo "msolve"
def error_relativo(A, b, x0):
    """El error realativo es: 
    ||R||/(cond(A)*||b|| < ||E||/||x|| < cond(A)*||R||/||b||"""
    
    cota = []
    fil, col = A.shape
    b = b.reshape(fil,1)
    x0 = x0.reshape(fil,1)
    r = np.dot(A, x0) - b   #vector residuo
    
    cotinf = max_norm(r)/(cond(A)*max_norm(b)) # cota inferior del error relativo
    cota.append(cotinf)
    
    cotsup = cond(A)*max_norm(r)/max_norm(b) # cota superior del error relativo
    cota.append(cotsup)
    
    return cota

#--------------------------------------------------------------
def redondeo(A, vr):
    "Redondea el los elemento de A, a vr cifras"
    fil,col = A.shape
    for f in range(fil):
        for c in range(col):
            A[f][c] = round(A[f][c], vr)
    return A

#--------------------------------------------------------------
def resultado_real(A, b):
    "Calcula el resultado real de la ecuacion Ax=b"
    return np.linalg.solve(A,b)

#--------------------------------------------------------------
def mostrar_Sistema(A,b):
    print("\nSistema de ecuaciones:")
    fil,col = A.shape
    for i in range(fil):
        row = ["{0:3g}*x{1}".format(A[i, j], j + 1) for j in range(col)]
        print("[{0}] = [{1:3g}]".format(" + ".join(row), b[i]))

#--------------------------------------------------------------
def comprobar_Metodo(A,b,x):
    """Recibe la matriz A, b, x solucion del metodo
    Imprime el resultado del metodo, resultado real, cota de error, condicion """
    print("Resultado del metodo:\n{}\n".format(x))
    print("Resultado real:\n{}\n".format(resultado_real(A, b)))
    print("Cotas del error relativo:\n{}\n".format(error_relativo(A, b, x)))
    print("Cond(A) es: {}".format(cond(A, v=True)))
    
#--------------------------------------------------------------    
if __name__ == 'toolNick' :
    print("toolNick se ha importado correctamente.")

