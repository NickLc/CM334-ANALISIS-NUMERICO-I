#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Library_metod.toolNick import *
import numpy as np


def elimGauss(a, b, pti = False, v = False, pp = False, r = False, vr = 0):
    """a matriz de coeficientes, b matriz independiente 
    v nos muestra el procedimiento detAllAdo de la eliminacion
    pti se utiliza si requiere pivotacion total inicio
    pp siempre se realiza el pivoteo parcial
    r si se quiere realizar un redondeo de vr cifras"""
    
    A_b = np.c_[a, b] # Matriz aumentada
    (fil, col) = A_b.shape #Guarda el #filas y #columnas
    line = "=============================================" 
    if r : A_b = redondeo(A_b, vr)
    if v : print("Matriz aumentada [A|b] al inicio:{} \n{}".format(A_b,line))
    if pti: A_b = pivoteo_Total(A_b, v)

    for i in range(fil): # se condiera dim(a)
        #El pivoteo parcio no el solo cuando se tiene 0 en la diagonal
        #Para realizar el pivoteo parcial eliminar la condicion A_b[i,i]==0
        
        #Si se quiere un pivoteo parcial quitar esta condicion
        if A_b[i, i] == 0 :
            P = pivoteo(A_b, i)
            if v : print("P_{}:\n{}".format(i+1,P))
            A_b = np.matmul(P, A_b)
            if r : A_b = redondeo(A_b, vr)
        
        if pp:
            P = pivoteo(A_b, i)
            if v : print("P_{}:\n{}".format(i+1,P))
            A_b = np.matmul(P, A_b)
            if r : A_b = redondeo(A_b, vr)
        #Obtenemos el L_i
        L_i = get_L(A_b, i, v) #Ya lo obtiene redondeado
        A_b = np.matmul(L_i, A_b)
        if r : A_b = redondeo(A_b, vr) 
        
        #mostrar L(i) * (A|b)
        if v : print("L_{} * A:\n{}\n{}".format(i+1,A_b,line))
            
    if v : print("Matriz aumentada [A|b] al final: \n", A_b)
    #U sera todo menos la ultima columna
    U = A_b[:, :col - 1]
    #b_sol sera solo la ultima columna
    b_sol = A_b[:, col - 1] 
    x = resolverMTriangularSup(U, b_sol)
    if r :
        x = x.reshape(1,col)
        print("Redondeando la solución")
        x = redondeo(x, vr)
    return x


def gauss_Jordan(a, b, v = False, pti = False, pp = False, inv = False):
    """a matriz de coeficientes, b matriz independiente 
    v nos muestra el procedimiento detAllAdo de la eliminacion
    pti se utiliza si requiere pivotacion total inicio
    pp siempre se realiza el pivoteo parcial
    inv retorna la inversa de la matriz A"""
     
    A_b = np.c_[a, b] # Matriz aumentada
    (fil, col) = A_b.shape #Guarda el #filas y #columnas
    line = "=============================================" 
    
    # Creamos una lista para almacenar las matrices T
    T_list = []
    if v: print("Matriz aumentada [A|b] al inicio:\n{}\n{}".format(A_b,line))
    
    if pti: A_b = pivoteo_Total(A_b, v)

    for i in range(fil): # se condiera dim(a)
        if A_b[i, i] == 0 :
            P = pivoteo(A_b, i)   
            if v : print("P_\n{}:".format(P))
            A_b = P @ A_b
            
        if pp:
            P = pivoteo(A_b, i)   
            if v : print("P_\n{}:".format(P))
            A_b = P @ A_b
        
        # Obtenemos el T_i
        T_i = get_T(A_b, i, v)
        T_list.append(T_i) # Agregamos para hallar la inversa
        
        A_b = T_i @ A_b
        
        # Mostrar T(i) * (A|b)
        if v : print("T_{} * [A|b]:\n{}\n{}".format(i+1,A_b,line))
        
    if inv:
        A_inv = np.identity(fil)
        for i in range(fil):
            A_inv = T_list[i] @ A_inv 
        if v: print("La matriz inversa de A es:\n{}\n{}".format(A_inv,line))
    
    if v: print("Matriz aumentada [A|b] al final: \n", A_b)
    x = A_b[:, col - 1]
    return x


def elim_LU(a, b, v = False, pti = False, pp = False):
    """a matriz de coeficientes, b matriz independiente 
    v nos muestra el procedimiento detAllAdo de la eliminacion
    pti se utiliza si requiere pivotacion total inicio
    pp siempre se realiza el pivoteo parcial"""

    A_b = np.c_[a, b] # Matriz aumentada
    (fil, col) = A_b.shape #Guarda el #filas y #columnas
    line = "=============================================" 
    if v: print("Matriz aumentada [A|b] al inicio:\n{}\n{}".format(A_b,line))
    
    if pti: A_b = pivoteo_Total(A_b, v)

    A = A_b[:, :col - 1] #A sera todo menos la ultima columna
    b = A_b[:, col - 1] #b sera solo la ultima columna
    
    #Creamos una lista para almacenar las matrices L y P
    L_list = []
    P_list = []  #Matrices de permutacion
    
    for i in range(fil): # se condiera dim(a)
        if A[i, i] == 0 :
            P = pivoteo(A, i)
            if v : print("P_{}:\n{}".format(i+1,P))
            A = P @ A
            P_list.append(P)
        if pp:
            P = pivoteo(A, i)
            if v : print("P_{}:\n{}".format(i+1,P))
            A = P @ A
            P_list.append(P)
            
        else: #En caso de no permutar agragamos la identidad
            P_list.append(np.identity(fil))  
            
        #Obtenemos el L_i
        L_i = get_L(A, i, v)
        L_list.append(L_i)
        
        A = L_i @ A
        
        #mostrar L(i) * (A)
        if v : print("L_{} * A:\n{}\n{}".format(i+1,A,line))
            
            
    U = A
    L,P = get_L_and_P_LU(L_list, P_list)
    if v : print("L_:\n{}\nU_:\n{}\nP_:\n{}".format(L,U,P))
    
    Pb = P @ b 
    #Tenemos L, U, P, b
    
    # Ax = b  &  PA = LU
    # PAx = Pb -> LUx = Pb 
    # Ly = Pb  hallamos y  ->  Ux = y  hallamos x

    y = resolverMTriangularInf(L, Pb)
    x = resolverMTriangularSup(U, y)
    return x


def fac_Grout_LU1(A, b, v = False, pti = False):
    """a matriz de coeficientes, b matriz independiente 
    v nos muestra el procedimiento detAllAdo de la eliminacion
    pti se utiliza si requiere pivotacion total inicio"""
    
    A_b = np.c_[A, b] # Matriz aumentada
    fil, col = A_b.shape #Guarda el #filas y #columnas
    #assert fil == col, "La matriz de coeficienes no es cuadrada."
    line = "==================================="    
    if v: print("Matriz aumentada [A|b] al inicio: \n", A_b)
    
    if pti:
        A_b = pivoteo_Total(A_b, v)
    
    # Separamos A y b de [A|b]
    A = A_b[:, :col - 1] #A sera todo menos la ultima columna
    b = A_b[:, col - 1] #b sera solo la ultima columna
    
    #Creamos la matrices L Y U1_ identidad
    L = np.zeros((fil, col))
    U = np.identity(fil)

    for k in range(fil):
        # Hallando L
        for i in range(k, fil):
            L[i, k] = A[i, k] - L[i, :k]@U[:k, k]
        
        #Hallando U
        for j in range(k + 1, fil):
            U[k, j] = (A[k, j] - L[k, :k]@U[:k, j]) / L[k, k]

    if v : print("{}\nL: \n{}\n{}\nU:\n{}\n".format(line,L,line,U))
        
    
    y = resolverMTriangularInf(L, b)
    x = resolverMTriangularSup(U, y)
    
    return x
    

def fac_Grout_L1U(A, b, v = False, pti = False):
    """a matriz de coeficientes, b matriz independiente 
    v nos muestra el procedimiento detAllAdo de la eliminacion
    pti se utiliza si requiere pivotacion total inicio"""
    
    A_b = np.c_[A, b] # Matriz aumentada
    fil, col = A_b.shape #Guarda el #filas y #columnas
    #assert fil == col, "La matriz de coeficienes no es cuadrada."
    line = "==================================="    
    if v: print("Matriz aumentada [A|b] al inicio: \n", A_b)
    
    if pti: A_b = pivoteo_Total(A_b, v)
    
    # Separamos A y b de [A|b]
    A = A_b[:, :col - 1] #A sera todo menos la ultima columna
    b = A_b[:, col - 1] #b sera solo la ultima columna
    
    #Creamos la matrices L Y U1_ identidad
    U = np.zeros((fil, col))
    L = np.identity(fil)

    for k in range(fil):
        # Hallando U
        for j in range(k, fil):
            U[k, j] = A[k,j] - L[k, :k]@U[:k, j]
                       
        #Hallando L
        for i in range(k + 1, fil):
             L[i, k] = (A[i, k] - L[i, :k]@U[:k, k]) /U[k, k]

    if v : print("{}\nL: \n{}\n{}\nU:\n{}\n".format(line,L,line,U))
    
    y = resolverMTriangularInf(L, b)
    x = resolverMTriangularSup(U, y)
    x = x[:-1]
    return x



def LDLt(A, b, v = False, pti = False, pp = False):
    """A matriz de coeficientes, b matriz independiente 
    v nos muestra el procedimiento detAllAdo de la eliminacion
    pti se utiliza si requiere pivotacion total inicio"""
    
    (fil, col) = A.shape
    assert fil == col, "La matriz de coeficientes no es cuadrada."
    assert np.all(np.linalg.eigvals(A) > 0), "Matriz de coeficientes no definida \
    positiva"
    #assert esSimetrica(A) == True, "La matriz no es simetrica, los resultados no seran buenos"

    #Si no es simetrica la multiplicamos por su transpuesta
    if not(esSimetrica(A)) : 
        if v : print("La matriz no es simetrica\n",A)
        A = A.T @ A;   b = A.T @ b
        if v : print("Convirtiendo a simetrica A = A.t*A\n",A)
            
    D = np.zeros((fil, col))
    L = np. identity(fil)
    P = np.identity(fil)   #Matriz de pivotacion
    line = "=================================="
    
    if v : print("Matriz aumentada [A|b] al inicio:\n{}\n{}".format(np.c_[A,b],line))
    
    ##==========================================================================
    # Pivoteamos a inicio
    if pti:
        P = np.identity(fil)
        for i in range(fil):
            P_i = pivoteo(A,i)
            P = P @ P_i
        
        A = P @ A ; A = A @ P.T
        b = P @ b ; b = b @ P.T 
        
    ##==========================================================================
    
    for k in range(fil):
        pre_sum = 0
        for p in range(k):
            pre_sum += D[p,p] * L[k, p]*L[k,p]
        D[k,k] = A[k,k] - pre_sum
        
        if D[k,k] == 0:
            break
            
        for i in range(k + 1, fil):
            pre_sum_2 = 0
            for p in range(k):
                pre_sum_2 += L[i, p]*L[k, p]*D[p,p]
           
            L[i,k] = (A[i,k] - pre_sum_2)/D[k,k]
        
        if v:
            print("D_{}.\n{}\n".format(k+1,D))
            print("L_{}.\n{}\n{}".format(k+1,L,line))
          
        
    # Se tiene L y D
    if v : print("D_\n{}\nL_\n{}".format(D,L))
    # Ax = b  &  A = LD(L_t)
    
    # Ax = b -> LD(L_t)x = b 
    # Ly = b  hallamos y  ->  D(L_t)x = y   
    # Dw = y  hallamos w  
    #(L_t)x = w  hallamos x
    
    y = resolverMTriangularInf(L, b)
    w = resolverMDiagonal(D, y)
    x = resolverMTriangularSup(L.T, w)
    
    return x


def cholesky(A, b, v = False, pti = False):
    """A matriz de coeficientes, b matriz independiente 
    v nos muestra el procedimiento detAllAdo de la eliminacion
    pti se utiliza si requiere pivotacion total inicio"""
    
    fil, col = A.shape
    assert fil == col, "La matriz de coeficientes no es cuadrada."
    assert np.all(np.linalg.eigvals(A) > 0), "Matriz de coeficientes no definida \
    positiva"
    ##assert esSimetrica(A) == True, "La matriz no es simetrica, los resultados no seran buenos"
    
    #Si no es simetrica la multiplicamos por su transpuesta
    if not(esSimetrica(A)) : 
        if v : print("La matriz no es simetrica\n",A)
        A = A.T @ A;   b = A.T @ b
        if v : print("Convirtiendo a simetrica A = A.t*A\n",A)
         
    line = "================================="
    if v : print("{}\nLa matriz aumentada es:\n{}".format(line,np.c_[A,b]))

    G = np.zeros((fil, col))
    
    ##==========================================================================
    # Pivoteamos a inicio
    if pti:
        P = np.identity(fil)
        for i in range(fil):
            P_i = pivoteo(A,i)
            P = P @ P_i
        
        A = P @ A ; A = A @ P.T
        b = P @ b ; b = b @ P.T 
        
    ##==========================================================================
    
    for k in range(fil):
        if A[k, k] > 0 :
            G[k, k] = np.sqrt(A[k, k] - G[k, :k].dot(G[k, :k]))
        for  i in range(k + 1, fil):
            G[i, k] = (A[i, k] - G[i, :k].dot(G[k, :k])) / G[k, k]
        if v : print("{}\nG_{}\n{}\n".format(line,k+1,G))
   
    # Ax = b  &  A = G * G_t
    if v : print("{}\nG_\n{}\n".format(line,G))
    # Ax = b -> G*G_t x = b 
    # Gy = b  hallamos y  ->  G_t * x = y   
    # G_t x = y  hallamos x  
    
    y = resolverMTriangularInf(G, b)
    x = resolverMTriangularSup(G.T, y)

    return x


def parlett_Reid(A, b, pp = False, v = False):
    """Recibe matriz A de coeficientes, b matriz de términos independientes.
    Resuelve el sistema utilizando el método de Parllet Reid, pp para pivotear parcial
    Retorna x vector solución.
    Para una solución detallada, pasar v = True como parámetro."""
    (fil, col) = A.shape
    assert fil == col, "La matriz no es cuadrada."
    assert np.all(A == A.T), "La matriz no es simétrica."

    P = np.identity(fil) # Matriz de pivotación
    L = np.identity(fil) # Matriz triangular superior
    T = np.empty(fil)
    line = "==============================================" 
    line2 = "------------------------------------------------"
    A_b = np.c_[A, b] # Matriz aumentada
    if v: print("Matriz aumentada [A|b] al inicio:\n{}\n{}".format(A_b,line))
    
    for i in range(fil - 2):
        P_i = np.identity(fil)
        
        if A[i, i] == 0 or pp : P_i = pivoteo_parlett(A, i)  
        if v : print("P_{} = \n{}".format(i+1, P_i))
        # P = P(n-2)*.......*P(1)
        P = P_i @ P    ;   L = P_i @ L 
        A = P_i @ A    ;   A = A @ P_i.T
        if v: print("P_{} * A * P_t_{} =\n{}\n{}".format(i+1, i+1, A,line2)) 
        
        #Obtenemos el L_i
        L_i = get_L_parlet(A, i, v=True)
        A = L_i @ A    ;   A = A @ L_i.T
        
        #mostrar L(i) * (A) * L_t(i)
        if v : print("L_{} * A * L{}_t:\n{}\n{}".format(i+1,i+1,A,line))
      
        L = L_i @ L
    
    T = A
    if v : print("T:\n{}\n".format(T))
    # L = (L * P_t)^-1
    L = np.linalg.inv(L @ P.T)
    
    #Tenemos L, T, P, b
    
    # Ax = b  &  P*A*P_t = L*T*L_t  => L*T*L_t*Px = Pb 
    # Lw = Pb -> w = T*L_t*Px  
    # T*y = w  -> y = L_t*Px
    # L_t*z = y -> z = Px
    # hallamos x
    
    w = resolverMTriangularInf(L, P @ b)
    if v : print("Lw = Pb\nw =\n{}.".format(w))
        
    y = gauss_Jordan(T, w)
    if v : print("Ty = w;\ny = {}.".format(w))
    
    z = resolverMTriangularSup(L.T, y)
    if v : print("L_tz = y;\nz = \n{}".format(y))
    
    x = P.T @ z
    if v : print("z = Px;\nx = \n{}.".format(x))

    return x



if __name__ == '__metodos_Lineales' :
    print("__metodos_Lineales se ha importado correctamente.")