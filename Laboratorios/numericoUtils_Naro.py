#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
numericoUtils

Módulo para el curso Análisis Numérico.
Se encuentran implementados métodos numéricos para resolver ecuaciones lineales,
ecuaciones no lineales, cálculo de la ecuación característica de una matriz
cuadrada dada, métodos iterativos.

El módulo está desarrollado con verbosidad para mostrar soluciones detalladas.
Sin embargo, para comprender íntegramente las respuestas que los algoritmos
proporcionan, se sugiere estudiar y comprender la teoría del curso de análisis
númerico.

@autor: gwynplaine
"""

import numpy


def impMatrix(a):
    """recibe una lista e imprime sus componentes"""

    for x in a:
        print(x)


def pivoteo(a, i):
    """Calcula la matriz de pivotación para a en un índice dado. La pivotación
    se realiza por filas.

    Recibe una matriz a y el indice a pivotear i.
    Retorna la matriz de pivotación."""
        
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
    """recibe una matriz de coeficientes a, la matriz de resultado b y x
    virtual-solución.
    Retorna r vector error residual."""
    
    ax = numpy.matmul(A, x)
    r = ax - b
    return r

def Cond(a, b, x, x_i, v = False):
    """Recibe el vector de coeficientes a, el vector resultado b, el vector 
    solución x y el vector de solución aproximada x_i.
    Para obtener detalle del proceso, pasar el v = True como parámetro.
    Retorna True si es que la solución es aceptable y False si la solución no 
    es aceptable según la norma infinita."""
    
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
    rb = norm_r / norm_b
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

def resolverMTriangularInf(L, b):
    """recibe un array L triangular inferior, b vector de
    coeficientes.
    retorna x vector solución"""
    (fil, col) = L.shape
    x = numpy.zeros(col)
    for i in range(fil):
        suma = L[i, :i + 1].dot(x[:i + 1])
        x[i] = (b[i] - suma) / L[i, i]
    return x

def inversa(a, v = False):
    """recibe una matriz a
    retorna la matriz inversa de a si esta existe
    Para una solución detallada, pasar v = True como parámetro."""
    (fil, col) = a.shape
    assert fil == col, "La matriz no es cuadrada."
    assert numpy.linalg.det(a) != 0, "Matriz singular! Determinante \
    igual a cero."
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
        Para realizar un pivoteo total, pasar p = True como parámetro.
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
    Para realizar un pivoteo total, pasar p = True como parámetro.
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

def L1U_fact(a, b, p = False, v = False):
    """recibe una matriz  cuadrada a de coeficientes y la matriz b
    retorna x vector resultado<aún no implementado>
    para realizar un pivoteo total pasar p = True como parametro
    para una solución detallada pasar v = True como parametro."""
    
    (fil, col) = a.shape
    assert fil == col, "La matriz no es cuadrada."
    a_b = numpy.c_[a, b] # Matriz aumentada
    
    if v :
        print("Matriz aumentada: ")
        print(a_b)

    if p :
        if v :
            print("Se pivotará antes de utilizar el método")
        P_t = numpy.identity(fil)

        for i in range(fil):
            P_i = pivoteo(a_b, i)

            if v :
                print("P_%i:" %i)
                print(P_i)

            P_t = numpy.matmul(P_i, P_t)

        a_b = numpy.matmul(P_t, a_b)

        if v :
            print("Fin del pivoteo. Matriz a_b después del pivoteo inicial:")
            print(a_b)
            print("Matriz P de pivotación total: ")
            print(P_t)
    L_inv = numpy.identity(fil)
    for i in range(fil):
        if a_b[i, i] == 0:
            P_i = pivoteo(a_b, i)
            
            if v :
                print("P_%d" %i)
                print(P_i)

            L_inv = numpy.matmul(P_i, L_inv)
        alpha_i = numpy.zeros(fil) # declarado como vector fila
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
        L_inv = numpy.matmul(L_i, L_inv)
    print("a_b final: \n", a_b)
    L = numpy.linalg.inv(L_inv)
    U = a_b[:, :col]
    if v :
        if p :
            print("P:")
            print(P_t)
            print()
        print("L: ")
        print(L)
        print()
        print("U: ")
        print(U)
    y = resolverMTriangularInf(L, b)
    x = resolverMTriangularSup(U, y)
    return x

def LU1_fact(a, b, p = False, v = False):
    """recibe una matriz  cuadrada a de coeficientes y la matriz b
    retorna x vector resultado
    para realizar un pivoteo total pasar p = 'total' como parametro
    para una solución detallada pasar v = True como parametro."""
    
    (fil, col) = a.shape
    assert fil == col, "La matriz de coeficienes no es cuadrada."
    l = numpy.zeros((fil, col))
    u = numpy.identity(fil)

    a_b = numpy.c_[a, b] # Matriz aumentada
    
    if v :
        print("Matriz aumentada: ")
        print(a_b)

    if p :
        print("Se pivotará antes de utilizar el método.")
        P_t = numpy.identity(fil)
        for i in range(fil):
            P_i = pivoteo(a_b, i)

            if v :
                print("P_%i:" %i)
                print(P_i)

            P_t = numpy.matmul(P_i, P_t)
        a_b = numpy.matmul(P_t, a_b)
        if v :
            print("Fin del pivoteo. Matriz a_b después del pivoteo inicial:")
            print(a_b)
            print("Matriz P de pivotación total: ")
            print(P_t)
    
    a_p = a_b[:, :col] # actualizamos a la matriz a_pivotada

    for k in range(fil):
        for i in range(k, fil):
            l[i, k] = a_p[i, k] - l[i, :k].dot(u[:k, k])
        for j in range(k + 1, fil):
            u[k, j] = (a_p[k, j] - l[k, :k].dot(u[:k, j])) / l[k, k]

    if v :
        print("L:")
        print(l)
        print()
        print("U: ")
        print(u)

    y = resolverMTriangularInf(l, b)
    x = resolverMTriangularSup(u, y)

    return x
    
def fragmentar(A, v = True):
    """Fragmenta la matriz A en su matriz diagonal, triangular inferior y 
    triangular superior.

    Recibe una matriz cuadrada A
    Retorna L, D, U matrices triangular inferior, diagonal y triangular
    superior respectivamente de manera que  L + D + U = A.
    Para una solución detallada, pasar v = True como parámetro."""
    (fil, col) = A.shape
    assert fil == col, "La matriz no es cuadrada."

    L = numpy.tril(A, -1)
    D = numpy.diag(numpy.diag(A)) # No preguntes... Es la forma que encontré.
    U = numpy.triu(A, 1)

    if v :
        print("L:")
        print(L)
        print("D:")
        print(D)
        print("U:")
        print(U)

    return L, D, U

def gram_schmidt(A, b, p = False, v = False):
    """Soluciona el sistema Ax = b utilizando el método de Gram-Schmidt.

    Recibe una matriz cuadrada A y la matriz b.
    Retorna x vector resultado.
    Para realizar un pivoteo total, pasar p = True como parámetro.
    
    Para una solución detallada pasar v = True como parámetro.
    """

    (fil, col) = A.shape
    tam = b.shape[0]
    assert fil == tam, "El sistema no es posible."
    a_b = numpy.c_[A, b] # Matriz aumentada.

    if v :
        print("A_b: {}".format(a_b))

    e = numpy.zeros((fil, col))
    u = numpy.zeros((col, col))

    for j in range(col):
        e[:, j] = A[:, j]
        for i in range(j - 1):
            u[i, j] = e[:, i].dot(A[:, j])
            e[:, j] = e[:, j] - u[i, j] * e[:, i]
        u[j, j] = numpy.linalg.norm(e[:, j])
        e[:, j] = e[:, j] / u[j, j]

    if v :
        print("U:")
        print(u)
        print("E:")
        print(e)

    E_b = numpy.matmul(numpy.linalg.inv(e), b)

    if v :
        print("E^(-1) * b:")
        print(E_b)

    x = resolverMTriangularSup(u, E_b)
    return x


def expandirQ(A, n):
    (fil, col) = A.shape
    u = numpy.zeros((n, col))
    l = numpy.vstack((numpy.identity(n), numpy.zeros((fil, n))))
    A = numpy.r_[u, A]
    A = numpy.c_[l, A]
    return A


def Householder(A, v = False):
    """Recibe una matriz A y la factoriza utilizando proyecciones de
    Householder.

    Retorna las matrices Q.T y R tal que QR = A."""
    n_row, n_col = A.shape
    A_k = numpy.copy(A)
    n_row_k, n_col_k = A_k.shape
    R = numpy.zeros((n_row, n_col))
    Q = numpy.identity(n_row)

    for k in range(n_col):

        if v:
            print("A_{}:".format(k))
            print(A)   

        current_col = numpy.asmatrix(A_k[:, 0]).T
        e_k = numpy.asmatrix(numpy.repeat(0, A_k.shape[0])).T
        e_k[0, 0] = 1

        w = numpy.sign(current_col[0,0]) * numpy.linalg.norm(current_col, 2) * e_k + current_col

        if numpy.linalg.norm(w, 2) != 0 :
            w = w / numpy.linalg.norm(w, 2)

        if v :
            print("w_{}\n{}".format(k, w))

        H_k = numpy.eye(A_k.shape[0]) - 2 * numpy.matmul(w, w.T)

        if v :
            print("H_{}\n{}".format(k, expandirQ(H_k, k)))

        Q = numpy.matmul(expandirQ(H_k, k), Q)

        HxA = numpy.matmul(H_k, A_k)
        R[k:n_row,k:n_col] = numpy.copy(HxA)

        if v :
            print("R_{}\n{}".format(k, R))

        A_k = numpy.zeros((n_row_k-1, n_col_k-1))
        A_k[:,:] = R[k+1:,k+1:]
        n_row_k -= 1
        n_col_k -= 1


    if v :
        print("Q final:\n{}".format(Q.T))
        print("R final:\n{}".format(R))

    return Q.T, R


def SOR(A, b, x_0, w, tau, v = False):
    """Utiliza el método de relajación SOR para resolver el sistema Ax = b.

    Recibe A matriz de coeficientes, b matriz de términos independientes,
    x_0 es el primer resultado aproximado, w el parámetro de relajación,
    tau es el error admitido.
    Retorna x vector solución y iter número de iteraciones.

    Para una solución detallada, pasar v = True como parámetro."""
    (fil, col) = A.shape
    assert fil == col, "La matriz no es cuadrada."
    tam = b.shape[0]
    assert fil == tam, "El sistema no es posible."
    A_b = numpy.c_[A, b]
    
    if v :
        print("A_b:")
        print(A_b)
        print()
        print("x_0:")
        print(x_0)
    
    x_prev = x_0
    num_iter = 1
    x = numpy.zeros(tam)

    maxIter = 100 # Como precaución, por si el resultado diverge.

    while numpy.linalg.norm(b - A@x) / numpy.linalg.norm(b - A@x_0) >= tau :
        for i in range(tam):
            sum1 = A[i, :i].dot(x[:i])
            sum2 = A[i, i + 1:].dot(x_prev[i + 1:])
            x[i] = w * (b[i] - sum1 - sum2) / A[i, i] + (1 - w) * x_prev[i]
        
        if v :
            print("x_%d:" %num_iter)
            print(x)

        if num_iter == maxIter :
            print("Se alcanzó el número máximo de iteraciones.")
            break

        x_prev = x
        num_iter = num_iter + 1

    print("Se concluyó después de %d iteraciones." %num_iter)
    return x, num_iter

def gauss_Seidel(A, b, x_0, tau, v = False):
    """Utiliza el método Gauss-Seidel para resolver el sistema Ax = b.

    Recibe A matriz de coeficientes, b matriz de términos independientes,
    x_0 es el primer resultado aproximado, tau es el error admitido.
    Retorna x vector solución y iter número de iteraciones.
    Para una solución detallada, pasar v = True como parámetro."""
    (fil, col) = A.shape
    assert fil == col, "La matriz no es cuadrada."
    tam = b.shape[0]
    assert fil == tam, "El sistema no es posible."
    A_b = numpy.c_[A, b]
    
    if v :
        print("A_b:")
        print(A_b)
        print()
        print("x_0:")
        print(x_0)
    
    x_prev = x_0
    num_iter = 1
    x = numpy.zeros(tam)

    maxIter = 100 # Como precaución, por si el resultado diverge.

    while numpy.linalg.norm(b - A@x) / numpy.linalg.norm(b - A@x_0) >= tau :
        for i in range(tam):
            sum1 = A[i, :i].dot(x[:i])
            sum2 = A[i, i + 1:].dot(x_prev[i + 1:])
            x[i] = (b[i] - sum1 - sum2) / A[i, i]
        
        if v :
            print("x_%d:" %num_iter)
            print(x)

        if num_iter == maxIter :
            print("Se alcanzó el número máximo de iteraciones.")
            break

        x_prev = x
        num_iter = num_iter + 1

    print("Se concluyó después de %d iteraciones." %num_iter)
    return x, num_iter

def gradConjugado(A, b, x_0, maxIter, tol, v = False):
    """Resuelve el sistema de ecuaciones Ax = b utilizando el método del
    gradiente conjugado.

    Recibe una matriz A de coeficientes, b matriz de términos independientes,
    x_0 solución aproximada inicial, maxIter número máximo de iteraciones,
    tol tolerancia.
    Retorna x vector solución, numIter número de iteraciones.
    Para una solución detallada, pasar v = True como parámetro."""

    (fil, col) = A.shape
    assert fil == col, "Matriz no cuadrada."
    assert numpy.all(A == A.T), "Matriz no simétrica."

    x = x_0
    r = b - A.dot(x_0) # Error
    vr = r
    c = r.dot(r)
    convergencia = False

    for numIter in range(maxIter):
        error = numpy.linalg.norm(vr)        
        
        if error < tol :
            convergencia = True
            break
        else:
            z = A.dot(vr)
            t = c / vr.dot(z)
            x = x + t * vr
            r = r - t * z
            d = r.dot(r)
            
            if v :
                print()
                print("x_{}: {}.".format(numIter, x))
                print("r_{}: {}.".format(numIter, r))
                print()

            if d ** 2 < tol :
                convergencia = True
                break
            else:            
                vr = r + (d / c) * vr
                c = d
    
    if not convergencia :
        print("No se alcanzó la convergencia.")
        return None, None
    else :
        return x, numIter

def Newton(x_0, f, df, max_iter = 100, E = 1E-3, v = False):
    """Encuentra la raíz de una función f utilizando el método de Newton.

    Reciba la primera aproximación x_0, la función f, la derivada df,
    maxIter número máximo de iteraciones, E error permitido.
    Retorna x solución y numIter número de iteraciones.
    Para una solución detallada, pasar v = True como parámetro"""
    x_old = numpy.copy(x_0)
    x_new = numpy.copy(x_0)

    if v :
        print("i\tx\t\t\t\tf(x)\t\t\t\tdf(x)")
    for numIter in range(0, max_iter):
        if v :
            print("{}\t{}\t\t\t\t{}\t\t\t\t{}".format(numIter, x_new, f(x_new), df(x_new)))

        x_new = x_old - f(x_old)/df(x_old)
        if numpy.abs(x_new-x_old)/numpy.abs(x_new) < E:

            if v :
                print("Converge en la iteración:{}".format(numIter))
                print("Solución aproximada:{}".format(x_new))

            return x_new, numIter
        
        x_old = numpy.copy(x_new)

    if v :
        print("Max. iteraciones alcanzado.")
        print("Solución aproximada:{}".format(x_new))

    return x_new, numIter

def metBisecc(f, a, b, tol = 1E-3, maxIter = 100, v = False):
    """Encuentra la raiz de una función f utilizando el método de la bisección.
    
    Recibe f función objetivo, a y b puntos iniciales, tol tolerancia, maxIter 
    número máximo de iteraciones.
    Retorna x solución y numIter número de iteraciones.
    
    Para una solución detallada, pasar v = True como parámetro."""
    fa = f(a)
    fb = f(b)

    if fa * fb < 0 : # solucion en <a,b>
        for numIter in range(maxIter):
            if v :
                print("Analizando intervalo <{};{}>".format(a, b))
            x_m = 0.5 * (a + b) # la mitad del intervalo
            error = (b - a) / 2 ** (i + 1)
            fx_m = f(x_m)

            if v :
                print("\tpunto medio del intervalo: {}".format(x_m))
                print("\tf evaluada en el punto medio: {}".format(fx_m))
                print("\terror: {}".format(error))

            if fx_m == 0 or error < tol :
                x = x_m
                break
            else:
                if fa * fx_m < 0 :
                    b = x_m
                    fb = fx_m
                else:
                    a = x_m
                    fa = fx_m

        if v :
            print("Número de iteraciones = {}".format(numIter))
            print("Solución aproximada = {}".format(x))

    else:
        print("La función no cambia de signo.")
        print("Ergo, no hay raiz en el intervalo <{}; {}>.".format(a, b))
    return x, numIter

def jacobi_nolineal(f, Jf, tau, x_0, v = False):
    """Encuentra la raiz de una función f utilizando el método no lineal de 
    jacobi.

    Recibe f función objetivo, Jf jacobiano de f, tau error,
    x_0 vector inicial.
    Retorna (x, numIter) vector solución y número de iteraciones necesarias
    para la convergencia.
    para una solución detallada pasar v = True como parámetro."""
    x_prev = x_0
    maxIter = 1000
    numIter = 1

    while numIter < maxIter:
        
        f_prev = f(x_prev)
        Jf_prev = Jf(x_prev)
        D_prev = numpy.diag(numpy.diag(Jf_prev))

        if v :
            print("f_{}:\n {}".format(numIter, f_prev))
            print("D:\n {}".format(D_prev))
        
        Japrox = D_prev

        if v :
            print("Japrox:\n{}".format(Japrox))
        
        Japrox = numpy.linalg.inv(Japrox)

        if v :
            print("Japroxinv:\n {}".format(Japrox))
        
        x = x_prev - Japrox.dot(f_prev)
        
        if v :
            print("x_{}: {}, {}".format(numIter, x[0], x[1]))

        if abs(numpy.linalg.norm(f(x))) <= tau :
            break

        numIter = numIter + 1
        x_prev = x

    return x, numIter


def SOR_noLineal(f, Jf, w, tau, x_0, v = False):
    """Encuentra la raiz de f función no lineal de varias variables utilizando
    el método de relajación SOR.

    Recibe f función objetivo, Jf jacobiano de f, w parámetro de relajación,
    tau error, x_0 vector inicial.
    Retorna (x, numIter) vector solución y número de iteraciones necesarias
    para la convergencia.   
    para una solución detallada pasar v = True como parámetro."""
    x_prev = x_0
    maxIter = 1000
    numIter = 1

    while numIter < maxIter:
        
        f_prev = f(x_prev)
        Jf_prev = Jf(x_prev)
        D_prev = numpy.diag(numpy.diag(Jf_prev))
        L_prev = numpy.tril(Jf_prev, -1)

        if v :
            print("f_{}:\n {}".format(numIter, f_prev))
            print("D:\n {}".format(D_prev))
            print("L:\n {}".format(L_prev))
        
        Japrox = ((1 - w) / w) * D_prev + L_prev

        if v :
            print("Japrox:\n{}".format(Japrox))
        
        Japrox = numpy.linalg.inv(Japrox)

        if v :
            print("Japroxinv:\n {}".format(Japrox))
        
        x = x_prev - Japrox.dot(f_prev)
        
        if v :
            print("x_{}, {}, {}".format(numIter, x[0], x[0]))

        if abs(numpy.linalg.norm(f(x))) <= tau :
            break

        numIter = numIter + 1
        x_prev = x

    return x, numIter

def krylov(A, y = None, v = True):
    """Obtiene los coeficientes de la ecuación característica de una matriz
    cuadrada A utilizando el método de Krylov.
    
    Recibe A matriz cuadrada, y vector arbitrario multiplicable con A.
    Retorna coeficientes de la ecuación característica de A.
    
    Para una solución detallada, pasar v = True como parámetro."""
    (fil, col) = A.shape
    
    assert fil == col, "Matriz no cuadrada."

    if y is None :
        y = numpy.zeros(fil)
        y[0] = 1

        if v :
            print("No se escogió y inicial. Se utilizará y = [1 0 0 ... 0].")

    else :
        dim_y, = y.shape # Solo tomamos el primer término del 'shape'
        assert fil == dim_y, "Dimensión del vector y no compatible con la matriz."
        if v :
            print("y: {}.".format(y))

    matCoef = y
    Ai_y = y # A ** 0 @ y

    for i in range(1, fil):
        Ai_y = A @ Ai_y

        if v :
            print("A^{} @ y : {}.".format(i, Ai_y))

        matCoef = numpy.c_[Ai_y, matCoef]
    
    b = -(A @ Ai_y)
    
    if v :
        print("Nuestra matriz de coeficientes quedó como:")
        print(matCoef)
        print("Nuestra matriz de términos independientes es:")
        print(b)
    
    # Se escoge eliminación de Gauss debido a sus pocos requerimientos
    # en relación a la matriz de coeficientes y a su velocidad de cómputo.
    return elimGauss(matCoef, b) # No es la única forma... Puede cambiarse.

def autovect_leverrierFaddev(A, v = False):
    """Obtiene los coeficientes de la ecuación característica de una matriz
    cuadrada A utilizando el método de Leverrier-Faddev.
    
    Recibe A matriz cuadrada, y vector arbitrario multiplicable con A.
    Retorna coeficientes de la ecuación característica de A.
    
    Para una solución detallada, pasar v = True como parámetro."""
    fil, col = A.shape
    assert fil == col, "La matriz no es cuadrada."

    coef = [] # Aquí se almacenan los coeficientes de la ecuación característica
    b_i = -(A.trace())
    coef.append(b_i)
    B_i = A # No se tratará como referencia

    if v :
        print("B_0:")
        print(B_i)
        print("b_0: {}.".format(b_i))

    for i in range(2, fil + 1):
        B_i = A @ (B_i + b_i * numpy.identity(fil))
        b_i = -1 * B_i.trace() / i

        if v :
            print("B_{}:".format(i))
            print(B_i)
            print("b_{}: {}.".format(i, b_i))
        
        coef.append(b_i)
    
    return coef

def autovect_Potencia(A, w = None, E = 1E-5, v = False):
    """Calcula el autovalor dominante de una matriz cuadrada A usando el
    método de potencia.
    
    Reciba matriz cuadrada A, w vector arbitrario multiplicable con A.
    Retorna autovalor dominante.

    Para una solución detallada, pasar v = True como parámetro.
    """
    fil, col = A.shape
    assert fil == col, "Matriz no cuadrada."

    if w is None :
        
        if v :
            print("Vector arbitrario no escogido.")
            print("Se usará w = [1 1 ...1]")
        
        w = numpy.ones(fil)
    else:
        dim_w, = w.shape # Solo tomamos el primer término del 'shape'
        assert dim_w == fil, "Dimensión del vector w no compatible con A."
    
    maxIter = 100
    y_i = A.dot(w)
    lambda_i = numpy.linalg.norm(y_i)
    w = (1 / lambda_i) * y_i

    if v :
            print("y_1: {}.".format(y_i))
            print("lambda_1: {}.".format(lambda_i))
            print("x_1: {}".format( w))

    for i in range(2, maxIter + 1):
        y_i = A @ w
        lambda_o = lambda_i
        lambda_i = numpy.linalg.norm(y_i)
        w = (1 / lambda_i) * y_i
        
        if v :
            print("y_{}: {}.".format(i, y_i))
            print("lambda_{}: {}.".format(i, lambda_i))
            print("x_{}: {}".format(i, w))
        
        if abs(lambda_i - lambda_o) < E :
            break 
    
    return lambda_i


def autovect_PotenciaInver(A, w = None, E = 1E-5, v = False):
    """Calcula el autovalor dominante de una matriz cuadrada A usando el
    método de potencia.
    
    Reciba matriz cuadrada A, w vector arbitrario multiplicable con A.
    Retorna autovalor dominante.

    Para una solución detallada, pasar v = True como parámetro.
    """
    pass

def interp_Lagrange(xs, ys , v = False):
    """Calcula la aproximación del polinomio interpolador en la forma de 
    Lagrange según los datos de entrada."""
    pass

def interp_Newton(xs, ys, v = False):
    pass


if __name__ == '__main__' :
    print("Se ha ejecutado numericoUtils.py como archivo principal.")
    print("Esto... No se ha programado nada para una situación como esta...")
    print("Sugiero que ejecutes python y luego importes el modulo.")
    print("Ten un buen día!")

    print("Problema 4:")
    A = numpy.array([[0, 0, 0.33], [0.18, 0, 0], [0, 0.71, 0.94]])
    l = krylov(A, y = numpy.array([1,0,0]),v = True)
    print("Coeficientes: {}".format(l))

    eig = autovect_Potencia(A, v = True)
    print("eigen: {}.".format(eig))

    print("vector propio asociado:")

    eig = autovect_Potencia(numpy.linalg.inv(A), v = True)
    print("eigen^-1: {}.".format(eig))
    print("eigen : {}".format(1 / eig))




    print("--------------------------------------------------------------")

if __name__ == 'numericoUtils' :
    print("numericoUtils se ha importado correctamente.")

