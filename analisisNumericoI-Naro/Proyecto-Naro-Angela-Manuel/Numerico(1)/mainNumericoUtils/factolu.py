#Usar para llamar a la funcion 
# factolu.sol(A,B)
#Donde A es la matriz
#B es el resultado
# -*- coding: utf-8 -*-
from numpy.linalg import inv
import numpy as np
def pivot(A,i,n):
    piv=np.argmax(abs(A[i:,i]))+i
    P=np.identity(n)
    #Intercambiamos filas.
    #A[[i,piv]]=A[[piv,i]]
    P[[i,piv]]=P[[piv,i]]
    #print(A)
    return P
def alpha(A,i,n):
    ap=np.identity(n)
    for j in range(i+1,n):
        ap[j][i]=-1*A[j][i]/A[i][i]   
    print(ap)
    return ap
#Para hallar el L debemos primero hacer n-1 iteraciones.
def factlu(A,n):
    #Ln-1Ln-2..L2L1
    Lx=np.identity(n)
    #Pn-1Pn-2....P2P1
    Px=np.identity(n)
    #Ln-iP-1..
    Ux=np.identity(n)
    for j in range(n-1):
        #Se va acumulando las matrices de permutacion
        print("Pivote")
        P=pivot(A,j,n)
        print(P)
        #Acumula.
        Px=np.dot(Px,P)
        #Halla Li
        Lp=alpha(np.dot(P,A),j,n)
        #Acumula.
        Lx=np.dot(Lx,Lp)
        #Acumula el U
        Ux=np.dot(Lp,np.dot(P,Ux))
        A=np.dot(Lp,np.dot(P,A))
        print(j)
        print(A)
    U=A
    print("--------------------------")
    print("Resultado final")
    print(U)
    P=Px
    L=np.dot(P,inv(Ux))
    return P,L,U
def sol(A,B):
    n=len(A)
    Ax=A
    #print(Ax)
    P, L, U = factlu(A,n)
    print(P)
    print(L)
    print(U)
    print("Solucion del problema")
    print("Matriz A")
    print(Ax)
    Pb=np.dot(P,B)
    y=np.linalg.solve(L,Pb)
    print(y)
    x=np.linalg.solve(U,y)
    print("Solución")
    print(x)
    return x.T