# -*- coding: utf-8 -*-

import numpy


"""
   imprime la matriz dada como una
   lista de numeros
"""
def impMatrix(a):
	for x in a:
		print(x)


"""
   recibe una matriz a y el indice a pivotear i
   devuelve la matriz de pivotación
"""

def pivoteo(a, i):
   fil = a.shape[0]
   aux_i = i
   
   for k in range(i + 1, fil):
      if abs(a[k, i]) > abs(a[aux_i, i]) :
         aux_i = k
   P = numpy.eye(fil)
   if aux_i != i :
      P[[aux_i, i]] = P[[i, aux_i]]
   return P

"""
   IN: U matriz superior, b vector de coeficientes
   OUT: x vector solución
"""
def resolverMTriangularSup(U, b):
   (fil, col) = U.shape
   x = numpy.zeros(col)
   for i in reversed(range(fil)):
      suma = U[i, i + 1:].dot(x[i + 1:]) # sumatoria
      x[i] = (b[i] - suma) / U[i, i]
   return x

"""
   Calcula la inversa de la matriz utilizando
   el metodo de gauss Jordan
   IN: a matriz cuadrada
   OUT: inversa de a
"""
def inversa(a, v = False):
   (fil, col) = a.shape
   assert fil == col, "La matriz no es cuadrada."
   assert numpy.linalg.det(a) != 0, "Matriz singular!"
   E = numpy.eye(fil)
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
            
      
      E_i = numpy.eye(fil) - alpha_i * e_i
      
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
   
"""
   La siguiente función recibe una matriz cuadrada A
   y la matriz b, junto a dim(A) y la trianguliza superior-
   mente mediante el método de eliminación de Gauss.
   Utiliza pivoteo parcial. Devuelve el vector de resultados
"""
def elimGauss(a, b, p = 'no', v = False):
   a_b = numpy.c_[a, b] # Matriz aumentada
   (fil, col) = a_b.shape
   
   if v :
      print("Matriz aumentada: ")
      print(a_b)

   if p == 'total' :
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
         print("Matriz a_b despupes del pivoteo inicial:")
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
            
      
      L_i = numpy.eye(fil) - alpha_i * e_i
      
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
