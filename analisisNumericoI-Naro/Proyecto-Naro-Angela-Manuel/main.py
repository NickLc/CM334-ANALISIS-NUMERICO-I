# -*- coding: utf-8 -*-
import numericoUtils
import numpy
from time import time
"""
   Este es el programa principal del proyecto
"""

P = ['A', 'B', 'C', 'D', 'E', 'F'] # Vector de paises
c = [2, 1, 3, 2, 1, 1] # Nuestro vector de costos

print("Bienvenido.")

# Luego se hará más automático y dinámico.
print("Tenemos:")
print("+------+----------------+")
print("| País | Costo de envío | ")
print("+------+----------------+")
print("|   A  |        2       |")
print("+------+----------------+")
print("|   B  |        1       |")
print("+------+----------------+")
print("|   C  |        3       |")
print("+------+----------------+")
print("|   D  |        2       |")
print("+------+----------------+")
print("|   E  |        1       |")
print("+------+----------------+")
print("|   F  |        1       |")
print("+------+----------------+")
print()
print("Inserte la matriz asociada al sistema de ecuaciones:")
while True :
   try:
      n = int(input("Ingrese la dimensión de las filas: "))
      m = int(input("Ingrese la dimensión de las columnas: "))
      assert n >= 0 or m >= 0
      assert n <= 6 or m <= 6
      break
   except AssertionError:
      print("Oh! Dimensiones no válidas. Intente nuevamente.")
      print("Para este problema, existen como máximo 6 variables.")
a = numpy.empty((n, m)); # Esta será nuestra matriz
print("Ahora ingresará los coeficientes:")
for i in range(n):
   for j in range(m):
      a[i, j] = float(input("[%d, %d]" %(i, j)))

print("La matriz ingresada es:")
print(a)
print()
print("Introducir la matriz de resultados(b):")
b = numpy.empty(n)
for i in range(n):
   b[i] = float(input("b[%d]: " %i))
print("Escoge el método que creas conveniente para \nresolver el problema.")
print("1. Eliminación de Gauss     2. Factorización LU")
print("3. Cholesky")

while True:
   try:
      op = int(input("opción: "))
      assert op <= 3 and op > 0
      break
   except AssertionError:
      print("Inserte una opción válida.")
   except ValueError:
      print("No seas ganso.")

if op == 1 : # Eliminación de Gauss
   if n != m :
      print("La matriz no es cuadrada.")
      print("Se hará el producto por su transpuesta\
      por la izquierda.")
      a_s = a
      a = numpy.matmul(a.transpose(), a)
      print(a_s)
      print(a)
   else:
      t_inicio = time()
      
      t_final = time()
      t_p = t_final - t_inicio
      print("Tiempo usado: %f." %(t_p))
elif op == 2 : # Factorización LU
   if n !=m : 
      print("Matriz no cuadrada.")
      print("Se hará el producto por su transpuesta\
      por la izquierda.")
   else:
      t_inicio = time()
      x = elimGauss(a, b)
      t_final = time()
      t_p = t_final - t_inicio
      print("El vector respuesta es x:")
      print(x)
      print("Tiempo usado: %f." %(t_p))
elif op == 3 : # Cholesky
   if n != m :
      print("La matriz no es cuadrada.")
      print("Se trabajará el producto por su transpuesta\
      por la izquierda.")
   else:
      print("Cholesky")

print("Fin del programa")





"""
   Escrito por Los Magníficos
"""
