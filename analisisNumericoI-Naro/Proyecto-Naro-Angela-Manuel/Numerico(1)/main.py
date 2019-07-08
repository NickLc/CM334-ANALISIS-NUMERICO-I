# -*- coding: utf-8 -*-
import numericoUtils as eg
import factolu as lu
import numpy as np
import Cholesky_mod as ch
from time import time
import os
"""
   Este es el programa principal del proyecto
"""



def pressEnter():
    input("Presione 'Enter' para continuar...")
    os.system('clear')


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

pressEnter()

print("Cargando datos...")
a = np.array([[1, 1, 1, 1,],[2, 1, 3, 2],[1, 1, 1, 0],[1, 1, 0, 0]])
b1 = np.array([[-1], [-1], [0], [-1]])
b2 = np.array([[170], [320], [70], [70]])
x1 = np.array([-1, 0, 1, -1])
x2 = np.array([50, 20, 0, 100])
##SISTEMA 1
da1=np.array([[0.05,0.,0.,0.03],[0.,0.001,0.,0.],[0.,0.01,0.,0.],[0.,0.,0.,0.]])
db1=np.array([[0.],[0.004],[0.],[0.003]])
print(b1)
print(db1)
print(b1 + db1)
da2=np.array([[0.,0.,0.,0.],[0.3,0.,0.,0.],[0.5,0.1,0.,0.],[0.,0.,0.2,0.]])
db2=np.array([[0.01],[0.],[0.],[0.]])
da = [da1, da2]
db = [db1, db2]
b = [b1, b2]


print("Después de un arduo análisis, nuestros matemáticos(1) llegaron a la")
print("conclusión de que era necesario solo resolver dos sistemas de ecuaciones")
print("Estos son:")
print("* Sistema 1 [A|b1] :")
print(np.c_[a, b1])
print()
print("* Sistema 2 [A|b2] :")
print(np.c_[a, b2])
print()
print("Las soluciones exactas de los sistemas antes mencionados son, respectivamente: ")
print("Solución 1:")
print(x1)
print()
print("Solución 2:")
print(x2)

pressEnter()

print("Para probar la eficiencia de nuestros algoritmos, introducimos error")
print("en nuestras matrices. Las matrices de perturbación escogidas fueron:")
print("Para el sistema 1: ")
print("[A + dA|b + db]:")
print(np.c_[a + da1, b1 + db1])
print()
print("Para el sistema 2: ")
print("[A + dA|b + db]:") 
print(np.c_[a + da2, b2 + db2])
print()

print("Vamos a analizar este problema desde diferentes métodos.")  
print("A continuación escoja el método que creas conveniente para resolver el problema.")
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
for i in range(2):
    if op == 1 : # Eliminación de Gauss
        print("Eliminación de Gauss escogido.")
        print("Solución del sistema %d perturbado:" %(i + 1))
        t_inicio = time()
        print(eg.elimGauss(a + da[i], b[i] + db[i]))
        t_final = time()
        t_p = t_final - t_inicio
        print("Tiempo usado: %f." %(t_p))
    elif op == 2 : # Factorización LU
        print("Factorización LU escogido.")
        print("Solución del sistema %d perturbado:" %(i + 1))
        t_inicio = time()
        print(lu.sol(a + da[i], b[i] + db[i]))
        t_final = time()
        t_p = t_final - t_inicio      
        print("Tiempo usado: %f." %(t_p))
    elif op == 3 : # Cholesky
       print("Cholesky escogido.")
       print("Solución del sistema %d perturbado: " %(i + 1))
       if ch.posdef(a + da[i])==0:
           print("Matriz no definida positiva\nSe trabajará con la matriz transpuesta que multiplique al sistema de ecuaciones");
           t_inicio = time()
           print(ch.choleskyporTrans(np.mat(a + da[i]), np.mat(b[i] + db[i])))
           t_final = time()
           t_p = t_final - t_inicio
           print("Tiempo usado: %f." %(t_p))
       else:
           print("Matriz Definida Positiva")
           if ch.simetric(a + da[i]) == 0:
              print("Matriz no simetrica. Se trabajará el producto por su transpuesta\
              por la izquierda.")
              t_inicio = time()
              print(ch.choleskyporTrans(a + da[i], b[i] + db[i]))
              t_final = time()
              t_p = t_final - t_inicio
              print("Tiempo usado: %f." %(t_p))
           else:
              print("Matriz Simetrica")
              print(ch.choleskyTrans(a + da[i], b[i] + db[i]))
              t_final = time()
              t_p = t_final - t_inicio
              print("Tiempo usado: %f." %(t_p))
    pressEnter()

print("Fin del programa. Tenga un buen día :)")




