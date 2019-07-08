# -*- coding: utf-8 -*-
"""Programa principal para presentación en la expociencia"""

import numpy, numericoUtils, os, Cholesky_mod as ch, factolu as lu
from time import time


P = ['A', 'B', 'C', 'D', 'E', 'F'] # Vector de paises
c = [2, 1, 3, 2, 1, 1] # Nuestro vector de costos
A = numpy.array([[1, 1, 1, 1, 1, 1],
                [2, 1, 3, 2, 1, 1],
                [1, 1, 1, 0, 0, 0],
                [1, 1, 0, 0, 1, 1]])
B = A[:, :4]
C = A[:, 4:]

b = numpy.array([[170],
                [320],
                [70],
                [70]])

def pausa():
    print("Presione una tecla para continuar:")
    input()



if __name__ == '__main__' :
    os.system('clear')
    print("Bienvenido.")
    print()
    print("Presentación del problema:")
    print("Una empresa exporta su producto a seis países A, B, C, D, E y F en \
cantidades anuales determinadas por el vector \
(x, y, z, u ,v, w) = (20, 20, 30, 70, 10, 20)")
    print("La empresa desea modificar un política de exportaciones de modod \
que la cantidad u exportada al país D sea lo menor posible. Ahora bien, \
la capacidad de producción de la empresa está limitada a 170 unidades de \
producto; los costos de producción son distintos según el país de destino \
(a causa de transporte), y son, respectivamente, de 2, 1, 3, 2, 1 y 1 \
unidad monetaria por unidad producida, el presupuesto de la empresa es de \
320 unidades monetarias; por último, los países A, B y C forman parte de \
una asociación económica que impone a la empresa una cuota de importación \
(máxima) de 70 unidades de producto anuales, y lo mismo sucede con los \
países A, B, E y F, con la misma cuota.")
    print("Determina las posibilidades de la empresa indicando los valores \
admisibles para los parámetros.")
    print()
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
    pausa()
    print("De los datos, tenemos la matriz de coeficientes:")
    print(A)
    print("La matriz de resultados será:")
    print(b)
    print("Al tener un sistema de 4 ecuaciones con 6 incógnitas, debemos \
        encontrar 4 variables en función de otras dos variables(parámetros).")
    print("Elegiremos a u y v como parámetros.")
    print()
    print("Así, planteamos las siguientes ecuaciones:")
    print("x = x_1 * v + x_2 * w + x_3")
    print("y = y_1 * v + y_2 * w + y_3")
    print("z = z_1 * v + z_2 * w + z_3")
    print()
    print("Así, tenemos los vectores")
    print("p_1 = (x_1, y_1, z_1, u_1)")
    print("p_2 = (x_3, y_3, z_3, u_3)")
    print("p_3 = (x_3, y_3, z_3, u_3)")
    print()
    pausa()

    print("Dvidiremos nuestra matriz de codeficientes de manera que \
    A = [B|C], donde ")
    print("B : \n", B)
    print("C : \n", C)
    print("En la ecuación inicial:")
    print("B * p_1 + C * e_1 + B * x_2 + C * e_2 + B * p_3 = b")




    print("De esta manera, conseguimos 3 sistemas con 4 ecuaciones y 4 incógnitas.")
    
    
    # Preparando sistema 1
    c1 = numpy.matmul(C, numpy.array([[1], [0]]))
    c1 = -c1
    print("Sistema 1:")
    print(numpy.c_[B, c1])
    
    # Preparando sistema 2
    c2 = numpy.matmul(C, numpy.array([[0], [1]]))
    c2 = -c2
    print("Sistema 2:")
    print(numpy.c_[B, c2])
    
    # Preparando sistema 3
    print("Sistema 3: ")
    print(numpy.c_[B, b])
    print()
    print("Escoja un método para solucionar los sistemas de ecuaciones.")
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

    os.system('clear') # aquí debería ir la transición.
    if op == 1 : # Eliminación de Gauss
        print("Has escogido eliminación de Gauss:")
        t_inicio = time()
        x1 = numericoUtils.elimGauss(B, c1);
        t_final = time()
        print("x1")
        print(x1)
        t_p = t_final - t_inicio      
        print("Tiempo usado: %f." %(t_p))
        print()
        print("x2")
        t_inicio = time()
        x2 = numericoUtils.elimGauss(B, c2);
        t_final = time()
        print(x2)
        t_p = t_final - t_inicio      
        print("Tiempo usado: %f." %(t_p))

        print("x3")
        t_inicio = time()
        x2 = numericoUtils.elimGauss(B, b);
        t_final = time()
        print(x2)
        t_p = t_final - t_inicio      
        print("Tiempo usado: %f." %(t_p))

    elif op == 2 : # Factorización LU
        print("Factorización LU escogido.")
        print("x1")
        t_inicio = time()
        print(lu.sol(B, c1))
        t_final = time()
        t_p = t_final - t_inicio      
        print("Tiempo usado: %f." %(t_p))

        print("x2")
        t_inicio = time()
        print(lu.sol(B, c2))
        t_final = time()
        t_p = t_final - t_inicio      
        print("Tiempo usado: %f." %(t_p))

        print("x3")
        t_inicio = time()
        print(lu.sol(B, b))
        t_final = time()
        t_p = t_final - t_inicio      
        print("Tiempo usado: %f." %(t_p))


    elif op == 3 : # Cholesky
       print("Cholesky escogido.")
       if ch.posdef(B) == 0:
            print("Matriz no definida positiva\nSe trabajará con la matriz transpuesta que multiplique al sistema de ecuaciones");
            print("x1")
            t_inicio = time()
            print(ch.choleskyporTrans(numpy.mat(B), numpy.mat(c1)))
            t_final = time()
            print("x2")
            t_inicio = time()
            print(ch.choleskyporTrans(numpy.mat(B), numpy.mat(c2)))
            t_final = time()
            t_p = t_final - t_inicio
            print("Tiempo usado: %f." %(t_p))

            print("x3")
            t_inicio = time()
            print(ch.choleskyporTrans(numpy.mat(B), numpy.mat(b)))
            t_final = time()
            t_p = t_final - t_inicio
            print("Tiempo usado: %f." %(t_p))
       else:
            print("Matriz Definida Positiva")
            if ch.simetric(B) == 0:
                print("Matriz no simetrica. Se trabajará el producto por su transpuesta\
                por la izquierda.")
                print("x1")
                t_inicio = time()
                print(ch.choleskyporTrans(B, c1))
                t_final = time()
                t_p = t_final - t_inicio

                print("x2")
                t_inicio = time()
                print(ch.choleskyporTrans(B, c2))
                t_final = time()
                t_p = t_final - t_inicio
                print("Tiempo usado: %f." %(t_p))

                print("x3")
                t_inicio = time()
                print(ch.choleskyporTrans(B, b))
                t_final = time()
                t_p = t_final - t_inicio
                print("Tiempo usado: %f." %(t_p))
            else:
                print("Matriz Simetrica")
                print("x1")
                t_inicio = time()
                print(ch.choleskyTrans(B, c1))
                t_final = time()
                t_p = t_final - t_inicio
                print("Tiempo usado: %f." %(t_p))

                print("x2")
                t_inicio = time()
                print(ch.choleskyTrans(B, c2))
                t_final = time()
                t_p = t_final - t_inicio
                print("Tiempo usado: %f." %(t_p))

                print("x2")
                t_inicio = time()
                print(ch.choleskyTrans(B, b))
                t_final = time()
                t_p = t_final - t_inicio
                print("Tiempo usado: %f." %(t_p))

    print("Fin del programa")
    
    
