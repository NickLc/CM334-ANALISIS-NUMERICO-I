#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import toolNick as tn
import numpy as np
import pandas as pd


# Interpolation methods: 

# Ln_i: Retorna el coeficiente de Lagrange L_ni
def Ln_i(n, i, x_vector, x):
	p1 = 1
	p2 = 1

	for j in range(0, n):
		if i != j:
			p1 *= (x - x_vector[j])
			p2 *= (x_vector[i] - x_vector[j])
	
	return p1/p2

# INT_Lagrange: Retorna el polinomio generado P_n_eq, listo para evaluar cualquier
# x y tomar el valor de la curva interpolada en cualquier punto
def INT_Lagrange(x_vector, y_vector, verbose=False):
	"""
		Metodo de Lagrange
		-------------------
		Este metodo genera una funcion P_n(x) 
		representando al polinomio interpolado por Lagrange.

		Argumentos
		-----------
		x_vector - Valores de la medición.\n
		y_vector - Valores encontrados a partir de x_vector, y = f(x).\n
		x_vector y v_vector son datos experimentales
		Devuelve
		----------

		"""
	def P_n_eq(x, verbose=False):

		if verbose:
			print("-----------------------------------")
			print("Evaluando x = {}".format(x))

		n = x_vector.size
		acc = 0
		for i in range(y_vector.size):

			if verbose:
				print("L_{},{}:{}".format(n, i, Ln_i(n, i, x_vector, x)))	

			acc += y_vector[i]*Ln_i(n, i, x_vector, x)

		if verbose:
			print("Finalmente, P_{}({}):{}".format(n, x, acc))
			print("-----------------------------------")

		return acc
	return P_n_eq

# INT_Lagrange_error: Retorna el polinomio de error generado P_n_eq,
# listo para evaluar cualquier c dentro del intervalo [a,b] en el
# que pertenecen los valores de x_vector
def INT_Lagrange_error(x_vector, dfn_plus_1, c, verbose=True):
	if verbose:
		print("METODO DE ERROR DE LAGRANGE")
		print("Este metodo genera una funcion E_n(x)\nrepresentando al polinomio de error de la funcion interpolada por Lagrange")
		print("Se usara los datos de \'x\':", x_vector)
		print("Se usara como c en el interv. [\'x\']:", c)
	if (c < x_vector[0]) or (c > x_vector[x_vector.size-1]):
			print("Advertencia: El valor de c ingresado no esta dentro de los limites de x")

	def E_n(x):

		if verbose:
			print("-----------------------------------")
			print("Evaluando x = {}".format(x))


		acc1 = 1 # product of (x-x_i)'s
		acc2 = 1 # (n+1)!
		n = x_vector.size

		for i in range(0,n):
			acc1 *= (x - x_vector[i])
			acc2 *= i+1
		
		result = (acc1 * dfn_plus_1(c))/acc2

		if verbose:
			print("Finalmente, (x-x_0)...(x-x_{})*f^({})({}) / {}! = {}".format(n-1, n, c, n,result))
		
		return result
	return E_n

def divdif(f, x_vector):
	n = x_vector.size
	if n is 1:
		return( f(x_vector) )
	else:
		x_copy_1 = np.copy(x_vector[1:])
		x_copy_2 = np.copy(x_vector[:n-1])
		return ( divdif(f, x_copy_1) - divdif(f, x_copy_2) )/(x_vector[x_vector.size-1] - x_vector[0])

def Newton_interpol_table(f, x_vector, verbose=False):
	n = x_vector.size
	table = np.empty((n, n)) # Tabla a generar
	table_aux = np.empty((n, 2))


	for i in range(table_aux.shape[0]):
		table_aux[i, 0] = i
		table_aux[i, 1] = x_vector[i]

	for i in range(table.shape[0]):
		for j in range(table.shape[1]):
			if (j <= i):
				table[i, j] = divdif(f, x_vector[i-j: i+1])
			else:
				table[i, j] = np.nan

	if verbose:
		print("Tabla de coeficientes via diferencia dividida:")
		table_head = "k\tx_k\t"
		for i in range(table.shape[1]):
			table_head += "f[x_k-{} -> x_k] | ".format(i)
		print(table_head)
		for i in range(table.shape[0]):
			print(table_aux[i, 0], '|', table_aux[i, 1], '|', end='')
			for j in range(table.shape[1]):
				print(table[i, j], "|\t", end='')
				if np.isnan(table[i, j]):
					print("\t\t", end='')
			print()

	return table

def INT_Newton(f, x_vector, verbose=False, my_coefs=None):

	if verbose:
		print("METODO DE NEWTON")
		print("Este metodo genera una funcion P_n(x)\nrepresentando al polinomio interpolado por Newton")
		print("Se usara los datos de \'x\':", x_vector)

	
	if my_coefs is None:
	
		if verbose:
			print("Coeficientes no proporcionados.")
			print("Se generaran coefs. via diferencias divididas.")
	
		coefs = Newton_interpol_table(f, x_vector, verbose)
		coefs = np.diag(coefs)
	
	else:
		coefs = np.copy(my_coefs)


	if verbose:
		print("Los coeficientes de la curva interpolada seran: ")
		print(coefs)


	str_eq = ""
	str_eq += str(coefs[coefs.size - 1])

	for i in range(coefs.size-1):
		str_eq = "(x - {}).[".format(x_vector[x_vector.size - i - 2]) + str_eq
		str_eq = "{} + ".format(coefs[coefs.size - i - 2]) + str_eq
	for i in range(coefs.size-1):
		str_eq += ']'
	str_eq = "P_{}(x) = ".format(coefs.size - 1) + str_eq


	def P_n_eq(x):
		acc = coefs[coefs.size - 1]

		for i in range(coefs.size-1):
			acc *= (x - x_vector[x_vector.size - i - 2])
			acc += coefs[coefs.size - i - 2]

		return acc
	
	if verbose:
		print("---------------------------------------------------------")
		print("Polinomio de interpolacion generado:\n{}".format(str_eq))
		print("---------------------------------------------------------")
		
	return P_n_eq


if __name__ == '__metodos_Iterpolacion' :
    print("__metodos_Interpolacion se ha importado correctamente.")


if __name__ == '__main__':

    f = lambda x : x**2 - 4
    df = lambda x : 2 * x
