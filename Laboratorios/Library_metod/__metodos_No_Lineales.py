#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import toolNick as tn
import numpy as np
import pandas as pd
import sympy as sy

import matplotlib.pyplot as plt
def solve_Newton(f, df, x_0, max_iter=100, tol=1e-6, v=False):
	"""
		Método Newton
		-------------------
		Halla una raíz de la función f usando derivadas mediante 
		el método de Newton por diferencias finitas.
		Argumentos
		----------
		f - Función, df - Derivada de la Funcion, max_iter = maximo numero de iteraciones, 
		tol (opcional) - Cota para el error absoluto de la x, v - verbosidad, ver el proceso
		Devuelve
		--------
		x - Raíz de f"""
	x_old = np.copy(x_0)
	x_new = np.copy(x_0)
	if v :
		name_col = ['x', 'f(x)', 'f\'(x)']	
		index = ['k_0']  
		data = [[x_0, f(x_0), df(x_0)]]

	for i in range(0, max_iter):

		x_new = x_old - f(x_old)/df(x_old)

		if np.abs(x_new-x_old)/np.abs(x_new) < tol:
			if v :
				df = pd.DataFrame(data,columns = name_col, index=index)
				print(df)
				print("\nConverge en iter:{}\nResultado:{}".format(i, x_new))
			return x_new

		if v:
			data.append([x_new, f(x_new), df(x_new)])
			index.append('k_{}'.format(i+1))
		
		x_old = np.copy(x_new) 	
	
	if v: print("Max. iteraciones alcanzado.\nResultado:{}".format(x_new))
	return x_new


def solve_NewtonMod(f, df, ddf, x_0, max_iter=100, tol=1e-6, v=False):
	if v:
		print("METODO DE RESOLUCION DE NEWTON: EC. NO LIN.")
	x_old = np.copy(x_0)
	x_new = np.copy(x_0)
	if v :
		name_col = ['x', 'f(x)', 'f\'(x)']	
		index = ['k_0']  
		data = [[x_0, f(x_0), df(x_0)]]

	for i in range(0, max_iter):

		x_new = x_old - (f(x_old)*df(x_old))/(df(x_old)**2 - f(x_old)*ddf(x_old))

		if np.abs(x_new-x_old)/np.abs(x_new) < tol:
			if v :
				df = pd.DataFrame(data,columns = name_col, index=index)
				print(df)
				print("\nConverge en iter:{}\nResultado:{}".format(i, x_new))
			return x_new

		if v:
			data.append([x_new, f(x_new), df(x_new)])
			index.append('k_{}'.format(i+1))
		
		x_old = np.copy(x_new) 	
	
	if v: print("Max. iteraciones alcanzado.\nResultado:{}".format(x_new))
	return x_new


def solve_Biseccion(f, a_0, b_0, tol=1e-6, v=False):
	"""
		Método de bisección
		-------------------
		Halla una raíz de la función f en el intervalo [a_0, b_0] mediante 
		el método de bisección.
		Argumentos
		----------
		f - Función, a_0 - Extremo inferior del intervalo, b_0 - Extremo superior del intervalo, 
		tol (opcional) - Cota para el error absoluto de la x, v - verbosidad, ver el proceso
		Devuelve
		--------
		x - Raíz de f en [a_0, b_0]"""

	a = a_0
	b = b_0
	c = (a*f(b) - b*f(a))/(f(b) - f(a))
	i = 0
	
	if v :
		name_col = ['a', 'c', 'b', 'f(c)']	
		index = [] 
		data = []

	while (b - a)/2 > tol:
				
		c = (a+b)/2
		if v :
			data.append([np.round(a, 10), np.round(c, 10), np.round(b, 10), (b - a)/2])
			index.append('k_{}'.format(i))

		if f(c) == 0:
			a = c ;  b = c
			break

		if np.sign(f(a)) == np.sign(f(c)) : a = c
		else : b = c

		i += 1

	if v: 
		df = pd.DataFrame(data, columns = name_col, index=index)
		print('\nData Frame - Metodo Biseccion\n\n', df)
		print("\nConverge en iter:{}\nResultado:{}".format(i, c))

	return np.round(c, 10)


def solve_Secante(f, x_0, max_iter=100, tol=1e-6, v=False):
	"""
		Método de Secante
		--------------------
		Halla una raíz de la función f a partir de un x_0 mediante 
		el método de bisección.
		Argumentos
		----------
		f - Función, x_0 - Valor inicial,  
		tol (opcional) - Cota para el error absoluto de la x, v - verbosidad, ver el proceso
		Devuelve
		--------
		x - Raíz de f """

	x_2old = np.copy(x_0)*1.01 # Jugar con este valor
	x_old = np.copy(x_0)
	x_new = np.copy(x_0)

	if v :
		name_col = ['x', 'f(x)']	
		index = [] 
		data = []

	for i in range(0, max_iter):
		if v :
			data.append([x_new, f(x_new)])
			index.append('k_{}'.format(i))

		x_new = x_old - (x_old-x_2old)*f(x_old)/(f(x_old)-f(x_2old))
		if np.abs(x_new-x_old) < tol:
			if v: 
				df = pd.DataFrame(data, columns = name_col, index=index)
				print('\nData Frame - Metodo Secante\n\n', df)
				print("\nConverge en iter:{}\nResultado:{}".format(i, x_new))
			return x_new

		x_2old = np.copy(x_old)
		x_old = np.copy(x_new)

	if v:
		print("Max. iteraciones alcanzado.\nResultado:{}".format(x_new))

	return x_new


def solve_ReguleFalsi(f, a_0, b_0, tol=1e-6, v=False):
	"""
		Método Regular Falsi
		-------------------
		Halla una raíz de la función f en el intervalo [a_0, b_0] mediante 
		el método de regular falsi o metodo posicion falsa.
		Argumentos
		----------
		f - Función, a_0 - Extremo inferior del intervalo, b_0 - Extremo superior del intervalo, 
		tol (opcional) - Cota para el error absoluto de la x, v - verbosidad, ver el proceso
		Devuelve
		--------
		c - Raíz de f en [a_0, b_0]"""

	a = a_0
	b = b_0
	c = (a*f(b) - b*f(a))/(f(b) - f(a))
	i = 0

	if v :
		name_col = ['a', 'c', 'b', 'f(c)']	
		index = [] 
		data = []

	while np.abs(f(c)) > tol:
		
		if v :
			data.append([a, c, b, f(c)])
			index.append('k_{}'.format(i))
			
		c = (a*f(b) - b*f(a))/(f(b) - f(a))

		if f(c) == 0:
			a = c ; b = c
			break

		if np.sign(f(a)) == np.sign(f(c)):
			b = c
		else:
			a = c
		i += 1

	if v: 
		df = pd.DataFrame(data, columns = name_col, index=index)
		print('\nData Frame - Metodo Regular falsi\n\n', df)
		print("\nConverge en iter:{}\nResultado:{}".format(i, c))

	return c


def solve_PuntoFijo(f, g, x_0, max_iter=100, tol=1e-6, v=False, graphic=False):
	"""
		Método de Punto Fijo
		--------------------
		Halla una raíz de la función f a partir de un x_0 mediante 
		el método de Punto Fijo.
		Argumentos
		----------
		f- Funcion, g - Función g(x) se obtiene al despejar una x de f(x) , x_0 - Valor inicial,  
		tol (opcional) - Cota para el error absoluto de la x, v - verbosidad, ver el proceso
		graphic - Mostrar grafica
		Devuelve
		--------
		x - Raíz de f 
		Ejemplo
		--------
		f(x) = g(x) - x\n
		f = lambda x: np.exp(-x) - x\n
		g = lambda x: np.exp(-x)\n
		"""
	
	if v :
		name_col = ['x', 'g(x)']	
		index = [] 
		data = []
	x = x_0
	x_new = g(x)
	values= [x]

	for i in range(0, max_iter):

		if v :
			data.append([x, x_new])
			index.append('k_{}'.format(i))

		if abs(x - x_new) < tol:
			if v: 
				df = pd.DataFrame(data, columns = name_col, index=index)
				print('\nData Frame - Metodo Punto Fijo\n\n', df)
				print("\nConverge en iter:{}\nResultado:{}".format(i, x_new))
			if graphic:
				a = -6. ; b = 10
				x = np.linspace(a, b)              # Vector de 50 elementos equiespaciados en (a,b)
				plt.plot(x, f(x),'g-',label='f(x)=x-g(x)')              
				plt.plot(x, g(x),'r-',label='g(x)')
				plt.plot(x, 0*x,'k-')              # Eje OX
				plt.plot(x, x,'b-',label='y = x')
				
				
				for i in range(0, len(values)):
					if i == 0:  # Valor inicial
						plt.plot(values[i],0,'ro',label=u'raíz')     
						plt.plot(values[i],values[i],'bo',label='punto fijo')

					plt.plot(values[i],0,'ro')     
					plt.plot(values[i],values[i],'bo')

				plt.legend(loc='best')
				plt.show()
			return x_new
		values.append(x_new)
		x = x_new
		x_new = g(x)

	if v:
		print("Max. iteraciones alcanzado.\nResultado:{}".format(x_new))

	return x_new	

# NL system solving methods
# Recordar poner antes de usar estos metodos:

def solve_Newton_Mod(x_0, f, J, max_iter, tol, v=False):
	"""
		Método Newton Modificado
		-------------------
		Halla una raíz de la función f usando metodo como 
		Jacobi, Gauss Seidel, SOR.
		Argumentos
		----------
		f - Función, J - Matriz Jacobiana de f, x_0 - Valor Inicial, max_iter - Maximo Numero de Iteraciones, 
		tol (opcional) - Cota para el error absoluto de la x, v - Verbosidad, ver el proceso
		Devuelve
		--------
		x - Raíz de f"""
	x_old = np.copy(x_0)
	x_new = np.copy(x_0)
	for i in range(0, max_iter):
		J_inv = tn.Inverse(J(x_old))
		x_new = x_old - np.matmul(J_inv, f(x_old))

		if np.linalg.norm(x_new-x_old)/np.linalg.norm(x_new) < tol:

			if v:
				print("Converge en iter:{}\nResultado:\n{}".format(i, x_new))

			return x_new
		x_old = np.copy(x_new)

	if v:
		print("Max. iteraciones alcanzado.\nResultado:\n{}".format(x_new))

	return x_new

# Calcular D: Tomar J -> D
def solve_NLS_Jacobi(x_0, f, D, max_iter, tol, v=False):
	return solve_Newton_Mod(x_0, f, D, max_iter, tol, v)

# Calcular L: Tomar J -> L
def solve_NLS_Gauss_Seidel(x_0, f, L, max_iter, tol, v=False):
	return solve_Newton_Mod(x_0, f, L, max_iter, tol, v)

# Calcular SOR_mat:
# Tomar J ->  L y D
# Tomar w
# SOR_mat = [ ( (1-w)/w ) D + L ]
def solve_NLS_SOR(x_0, f, SOR_mat, max_iter, tol, v=False):
	return solve_Newton_Mod(x_0, f, SOR_mat, max_iter, tol, v)


if __name__ == '__metodos_Iterativos' :
    print("__metodos_No_Lineales se ha importado correctamente.")

if __name__ == '__main__':
	f = lambda x : x**2 - x
	df = lambda x : 2 * x
	g = lambda x: x**0.5

	#print('Newton: {}'.format(solve_Newton(f, df, 4,v='True')))
	#print('Biseccion\n x = {}'.format(solve_Biseccion(f,0, 10,v='True')))
	#x = solve_Secante(f, 2,v='True')
	#x = solve_ReguleFalsi(f,0, 10,v='True')
	f = lambda x : x - np.sqrt(x+6) 
	g = lambda x : np.sqrt(x+6) 
	x = solve_PuntoFijo(f, g, 7, max_iter=100, tol=1e-6, v = True, graphic=True)
	print(x)