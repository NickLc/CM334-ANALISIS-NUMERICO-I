import numpy as np
import numericoUtils
from time import time

a=np.array([[1.,1.,1.,1.],[2.,1.,3.,2.],[1.,1.,1.,0.],[1.,1.,0.,0.]])
n=len(a)
#Sistema1
b1=np.array([[-1.],[-1.],[0.],[-1.]])
#Sistema2
b2=np.array([[170.],[320.],[70.],[70.]])

lb = [b1, b2]

#PERTURBACIONES
   

##SISTEMA 1 
da1=np.array([[0.05,0.,0.,0.03],[0.,0.001,0.,0.],[0.,0.01,0.,0.],[0.,0.,0.,0.]])
db1=np.array([[0.],[0.004],[0.],[0.003]])

da2=np.array([[0.,0.,0.5,0.3],[0.,0.05,0.,0.],[0.,0.7,0.,0.],[0.,0.6,0.,0.]])
db2=np.array([[0.],[0.1],[0.06],[0.005]])

##SISTEMA 2

da3=np.array([[0.,0.001,0.4,0.],[0.,0.,0.,0.],[0.,0.1,0.,0.],[0.,0.,0.1,0.]])
db3=np.array([[0.5],[0.],[0.6],[0.1]])

da4=np.array([[0.,0.,0.,0.],[0.3,0.,0.,0.],[0.5,0.1,0.,0.],[0.,0.,0.2,0.]])
db4=np.array([[0.],[0.],[0.],[0.]])

lda = [da1, da2, da3, da4]
ldb = [db1, db2, db3, db4]


def R(A, b, x):
   ax = np.matmul(A, x)
   r = ax - b
   return np.linalg.norm(r, np.inf)

def Cond(a, b, r, x, x_i):
   norm_b = np.linalg.norm(b, np.inf)
   condA = np.linalg.cond(a, np.inf)
   condA_inv = 1 / condA
   E = x - x_i
   E_r = np.linalg.norm(E, np.inf) / np.linalg.norm(x, np.inf)
   rb = r / norm_b
   return condA_inv * rb <= E_r and E_r <= condA * rb 


if __name__ == '__main__' :
   r = []
   for b in lb:
      print("-----------------------------------------------")
      print("Sistema integro:")
      print(np.c_[a, b])
      x0 = numericoUtils.elimGauss(a, b);
      print("Solución del sistema integro:")
      print(x0)
      print()
      print()
      r_temp = []
      for i in range(4):
         print("Sistema Perturbado %d:" %i)
         print(np.c_[a + lda[i], b + ldb[i]])
         t_inicio = time()
         x_i = numericoUtils.elimGauss(a + lda[i], b + ldb[i]);
         t_final = time()
         t_p = t_final - t_inicio
         print("Solución sistema perturbado")
         print("x_%d" %(i+1))
         print(x_i)
         print("Tiempo usado: %f." %t_p)
         r_temp.append(R(a, b, x_i))
         print("El residuo R es: ", r_temp[i])
         if Cond(a, b, r_temp[i], x0, x_i) :
            print("Es problema está bien condicionado.")
         else :
            print("El problema está mál condicionado.")
         print("------------------------------------------")
      r.append(r_temp)
   print("Perturbaciones:")
   print(r[0])
   print(r[1])

##Correr da1 db1 y da2 db2 para el sistema 1 y en cada caso enviar la norma de Ax-b donde x es la solucion perturbada
##Correr da3 db3 y da4 db4 para el sistema 2 y en cada caso enviar la norma de Ax-b donde x es la solucion perturbada
