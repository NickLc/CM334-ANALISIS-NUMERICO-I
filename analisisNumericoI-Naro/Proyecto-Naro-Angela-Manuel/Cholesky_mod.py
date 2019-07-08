import numpy as np

def simetric(A):
    M=0
    n = A.shape[0]
    for i in range(0,n-1):
        for j in range(0,n-1):
            if A[i,j]!=A[j,i]:
                M=M-1
    if M==0:
        return 1
    else:
        return 0
def Pivot_mat(A, c_index):
    dim = A.shape[0]
    P = np.eye(dim)

    index = c_index
    for i in range(c_index + 1, dim):
        if abs(A[i,c_index]) > abs(A[index,c_index]):
            index = i

    if index != c_index:
        P[[index, c_index]] = P[[c_index, index]]

    return P
def Inverse(A):
    n_row, n_col = A.shape[0], A.shape[1]
    if n_row != n_col or np.linalg.det(A) == 0:
        print("Singular matrix or dimensions error")
        return

    A_inv = np.eye(n_row)


    for i in range(0, n_row): 
       
        T_i = np.eye(n_row)

        alpha_i = np.empty(n_row)
        alpha_i[:] = (A[:,i].T/A[i,i])
        alpha_i[i] = 1/A[i,i]

        e_i = np.zeros(n_row)
        e_i[i] = 1


        T_i = np.eye(n_row) - np.outer(alpha_i, e_i)

        T_i[i,i] = 1/A[i,i]


        A = np.matmul(T_i, A)
        A_inv = np.matmul(T_i, A_inv)



    return A_inv

def posdef(A):
    eig_values = np.linalg.eig(A)[0]
    M=0
    for i in eig_values:
	    if i < 0:
	        M=M-1
	       
    if M<0:
        return 0;
    else:
        return 1;
    
def CholeskyG(A, pivot="none"):
	eig_values = np.linalg.eig(A)[0]
	n = A.shape[0]
	G = np.zeros((n,n))
	P = np.eye(n)
	process_string = "A"


	for k in range(0, n):
		P_k = np.eye(n)

		if (A[k,k] == 0 and pivot == "partial") or pivot=="total":
			P_k = Pivot_mat(A, k)

		A = np.matmul(P_k, A)
		A = np.matmul(A, P_k.T)
		P = np.matmul(P_k, P)
  
		pre_sum = G[k, 0:k].dot(G[k, 0:k])
		G[k,k] = np.sqrt(A[k,k] - pre_sum)

		for i in range(k+1, n):
		    pre_sum_2 = G[i, 0:k].dot(G[k,0:k])
		    G[i,k] = (A[i,k] - pre_sum_2)/G[k,k]
            
            		          	          
	return G


def choleskyporTrans(A,b):
    
    
    Atr=A.T
    valor=Atr*A
    #print( Atr)
    Atrb=Atr*b
    #print(Atrb)
    #print(Atr*A)
    
    G=CholeskyG(valor)
    Gtr=G.T
    #print(G)
    #print(Gtr)
    
    #print(G*Gtr)
     #E*Etr*x=Atr*b
     #E*y=Atr*b
     
     #Definimos Y
    Ginv=Inverse(G)
    Gtrinv=Inverse(Gtr)
    #print(Ginv)
    GinvAtrb=Ginv*Atrb
    GtrinvGinvAtrb=Gtrinv*GinvAtrb
    #print(GtrinvGinvAtrb)
    return GtrinvGinvAtrb

def choleskysinTrans(A,b):
     
    G=CholeskyG(A)
    Gtr=G.T
    #print(G)
    #print(Gtr)
    
    #print(G*Gtr)
     #E*Etr*x=Atr*b
     #E*y=Atr*b
     
     #Definimos Y
    Ginv=Inverse(G)
    Gtrinv=Inverse(Gtr)
    #print(Ginv)
    Ginvb=Ginv*b
    GtrinvGinvb=Gtrinv*Ginvb
    #print(GtrinvGinvAtrb)
    return GtrinvGinvb


###def main():
    #A=np.array([(2, 1, 0, 4),(-4, -2, 3, -7),(4, 1, -2, 8),(0, -3, -12, -1)])
    #b=np.array([(2,),(3,),(4,),(5,)])
    #Amat=np.mat(A)
    #bmat=np.mat(b)           
    #print(cholesky(Amat,bmat))
    
#main()!##