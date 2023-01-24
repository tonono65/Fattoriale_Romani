import numpy as np
from  math import *
from numpy.fft import *



def DFT_matrix(N):
    i, j = np.meshgrid(np.arange(N), np.arange(N))
    test = i*j
    print(i)
    print()
    print(j)
    print()
    print(test)
    print()
    omega = np.exp( - 2 * pi * 1J / N )
    W = np.power( omega, i * j ) / sqrt(N)
    return W


def trova_prima_potenzadi_due_non_minore(n):
    potenza_di_due = 1
    while potenza_di_due < n:
        potenza_di_due *= 2
    return potenza_di_due



def propaga_resti_che_eccedono_BASE(vettore, base):
    # * propaga i resti
    #  
	riporto = 0
	vettore_propagato = [0 for i in range(len(vettore))] 
	for i in range(len(vettore)):
		x = vettore[i] + riporto
		riporto = x // base
		vettore_propagato[i] = x % base
	if riporto == 0:
		return vettore_propagato
	else:
		return vettore_propagato + [riporto]
		




def main():
    # print(DFT_matrix(5))
    
    '''
    n = 8 #numero dei campionamenti
    a = np.array([4,3,2,1])
    b = np.array([8,7,6,5])
    '''
    
    
    
    a = np.array([9,8,7,6,5,4,3,2,1])
    b = np.array([8,8,8,8,8,8,8,8])
    
    # n = 16 #numero dei campionamenti
    n = trova_prima_potenzadi_due_non_minore(2 * max(len(a), len(b))) 
    

    
    ya = n * ifft(a, n)
    yb = n * ifft(b, n)
    yc = ya * yb

    
    
    c = np.real(1./n * fft(yc))
    c = np.around(c, 0)
    c = propaga_resti_che_eccedono_BASE(c, 10)    

    print("\n" * 5)
    print("*************************************")
    print("*                                   *")
    print("*************************************")
    print("n: ", n, end='')
    print()
    print("Vettore a: ", end='')
    print(str(a))
    print()
    print("Vettore b: ", end ='')
    print(str(b))
    print()
    print("Vettore c: ", end='')
    print(str(c))

    print(" \n" * 2)
    print("Vettore dei coefficuenti ya: ", end = '')
    print(str(ya))
    print(" \n" * 1)
    print("Vettore dei coefficuenti yb: ")
    print(str(yb))
    print(" \n" * 1)
    print("Vettore dei coefficuenti yc: ")
    print(str(yc))
    print(" \n" * 1)
    
    # print("\n Matrice DFT(n): ")
    # print(str(DFT_matrix(n)))
    



main()