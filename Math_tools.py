import time 
import math
from Costanti import *
import scipy as sp
from  math import *
from scipy.fft import *



def propaga_resti_che_eccedono_BASE(vettore):
    # * propaga i resti
    #  
	riporto = 0
	vettore_propagato = [0 for i in range(len(vettore))] 
	for i in range(len(vettore)):
		x = vettore[i] + riporto
		riporto = x // BASE
		vettore_propagato[i] = x % BASE
	if riporto == 0:
		return vettore_propagato
	else:
		return vettore_propagato + [riporto]
		


def testa_potenza_di_2(n):
    # * verifica se n è una potenza di 2
    # *
    k = n
    while (k > 1):
        if (k % 2) == 1:
            raise Exception(f"FFT: {n} not power of two ") 
        k /= 2



def trova_prima_potenzadi_due_non_minore(n):
    potenza_di_due = 1
    while potenza_di_due < n:
        potenza_di_due *= 2
    return potenza_di_due



def BASE_to_10(vettore):
    vettore = togli_zeri_in_testa_vettore(vettore)
    if vettore[-1] < 10:
        return len(vettore) * 2 - 1
    else:
        return len(vettore) * 2



def togli_zeri_in_testa_vettore(vettore):
    i = len(vettore) - 1
    while i >= 0 and vettore[i] == 0:
        i -= 1
    if i == -1:
        return []
    elif i == len(vettore) - 1:
        return vettore
    else:
        return vettore[0:i+1]



# Trasforma un intero in un vettore in base = BASE con le cifre da sx [0] a dx [len-1]
def int_to_vettore(x):
    vettore = []
    if x == 0:
        return [1]
    else:
        n = 1
        z = x
        while (z // BASE) > 0:
            z = z // BASE
            n +=1
        for i in range(n):
            vettore.append(x % BASE)
            x = x // BASE
        return vettore

def multiply(vettore1, vettore2):
    if isinstance(vettore1, int):
        vettore1 = int_to_vettore(vettore1)
    if isinstance(vettore2, int):
        vettore2 = int_to_vettore(vettore2)
        
    if len(vettore1) < LIMITE or len(vettore2) < LIMITE:
        return simple_multiply(vettore1, vettore2)
    else:
        return new_fft_multiply(vettore1, vettore2)
	


def simple_multiply(vettore1, vettore2):
    if len(vettore1) < len(vettore2):
        return simple_multiply(vettore2, vettore1)

    acc = [0 for i in range(len(vettore1) + len(vettore2))]
    for k in range(len(vettore2)):
        for i in range(len(vettore1)):
            acc[i + k] += vettore1[i] * vettore2[k]
            acc = propaga_resti_che_eccedono_BASE(acc)
            # print(str(acc))
    return acc


def new_fft_multiply(vettore1, vettore2):
    n = trova_prima_potenzadi_due_non_minore(2 * max(len(vettore1), len(vettore2))) 
    

    ya = n * ifft(vettore1, n)
    yb = n * ifft(vettore2, n)
    yc = ya * yb
    
    vettore3 = sp.real(1./n * fft(yc))
    vettore3 = sp.around(vettore3, 0)
    vettore3 = propaga_resti_che_eccedono_BASE(vettore3)    
    return(vettore3)
    


    

# *****************************************************************
# FUNZIONI DI ROMANI (non utilizzate)
# *****************************************************************

def realFFT(data, transform):
    n = len(data)
    try:
        testa_potenza_di_2(n)
    except Exception as e:
        print(e)
        # int i, i1, i2, i3, i4, np3;
        # double wr, wi, wpr, wpi, wtemp, theta, c1 = 0.5, c2, h1r, h1i, h2r, h2i;
        
        c1 = 0.5
        theta = math.pi / (n >> 1)
        if transform == DIRECT:
            c2 = -0.5
            four1(data, DIRECT)
        else:
            c2 = 0.5
            theta = -theta
        wtemp = math.sin(0.5 * theta)
        wpr = -2.0 * wtemp * wtemp
        wpi = math.sin(theta)
        wr = 1.0 + wpr
        wi = wpi
        np3 = n + 3
        for i in range(2, n >> 2 + 1):
            i1 = i + i - 1
            i2 = 1 + i1
            i3 = np3 - i2
            i4 = 1 + i3
            h1r = c1 * (data[i1 - 1] + data[i3 - 1])
            h1i = c1 * (data[i2 - 1] - data[i4 - 1])
            h2r = -c2 * (data[i2 - 1] + data[i4 - 1])
            h2i = c2 * (data[i1 - 1] - data[i3 - 1])
            data[i1 - 1] = (h1r + wr * h2r - wi * h2i)
            data[i2 - 1] = (h1i + wr * h2i + wi * h2r)
            data[i3 - 1] = (h1r - wr * h2r + wi * h2i)
            data[i4 - 1] = (-h1i + wr * h2i + wi * h2r)
            wtemp = wr
            wr = wtemp * wpr - wi * wpi + wr
            wi = wi * wpr + wtemp * wpi + wi
        
        if (transform == DIRECT):
            h1r = data[0]
            data[0] = h1r + data[1]
            data[1] = (h1r - data[1])
            for i in range(n):
                data[i] = (data[i] / math.sqrt(n))
        else:
            h1r = data[0]
            data[0] = (c1 * (h1r + data[1]))
            data[1] = (c1 * (h1r - data[1]))
            four1(data, INVERSE)
            for i in range(n):
                data[i] = (2 * data[i] / math.sqrt(n))



'''
     * performs a transform on a set of complex values
     *
     * @param data the set of data (the lenght must be a power of two)
     * @param direct the task to be accomplished
     * 	direct == true  means direct transform
     * 	direct == false means inverse transform
'''

def four1(data, transform):
    n = len(data)

    try:
        testa_potenza_di_2(n)
    except Exception as e:
        print(e)

# int mmax, m, j, istep, i, isign;
#	double wtemp, wr, wpr, wpi, wi, theta, tempi, tempr;
#	double temp;
    if transform == DIRECT:
        isign = 1
    else:
        isign = -1
    j = 1
    for i in range(1, n, 2):
    #for (i = 1; i < n; i += 2) {
        if j > i:
            temp = data[j - 1]
            data[j - 1] = data[i - 1]
            data[i - 1] = temp
            temp = data[j]
            data[j] = data[i]
            data[i] = temp
        m = n >> 1
        while m >= 2 and j > m:
            j -= m
            m = m >> 1
        j += m
    mmax = 2;

    while n > mmax:
        istep = mmax << 1;
        theta = isign * (2 * math.pi / mmax)
        wtemp = math.sin(0.5 * theta)
        wpr = -2.0 * wtemp * wtemp
        wpi = math.sin(theta)
        wr = 1.0
        wi = 0.0
        for m in range(1, mmax, 2):
            for i in range(m, n +1, istep):
                j = i + mmax
                tempr = wr * data[j - 1] - wi * data[j]
                tempi = wr * data[j] + wi * data[j - 1]
                data[j - 1] = (data[i - 1] - tempr)
                data[j] = (data[i] - tempi)
                data[i - 1] = (data[i - 1] + tempr)
                data[i] = (data[i] + tempi)
            wtemp = wr
            wr = wtemp * wpr - wi * wpi + wr
            wi = wi * wpr + wtemp * wpi + wi
        mmax = istep

def multiply_interna(AF, BF):
    n = len(AF)
    if len(AF) != len(BF):
        raise Exception(f"Lunghezza AF diversa da lunghezza BF") 
        c = math.sqrt(n)
        AF[1] *= (c * BF[1])
        AF[1] *= (c * BF[1])
        for i in range(2, n, 2):
            ai = c * AF[i]
            ai1 = c * AF[i + 1]
            AF[i] = ai * BF[i] - ai1 * BF[i + 1]
            AF[i + 1] = ai1 * BF[i] + ai * BF[i + 1]
        realFFT(AF, INVERSE)
        A = [0 for i in range(len(AF))]
        riporto = 0
        for i in range(len(AF)):
            x = round(AF[i]) + riporto
            riporto = x / BASE
            A[i] = x % BASE
		
        if riporto == 0:
            return A
        else:
            return A + [riporto]



def fft_multiply(vettore1, vettore2):
    n = makePowerOfTwo(max(len(vettore1), len(vettore2))) * 2
    AF = vettore1 + [0 for i in range(n-len(vettore1))]
    BF = vettore2 + [0 for i in range(n-len(vettore2))]
    realFFT(AF, DIRECT)
    realFFT(BF, DIRECT)
    return multiply_interna(AF, BF) 
