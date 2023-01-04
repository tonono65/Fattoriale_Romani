# import time 
import sys
from Costanti import *
from Utils import *
from Math_tools import * 

# in java funzione q()
def calcolo_fattoriale(a, b):
    """ 
    in java funzione q()
    a = 0
    b = valore di cui calcolare il fattoriale
    """
    x = a + 1
    y = a + 2
    while (y <= b) and (x < sys.maxsize / y):
        x *= y
        y += 1
    # se y 0 b+1 abbiamo calcolato tutto il fattoriale
    if (y == b + 1): 
        return x
    else:
        return multiply(calcolo_fattoriale(a, (a + b) // 2), calcolo_fattoriale((a + b) // 2, b))
        # return BigInt.multiply(q(a, (a + b) / 2), q((a + b) / 2, b));
        


def computeF(n):
    # print(n)
    t0 = get_time()
    fattoriale = calcolo_fattoriale(0, n)
    if isinstance(fattoriale, list):
        digits = BASE_to_10(fattoriale)
    else:
        digits = len(str(fattoriale))
    t1 = get_time()
    print("n = " + str(n) + "--> " + str(digits))
    # print("n = " + n + ", digits = " + lp+", time = " + Util.askTime() + "sec.");
    time_elapsed(t0["tms"], t1["tms"])


def main():
    n = 10
    while n <= 10000000:
        computeF(n)
        
        n *= 10



main()

'''
if __name__ == "__main__":
    n = 10
    while n <= 10000000:
        computeF(n)
        
        n *= 10
'''

# ============================= 8-proc ===================================================================
#	n =    10000, digits =   35660, time =    0.145sec., memory 0.0M, max memory 0M
#	n =   100000, digits =   456574, time =   0.533sec., memory 0.0M, max memory 0M
#	n =  1000000, digits =  5565709, time =  22.055sec., memory 2.0M, max memory 2M
#	n = 10000000, digits = 65657060, time = 564.378sec., memory 31.0M, max memory 31M
# ========================================================================================================
