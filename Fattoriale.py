# import time 
import sys
from Utils import *


def calcolo_fattoriale(a, b):
    x = a + 1
    y = a + 2
    while (y <= b and x < sys.maxsize / y):
        x *= y
        y += 1
    # se y 0 b+1 abbiamo calcolato tutto il fattoriale
    if (y == b + 1): 
        return x
    else:
        return 0
        # return BigInt.multiply(q(a, (a + b) / 2), q((a + b) / 2, b));


def computeF(n):
    # print(n)
    t0 = get_time()
    lp = calcolo_fattoriale(0, n)
    t1 = get_time()
    print("n = " + str(n) + "--> " + str(lp))
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
