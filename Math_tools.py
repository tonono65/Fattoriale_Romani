import time 





'''
public static void testN(int n) {
		int k = n;
		while (k > 1) {
			if (k % 2 == 1)
				throw new IllegalArgumentException("FFT: n not power of two "+ n);
			k /= 2;
		}
	}
'''
def testa_potenza_di_2(n):
    k = n
    while (k > 1):
        if (k % 2) == 1:
            pass
            # qui bisogna gestire l'ecczione
        k /= 2



'''
	/**
     * find the first power of two not smaller than n
     *
     * @param n the integer 
     * @return  the power of two
     */  
	protected static int makePowerOfTwo(int max) {
		int twoPower = 1;
		for (; twoPower < max; twoPower *= 2);
		return twoPower;
	}

'''

def trova_prima_potenzadi_due_non_minore(n):
    potenza_di_due = 1
    while potenza_di_due < n:
        potenza_di_due *= 2
    return potenza_di_due