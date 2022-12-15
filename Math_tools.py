import time 
import math
from Costanti import *


def testa_potenza_di_2(n):
    k = n
    while (k > 1):
        if (k % 2) == 1:
            raise Exception(f"FFT: {n} not power of two ") 
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




'''
	/**
     * performs a transform on a set of real values 
     *
     * @param data the set of data (the lenght must be a power of two)
     * @param transform the task to be accomplished<br>
     * 	transform == FFT.DIRECT  means direct transform<br>
     * 	transform == FFT.INVERSE means inverse transform
     * 
     */  
'''
def realFT(data, transform):
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
     * @param direct the task to be accomplished<br>
     * 	direct == true  means direct transform<br>
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


def togli_zeri_in_testa_vettore(vettore):
    """
    Bigint.normalize()
    """
    pass

'''
  {
		isZero = false;
		int i = A.length - 1;
		while (i >= 0 && A[i] == 0) i--;
		if (i == -1) {
			this.A = new byte[1];
			isZero=true;
		}
		else if (i == A.length - 1) return;
		byte[] res = new byte[i + 1];
		System.arraycopy(A, 0, res, 0, i + 1);
		A = res;
	}

'''

def simple_multyply(vettore1, vettore2):
    pass

'''
private static BigInt simpleMultiply(BigInt val1, BigInt val2) {
		if (val1.A.length < val2.A.length) return simpleMultiply(val2, val1);
		int[] acc = new int[val1.A.length + val2.A.length];
		for (int k = 0; k < val2.A.length; k++)
			for (int i = 0; i < val1.A.length; i++)
				acc[i + k] += val1.A[i] * val2.A[k];
		BigInt res = propagateCarries(acc);
		return res;
'''


def fft_multyply(vettore1, vettore2):
    pass


def multiply(vettore1, vettore2):
    if len(vettore1) < LIMITE or len(vettore2) < LIMITE:
        return simple_multiply(vettore1, vettore2)
    else:
        return fft_multiply.multiply(vettore1, vettore2)
	


'''
    /**
     * performs a transform on a set of complex values
     *
     * @param data the set of data (the lenght must be a power of two)
     * @param direct the task to be accomplished<br>
     * 	direct == true  means direct transform<br>
     * 	direct == false means inverse transform
     * 
     */  
	public static void four1(double[] data, boolean direct) {
		int n = data.length;
		testN(n);

		int mmax, m, j, istep, i, isign;
		double wtemp, wr, wpr, wpi, wi, theta, tempi, tempr;
		double temp;
		if (direct)
			isign = 1;
		else
			isign = -1;
		j = 1;
		for (i = 1; i < n; i += 2) {
			if (j > i) {
				temp = data[j - 1];
				data[j - 1] = data[i - 1];
				data[i - 1] = temp;
				temp = data[j];
				data[j] = data[i];
				data[i] = temp;
			}
			m = n >> 1;
			while (m >= 2 && j > m) {
				j -= m;
				m >>= 1;
			}
			j += m;
		}
		mmax = 2;

		while (n > mmax) {
			istep = mmax << 1;
			theta = isign * (2 * Math.PI / mmax);
			wtemp = Math.sin(0.5 * theta);
			wpr = -2.0 * wtemp * wtemp;
			wpi = Math.sin(theta);
			wr = 1.0;
			wi = 0.0;
			for (m = 1; m < mmax; m += 2) {
				for (i = m; i <= n; i += istep) {
					j = i + mmax;
					tempr = wr * data[j - 1] - wi * data[j];
					tempi = wr * data[j] + wi * data[j - 1];
					data[j - 1] = (data[i - 1] - tempr);
					data[j] = (data[i] - tempi);
					data[i - 1] = (data[i - 1] + tempr);
					data[i] = (data[i] + tempi);
				}
				wr = (wtemp = wr) * wpr - wi * wpi + wr;
				wi = wi * wpr + wtemp * wpi + wi;
			}
			mmax = istep;
		}
	}
	
'''

'''
    /**
     * performs a transform on a set of real values 
     *
     * @param data the set of data (the lenght must be a power of two)
     * @param transform the task to be accomplished<br>
     * 	transform == FFT.DIRECT  means direct transform<br>
     * 	transform == FFT.INVERSE means inverse transform
     * 
     */  
	public static void realFT(double[] data, int transform) {
		int n = data.length;
		testN(n);

		int i, i1, i2, i3, i4, np3;
		double wr, wi, wpr, wpi, wtemp, theta, c1 = 0.5, c2, h1r, h1i, h2r, h2i;

		theta = Math.PI / (n >> 1);
		if (transform == DIRECT) {
			c2 = -0.5;
			four1(data, true);
		} else {
			c2 = 0.5;
			theta = -theta;
		}
		wtemp = Math.sin(0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = Math.sin(theta);
		wr = 1.0 + wpr;
		wi = wpi;
		np3 = n + 3;
		for (i = 2; i <= (n >> 2); i++) {
			i4 = 1 + (i3 = np3 - (i2 = 1 + (i1 = i + i - 1)));
			h1r = c1 * (data[i1 - 1] + data[i3 - 1]);
			h1i = c1 * (data[i2 - 1] - data[i4 - 1]);
			h2r = -c2 * (data[i2 - 1] + data[i4 - 1]);
			h2i = c2 * (data[i1 - 1] - data[i3 - 1]);
			data[i1 - 1] = (h1r + wr * h2r - wi * h2i);
			data[i2 - 1] = (h1i + wr * h2i + wi * h2r);
			data[i3 - 1] = (h1r - wr * h2r + wi * h2i);
			data[i4 - 1] = (-h1i + wr * h2i + wi * h2r);
			wr = (wtemp = wr) * wpr - wi * wpi + wr;
			wi = wi * wpr + wtemp * wpi + wi;
		}
		if (transform == DIRECT) {
			data[0] = ((h1r = data[0]) + data[1]);
			data[1] = (h1r - data[1]);
			for (i = 0; i < n; i++)
				data[i] = (data[i] / Math.sqrt(n));
		} else {
			data[0] = (c1 * ((h1r = data[0]) + data[1]));
			data[1] = (c1 * (h1r - data[1]));
			four1(data, false);
			for (i = 0; i < n; i++)
				data[i] = (2 * data[i] / Math.sqrt(n));
		}
	}
	
'''

'''
ORA testa_potenza_di_2(n) 
public static void testN(int n) {
		int k = n;
		while (k > 1) {
			if (k % 2 == 1)
				throw new IllegalArgumentException("FFT: n not power of two "+ n);
			k /= 2;
		}
	}
'''
'''
public static BigInt multiply(BigInt val1, BigInt val2) {
		if(val1.A.length<LIM || val2.A.length<LIM)return simpleMultiply(val1, val2);
		return FFTMultiply.multiply(val1, val2);
	}
'''
'''
	private static BigInt simpleMultiply(BigInt val1, BigInt val2) {
		if (val1.A.length < val2.A.length) return simpleMultiply(val2, val1);
		int[] acc = new int[val1.A.length + val2.A.length];
		for (int k = 0; k < val2.A.length; k++)
			for (int i = 0; i < val1.A.length; i++)
				acc[i + k] += val1.A[i] * val2.A[k];
		BigInt res = propagateCarries(acc);
		return res;
	}
'''

'''
    public void normalize() {
		isZero = false;
		int i = A.length - 1;
		while (i >= 0 && A[i] == 0) i--;
		if (i == -1) {
			this.A = new byte[1];
			isZero=true;
		}
		else if (i == A.length - 1) return;
		byte[] res = new byte[i + 1];
		System.arraycopy(A, 0, res, 0, i + 1);
		A = res;
	}
'''

'''
**********************************************************************************************************	
public static BigInt add(BigInt val1, BigInt val2) {
		if(val1.A.length<val2.A.length)return add(val2,val1);
		int[] res = new int[val1.A.length+1];
		for (int i = 0; i < val2.A.length; i++) res[i]=val1.A[i]+val2.A[i];
		for (int i = val2.A.length; i < val1.A.length; i++) res[i]=val1.A[i];
		return propagateCarries(res);
	}

**********************************************************************************************************	
	private static BigInt propagateCarries(int[] AF) {
		BigInt res = null;
		byte[] A = new byte[AF.length];
		long carry = 0;
		for (int i = 0; i < AF.length; i++) {
			long x = ((long) AF[i]) + carry;
			carry = x / BASE;
			A[i] = (byte) (x % BASE);
		}
		if (carry == 0) {
			res = new BigInt(A);
 			res.normalize();
			return res;
		}		
		byte[] A1 = new byte[A.length + 1];
		System.arraycopy(A, 0, A1, 0,A.length);
		A1[A.length] = (byte) carry;
		return new BigInt(A1);
	}

}'''