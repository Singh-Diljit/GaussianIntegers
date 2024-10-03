"""Functions in support of the GaussianInteger class."""

def gcd_Z(n, k):
    """Return the GCD of two integers via the Euclidean algorithm."""
    while k:
        n, k = k, n%k

    return abs(n)

def signedGCD_Z(n, k):
    """Return the (signed) GCD of two integers via the Euclidean algorithm.

    Example(s)
    ----------
    >>> signedGCD_Z(-5, 10)
    >>> 5

    >>> signedGCD_Z(-5, -10)
    >>> -5

    Notes
    -----
    The signed GCD has the same magnitude as the usual GCD and is negative
    iff both inputs are negative.

    """
    return gcd_Z(n, k) if not(n < 0 and k < 0) else -gcd_Z(n, k)

def multGCD_Z(*els):
    """Return the GCD of multiple integers via the Euclidean algorithm.

    Example(s)
    ----------
    >>> multGCD_Z(2,6,4)
    >>> 2

    """
    nums = sorted(els, reverse=True)
    x = nums.pop()
    while nums:
        y = nums.pop()
        x = gcd_Z(x, y)
        
    return x

def euclidean_normSq(tup):
    """Return the sum of squared components for a tuple."""
    return sum(map(lambda x: x**2, tup))

def powsOfi_tup(k, scale=1):
    """Return scale*(i)**k as a tuple of the form: (real, imag)."""
    return [(scale, 0), (0, scale), (-scale, 0), (0, -scale)][k%4]

def nonTrivial_binaryPower(zVec, k):
    """Raise an imaginary number represented as a tuple to an integer power.
    
    Parameters
    ----------
    zVec: tuple : Represents zVec[0] + i*zVec[1], both components are non-zero.
    k   : int   : Greater than 1.

    Returns
    -------
    res: tuple : Associated tuple of resulting complex number.

    See Also
    --------
    'implementPower' is the same function with no assumptions on non-zero
    components or restricting k to be > 2.    

    """
    resRe, resIm = 1, 0
    curRe, curIm = zVec
    for i in bin(k)[::-1]:
        if i == '1':
            resRe, resIm = resRe*curRe - resIm*curIm, resRe*curIm+resIm*curRe
        curRe, curIm = curRe**2 - curIm**2, 2*curRe*curIm

    res = (resRe, resIm)
    return res

def implementPower(zVec, k):
    """Raise an imaginary number represented as a tuple to an integer power.

    Parameters
    ----------
    zVec: tuple : Represents zVec[0] + i*zVec[1].
    k   : int   : Power.

    Returns
    -------
    res: tuple : Associated tuple of resulting complex number.

    """
    if zVec[1] == 0:
        res = (zVec[0]**k, 0)
    elif k in {0, 1}:
        res = [(1, 0), zVec][k]
    elif zVec[0] == 0:
        mag_ = zVec[1]**k
        res = powsOfi_tup(k, scale=mag_)
    else:
        res = nonTrivial_binaryPower(zVec, k)

    return res

def roundNearest(x):
    """Round to the nearest integer, ties settled via min Z[i]-norm.

    Example(s)
    ----------
    >>> roundNearest(-2.5)
    >>> -2

    >>> roundNearest(2.5)
    >>> 2

    """
    res = round(x)
    if x < 0 and (res + .5 == x):
        res += 1
        
    return res

def roundVec(x):
    """Round components to the nearest integer, ties settled via min Z[i]-norm."""       
    return tuple([roundNearest(x_) for x_ in x])

def makeString_GaussianInteger(zVec):
    """Represent a Gaussian integer represented as a tuple as a string.

    Example(s)
    ----------
    >>> makeString_GaussianInteger((0, -1))
    >>> -i

    >>> makeString_GaussianInteger((1,1))
    >>> 1 + i

    >>> makeString_GaussianInteger((9,5))
    >>> 9 + 5i

    """
    re, im = zVec
    if im == 0:
        res = str(re)
    elif re == 0: #note im != 0
        if abs(im) == 1:
            res = 'i' if (im > 0) else '-i'
        else: res = str(im) + 'i'
    else:
        sgn = '-' if (im < 0) else '+'
        str_imag = abs(im) if (abs(im) != 1) else ''
        res = f'{str(re)} {sgn} {str_imag}i'
        
    return res

def proj_to_2tuple(X):
    """Return a 2-tuple of ints from data given.

    Parameters
    ----------
    X : * : Data to be projected to 2-tuple.

    Notes
    -----
    * Accepted objects are: int, float, or have a __getitem__ attribute.
    The choice to map an object with __getitem__ but an IndexError on its
    first element to (0, 0) is to encourage easy compatibility with other
    arithmetic-focused classes. Classes
    
    Example(s)
    ----------
    >>> proj_to_2tuple(-2.5)
    >>> (-2, 0)

    >>> proj_to_2tuple((3.0, 5))
    >>> (3, 5)

    >>> proj_to_2tuple((1,))
    >>> (1, 0)

    >>> proj_to_2tuple((0,0))
    >>> (0, 0)
    
    """
    res = 'failed'
    if isinstance(X, int):
        res = (X, 0)
    elif isinstance(X, float):
        res = (roundNearest(X), 0)
    elif X.__class__.__name__ == 'GaussianInteger':
        res = X.vec
    else:
        try:
            res = (roundNearest(X[0]), roundNearest(X[1]))
        except IndexError:
            try: res = (roundNearest(X[0]), 0)
            except IndexError: res = (0, 0)

    if res == 'failed':
        raise TypeError('Cannot project X to a 2-tuple.')
        
    return res
