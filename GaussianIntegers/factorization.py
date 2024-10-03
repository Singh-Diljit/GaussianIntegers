"""Prime factorization over Z[i]."""

from GaussianInteger import GaussianInteger
import pickle

with open('primeDataZi/SOSpr_trivial.pickle', 'rb') as f:
    SOSList = pickle.load(f)

with open('primeDataZi/SOSpr_trivial_dict.pickle', 'rb') as f:
    SOSDict = pickle.load(f)

def makeSOSprime_Z(tup):
    """Reconstruct a prime from an entry of SOSList.

    Returns
    -------
    - : int : A prime number over Z.

    Example(s)
    ----------
    >>> (3, 0)
    >>> 3

    >>> (2, 1)
    >>> 5

    Notes
    -----
    The 'SOSList' is an ordered list of primes in Z, each prime is represented
    by a tuple of the form: (a, b) with a > 0 and b >=0. If b == 0, then a is
    a prime number. If b > 0 then a**2 + b**2 is a prime number.
    
    """
    return tup[0] if (tup[1] == 0) else tup[0]**2 + tup[1]**2

def primeDecomposition(N):
    """Return the prime decomposition over Z.

    Parameters
    ----------
    N : int : Number to factor.
    
    Returns
    -------
    primeDecomp : dict : Prime divisors along with their multiplicity.
    
    Example(s)
    ----------
    >>> primeDecomposition(360)
    >>> {2: 3, 3: 2, 5: 1}

    """
    if N == 0:
        return {}
    
    N = abs(N)
    primeDecomp = {}
    for tup in SOSList:
        pr = makeSOSprime_Z(tup)
        if pr**2 > N:
            primeDecomp[N] = 1
            break

        q, r = divmod(N, pr)
        prMult = 0
        while not r:
            N = q
            prMult += 1
            q, r = divmod(N, pr)

        if prMult > 0:
            primeDecomp[pr] = prMult
            
        if N == 1:
            break

    return primeDecomp

def primeDecomp(z):
    """Return the prime decomposition over Z[i].
    
    Returns
    -------
    res : dict : Prime divisors with their multiplicity, entries are tuples.
    
    Example(s)
    ----------
    >>> el = GaussianInteger(26, 18)
    >>> primeDecomp(el)
    >>> {(1, 1): 3, (2, 1): 3}
    
    """
    decompZ = primeDecomposition(z.norm)
    res = {}
    for p in decompZ:
        m = decompZ[p]
        if p == 2:
            res[(1,1)] = m
        elif p % 4 == 3:
            res[(p, 0)] = m//2
        else:
            p = GaussianInteger(vec=SOSDict[p])
            expP = 0
            while z.gcd(p) == p:
                expP += 1
                z = z // p
            if expP != 0:
                res[p.vec] = expP
            if expP != m:
                res[p.conjugate.vec] = m - expP
    return res

def decompUnit(z, decomp=False):
    """Return the unit associated with a prime decomposition over Z[i]."""
    if not decomp:
        decomp = primeDecomp(z)

    prod = 1
    for p in decomp:
        p, m = GaussianInteger(vec=p), decomp[p]
        prod *= p ** m
        
    return prod//z

def divisorsZi(z, prDecomp=False, proper=False, sort=False):
    """Return all divisors over Z[i].

    Parameters
    ----------
    z : GaussianInteger : Gaussian integer whose divisors will be generated.
    
    Returns
    -------
    div : list : List of all divisors over Z[i], entries are 'GaussianInteger'.
    
    Example(s)
    ----------
    >>> el = GaussianInteger(5, 5)
    >>> [str(x) for x in divisorsZi(el)]
    >>> ['2 - i', '5 + 5i', '2 + i', '3 + i', '1 + i', '5', '1', '1 + 3i']

    >>> [str(x) for x in divisorsZi(el, proper=True)]
    >>> ['2 - i', '2 + i', '3 + i', '1 + i', '5', '1 + 3i']

    >>> [str(x) for x in divisorsZi(el, proper=True, sort=True)]
    >>> ['1 + i', '2 - i', '2 + i', '3 + i', '1 + 3i', '5']
    
    """
    if not prDecomp:
        prDecomp = primeDecomp(z)

    div = {(1,0)}
    for pr in prDecomp:
        p, m = GaussianInteger(vec=pr), prDecomp[pr]
        for _ in range(m):
            toAdd = {(factor * p).vec for factor in div}
            div.update(toAdd)

    if proper:
        div.remove((1,0))
        div.remove(z.vec)
        
    div = [GaussianInteger(vec=r_) for r_ in div]
    return div if (sort==False) else sorted(div, key=lambda x: x.norm)
