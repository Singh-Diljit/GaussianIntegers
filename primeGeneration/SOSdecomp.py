"""Generate prime bank for use in constructing primes over Z[i]."""

from smallPrimes.primes import genPrimes, sumOfSquares
import pickle

def picklePrimeSOS(N, folderName='primeDataZi'):
    """Save primes under N with their sum of squares representation if possible.

    Parameters
    ----------
    N          : int : Upper bound of prime generation, N > 5.
    folderName : str : Folder data will be written to (must already exist).
    
    Writes
    ------
    SOSpr_trivial.pickle     : list : Ordered primes via SOS decomp (dtype=int).
    SOSpr_trivial_dict.pickle: dict : Keys are primes values are SOS decomp. 

    Notes
    -----
    Let P be the i'th prime (counting from 0, so 2 is the 0th prime.)
    Then:
        SOSpr_trivial[i] = (a, b), where P = a if (b==0) else a**2 + b**2.
        SOSpr_trivial_dict[P] = SOSpr_trivial[i] = (a, b).
        
    Examples
    --------
    >>> SOSpr_trivial[:5]
    >>> [(1,1), (3,0), (2,1), (7,0), (3, 2)]

    >>> SOSpr_trivial_dict[5]
    >>> (2,1)
    
    """
    path_ = lambda fileName: f'{folderName}/{fileName}.pickle'
    prBank = list(map(int, genPrimes(N)))
    SOSexpression = lambda p: (p, 0) if (p%4 == 3) else sumOfSquares(p)

    SOSpr_trivial = []
    SOSpr_trivial_dict = {}
    for p in prBank:
        SOSexp = SOSexpression(p)
        SOSpr_trivial.append(SOSexp)
        SOSpr_trivial_dict[p] = SOSexp
    
    with open(path_('SOSpr_trivial'), 'wb') as f: 
        pickle.dump(SOSpr_trivial, f, pickle.HIGHEST_PROTOCOL)

    with open(path_('SOSpr_trivial_dict'), 'wb') as f: 
        pickle.dump(SOSpr_trivial_dict, f, pickle.HIGHEST_PROTOCOL)
