"""Implement 'GaussianInteger' class."""

from helperFuncs import *

class GaussianInteger:
    """Implements elements of the form: p + qi for p, q in Z and i**2 = -1."""

    def __init__(self, real=0, imaginary=0, vec=None):
        """Initialize a Gaussian integer.

        Parameters
        ----------
        real      : int   : Real part of Guassian integer.
        imaginary : int   : Imaginary part of Guassian integer.
        vec       : tuple : 2-tuple of ints, used for backend initialization.

        Initializes
        -----------
        self.vec  : tuple : 2-tuple of integers, denoting real and imag parts.
        self.real : int   : Real part of Guassian integer.
        self.imag : int   : Imaginary part of Guassian integer.

        Notes
        -----
        Initialization via 'vec' does not sanitize inputs.
        
        """
        if vec != None:
            vec_ = roundVec(vec)
        else:
            vec_ = roundVec((real, imaginary))
            
        self.vec = vec_
        self.real, self.imag = self.vec

    # Basic Classification, Representation, and Accessing
    
    @property
    def isReal(self):
        """Return if self is purely real."""
        return self.imag == 0

    @property
    def isPurelyImag(self):
        """Return if self is purely imaginary."""
        return self.real == 0 and self.imag != 0

    @property
    def isAdditiveID(self):
        """Return self is the additive identity."""
        return self.vec == (0, 0)

    @property
    def isMultiplicativeID(self):
        """Return self is the multiplicative identity."""
        return self.vec == (1, 0)

    def __bool__(self):
        """Return if self is not the additive identity."""
        return not(self.isAdditiveID)

    def __key(self):
        """Return a key uniquely associated with self.

        Returns
        -------
        - : int, tuple : Is uniquely associated with self.vec.

        Notes
        -----
        If self has a natural interpretation as an 'int', the key returns
        the associated int. In all other cases: self.__key() == self.vec.

        This allows a more dynamic implementation of hashes useful in
        mixed class settings e.g. primes over Z[i].
        
        See Also
        --------
        See the 'unitKey' method for to identify Gaussian integers up to unit.
        
        """
        return self.real if self.isReal else self.vec
    
    @property
    def __hash__(self):
        """Return hash(self)."""
        return hash(self.__key())

    def __repr__(self):
        """Return repr(self).

        Notes
        -----
        This is the minimal data required to initialize an equivalent
        instance of the GaussianInteger class.
        
        """
        return f'GaussianInteger(real={self.real}, imaginary={self.imag})'

    def __str__(self):
        """Return self as a (simplified) string.

        See Also
        --------
        See the 'makeString_GaussianInteger' function for implementation
        details.
        
        """
        return makeString_GaussianInteger(self.vec)

    def __int__(self):
        """Return self.real."""
        return self.real

    def __getitem__(self, index):
        """Allow access to real (index=0) and imaginary (index=1) components."""
        return self.vec[index]

    # Standard Complex Number Features

    @property
    def slope(self):
        """Return the slope when interpreting self as a vector.

        Returns
        -------
        res : float : Steepness of self w.r.t. the real-axis.

        Notes
        -----
        This mirrors the role of 'theta' in the polar-coordinate
        expression of a complex number.

        The slope of point 0 is defined to be 'inf'.

        See Also
        --------
        For an exact version see the 'slope_exact' method.

        """
        if self.real == 0:
            res = float('inf') if (self.imag >= 0) else -float('inf')
        else:
            res = self.imag/self.real

        return res

    @property
    def gcdComponents(self):
        """Return the (signed) GCD of self.real and self.imag.

        Returns
        -------
        - : int : The (signed) GCD.

        See Also
        --------
        See the 'signedGCD_Z' function for examples and implementation details.

        """
        return signedGCD_Z(self.real, self.imag)

    @property
    def slopeExact(self):
        """Return data for the exact slope when interpreting self as a vector.

        Returns
        -------
        res : tuple : Data for the steepness of self w.r.t. the real-axis.

        Notes
        -----
        If possible, output is given by two relatively prime integers.

        This mirrors the role of 'theta' in the polar-coordinate
        expression of a complex number.

        The slope of point 0 is defined to be inf.

        See Also
        --------
        This function gives data for an exact form of the 'slope' method.

        """
        if self.real == 0:
            res = (float('inf'), 1) if (self.imag >= 0) else (-float('inf'), 1)
        else:
            gcd_ = self.gcdComponents
            res = (self.imag//gcd_, self.real//gcd_)
            
        return res

    @property
    def quadrant(self):
        """Return the quadrant (1 to 4) self is on.

        Notes
        -----
        If on an axis, the lowest quadrant is returned. For instance,
        0 lies on an axis (and the border of all quadrants)
        the minimal quadrant '0' is on is 1, which is returned.

        See Also
        --------
        This function provides the same information (in a different
        framework) as the 'unitPos' method.

        """
        if self.imag >= 0:
            res = 2 if (self.real < 0) else 1
        else:
            res = 4 if (self.real > 0) else 3
            
        return res

    @property
    def conjugate(self):
        """Return the conjugate.

        Returns
        -------
        - : GaussianInteger : The conjugate Gaussian integer.

        """
        return GaussianInteger(vec=(self.real, -self.imag))

    # Units Over Z[i]

    @property
    def default(self):
        """Return the equal-up-to-unit Gaussian integer with non-neg components.

        Notes
        -----
        This is a rotation by k*pi/2 degrees in the complex plane for some
        k between 0 and 3 (inclusive). The choice for non-neg components
        is 'natural' but arbitrary. For certain prime-related features this
        choice may not be optimal.
        
        """
        if (self.quadrant % 2) == 1:
            res = (abs(self.real), abs(self.imag))
        else:
            res = (abs(self.imag), abs(self.real))
            
        return GaussianInteger(vec=res)

    @property
    def unitKey(self):
        """Return a key associated with self up to unit."""
        selfUnit = self.default
        return selfUnit.real if selfUnit.isReal else selfUnit.vec

    @property
    def isUnit(self):
        """Return if self is a unit in Z[i].

        Notes
        -----
        The units in Z[i] are +/- 1 and +/- i.
        
        """
        return self.unitKey in {1, (0, 1)}

    @property
    def upToUnit(self):
        """Return all elements of Z[i] equivalent up to unit to self.

        Returns
        -------
        - : list : Rotations in C of self by k*pi/2 for k in {0, 1, 2, 3}.

        Notes
        -----
        All elements besides 0 have 4 unique elements equivalent up to unit.
        
        """
        if self.isAdditiveID:
            res = [(0, 0)]
        else:
            re, im = self.vec
            res = [(re, im), (-re, -im), (-im, re), (im, -re)]

        return [GaussianInteger(vec=rot_) for rot_ in res]

    def __i(self):
        """Return the tuple representing self*i."""
        return (-self.imag, self.real)

    def __neg(self):
        """Return the tuple representing -self."""
        return (-self.real, -self.imag)

    def __negi(self):
        """Return the tuple representing -self*i."""
        return (self.imag, -self.real)

    @property
    def unitPos(self):
        """Return the unit, u, so that self*u has non-negative components.

        Returns
        -------
        u : GaussianInteger : Unit in Z[i] s.t. self*u == self.default.

        Notes
        -----
        If self is 0, the chosen unit is defined to be 1.

        See Also
        --------
        This function provides the same information (in a different
        framework) as the 'quadrant' method.
        
        """
        rot_ = (1 - self.quadrant) % 4
        u = GaussianInteger(vec=powsOfi_tup(rot_))
        return u
        
    # Norm Related Features
    
    @property
    def norm(self):
        """Return the (number-theoretic) norm (sum of squared components).

        Notes
        -----
        For z = p + qi in Z[i] the norm is: p**2 + q**2. This is the
        square of what would be the Euclidean norm if viewing z as an
        element of C.

        In number theory (so often in Z[i]) integer-valued norms are
        preferred while in analysis (so often C) it is favorable to have
        norm(x) = x for large subsets (in the case of C, the norm function
        is the identity on the positive reals.)

        See Also
        --------
        The 'modulus' method is the square-root of the norm
        and coincides with the Euclidean norm over C.
        
        """
        return self.real**2 + self.imag**2

    @property
    def modulus(self):
        """Return the Euclidean norm (magnitude)."""
        return (self.norm) ** .5
    
    def __abs__(self):
        """Return the Euclidean norm (magnitude). Copy of 'modulus' method."""
        return self.modulus

    def __Euclidean_metric(self, zVec, sq=False):
        """Return either the Euclidean metric or its square.

        Parameters
        ----------
        zVec : tuple : Represents a point in Z[i].
        sq   : bool  : If returning square of metric.

        Returns
        -------
        res : int, float : Euclidean distance (or its square).
        
        """
        metSq = (self.real-zVec[0])**2 + (self.imag-zVec[1])**2
        res = int(metSq) if sq else (metSq)**.5
        return res

    def distSq(self, other):
        """"Return the square of the distance (via the Euclidean metric.)

        See Also
        --------
        See the 'proj_to_2tuple' function for attempt at projecting other
        to a form compatible with the '__Euclidean_metric' method.

        The 'distance' method is the square-root of this.
        
        """
        return self.__Euclidean_metric(proj_to_2tuple(other), sq=True)
            
    def distance(self, other):
        """"Return the distance (via the Euclidean metric.)

        See Also
        --------
        See the 'proj_to_2tuple' function for attempt at projecting other
        to a form compatible with the '__Euclidean_metric' method.

        The 'distSq' method is the square of this.
        
        """
        return self.__Euclidean_metric(proj_to_2tuple(other))

    # Comparing Points in Z[i]

    def __eq__(self, other):
        """Return if mathematically equivalent."""
        return self.vec == proj_to_2tuple(other)
        
    def __ne__(self, other):
        """Return if not mathematically equivalent."""
        return not(self.__eq__(other))

    def __gt__(self, other):
        """Return if self has the greater Z[i] norm.

        Parameters
        ----------
        other : * : Point (element) being compared.

        Notes
        -----
        * Anything ameible with 'proj_to_2tuple' can be used as a parameter.

        With respect to the Z[i] norm: 4 < -5.

        """
        try:
            otherNorm = other.Norm
        except:
            otherVec = proj_to_2tuple(other)
            otherNorm = euclidean_normSq(otherVec)
            
        return (self.norm > otherNorm)

    def __lt__(self, other):
        """Return if self has the lesser Z[i] norm.

        Parameters
        ----------
        other : * : Point (element) being compared.

        Notes
        -----
        * Anything ameible with 'proj_to_2tuple' can be used as a parameter.

        With respect to the Z[i] norm: 4 < -5.

        """
        try:
            otherNorm = other.Norm
        except:
            otherVec = proj_to_2tuple(other)
            otherNorm = euclidean_normSq(otherVec)
            
        return (self.norm < otherNorm)
    
    def __ge__(self, other):
        """Return if self has a greater (or equal to) Z[i] norm.

        Parameters
        ----------
        other : * : Point (element) being compared.

        Notes
        -----
        * Anything ameible with 'proj_to_2tuple' can be used as a parameter.

        With respect to the Z[i] norm: 4 < -5.

        Oftentimes, __ge__ == (__gt__ or __eq__) this is false in this setting.
        (__ge__, __gt__, __le__, __lt__) are statements about concentric
        circles while __eq__ is about the exact location on the plane.    

        """
        try:
            otherNorm = other.Norm
        except:
            otherVec = proj_to_2tuple(other)
            otherNorm = euclidean_normSq(otherVec)
            
        return (self.norm >= otherNorm)

    def __le__(self, other):
        """Return if self has a lesser (or equal to) Z[i] norm.

        Parameters
        ----------
        other : * : Point (element) being compared.

        Notes
        -----
        * Anything ameible with 'proj_to_2tuple' can be used as a parameter.

        With respect to the Z[i] norm: 4 < -5.

        Oftentimes, __ge__ == (__gt__ or __eq__) this is false in this setting.
        (__ge__, __gt__, __le__, __lt__) are statements about concentric
        circles while __eq__ is about the exact location on the plane.    

        """
        try:
            otherNorm = other.Norm
        except:
            otherVec = proj_to_2tuple(other)
            otherNorm = euclidean_normSq(otherVec)
            
        return (self.norm <= otherNorm)

    # Arithmetic

    def __componentWise(self, zVec, signs=(1,1)):
        """Return the weighted component wise sum.

        Parameters
        ----------
        zVec : tuple : Represents a point in Z[i].
        signs: tuple : Has form (+/-1,+/-1) represents sign of self and zVec.

        """
        s1, s2 = signs
        return tuple(map(lambda x,y: s1*x+s2*y, self.vec, zVec))

    def __add__(self, other):
        """"Return self + other.

        Parameters
        ----------
        other : * : Point (element) being compared.

        Notes
        -----
        * Anything ameible with 'proj_to_2tuple' can be used as a parameter.
        
        """
        otherVec = proj_to_2tuple(other)
        return GaussianInteger(vec=self.__componentWise(otherVec))

    def __radd__(self, other):
        """"Return other + self."""
        return self.__add__(other)
    
    def __iadd__(self, other):
        """"Return self + other."""
        return self.__add__(other)

    def __sub__(self, other):
        """"Return self - other.

        Parameters
        ----------
        other : * : Point (element) being compared.

        Notes
        -----
        * Anything ameible with 'proj_to_2tuple' can be used as a parameter.
        
        """
        otherVec = proj_to_2tuple(other)
        return GaussianInteger(vec=self.__componentWise(otherVec, signs=(1,-1)))

    def __rsub__(self, other):
        """"Return other - self."""
        otherVec = proj_to_2tuple(other)
        return GaussianInteger(vec=self.__componentWise(otherVec, signs=(-1,1)))
    
    def __isub__(self, other):
        """"Return self - other."""
        return self.__sub__(other)

    def __neg__(self):
        """Return -self."""
        return GaussianInteger(vec=self.__neg())

    def __componentMul(self, zVec=(1,0), scaler=1):
        """Return the scaled complex product as a tuple.

        Parameters
        ----------
        zVec  : tuple : Represents a point in Z[i].
        scaler: *     : Scaling factor, usually an int.

        Notes
        -----
        * Scaler requires a __mul__ operation compatible with integers.

        """
        if zVec == (1, 0):
            re_, im_ = self.vec
        else:
            re_ = self.real*zVec[0] - self.imag*zVec[1]
            im_ = self.real*zVec[1] + self.imag*zVec[0]
            
        return (scaler*re_, scaler*im_)
        
    def __mul__(self, other):
        """"Return self * other.

        Parameters
        ----------
        other : * : Point (element) being compared.

        Notes
        -----
        * Anything ameible with 'proj_to_2tuple' can be used as a parameter.
        
        """
        otherVec = proj_to_2tuple(other)
        return GaussianInteger(vec=self.__componentMul(zVec=otherVec))
        
    def __rmul__(self, other):
        """Return other*self."""
        return self.__mul__(other)

    def __imul__(self, other):
        """Return self*other."""
        return self.__mul__(other)

    def __sq(self):
        """Return the tuple associated with self**2."""
        return (self.real**2-self.imag**2, 2*self.real*self.imag)

    def __pow__(self, k):
        """Return self**k for k a non-negative integer."""
        return GaussianInteger(vec=implementPower(self.vec, int(k)))

    def __ipow__(self, k):
        """Return self**k for k a non-negative integer."""
        return self.__pow__(k)

    def scale(self, factor):
        """Return factor*self."""
        return GaussianInteger(vec=self.__componentMul(scaler=factor))

    @property
    def square(self):
        """Return self**2."""
        return GaussianInteger(vec=self.__componentMul(zVec=self.vec))

    @property
    def cube(self):
        """Return self**3."""
        re_ = self.real**3 - 3*self.imag**2*self.real
        im_ = -self.imag**3 + 3*self.real**2*self.imag
        return GaussianInteger(vec=(re_, im_))

    def __componentQuot(self, zVec=(1,0), scaler=1,
                        nearest=True, recip=False):
        """Return the scaled complex quotient as a tuple.

        Parameters
        ----------
        zVec   : tuple : Represents a point in Z[i].
        scaler : *     : Scaling factor, usually an int.
        nearest: bool  : If rounding to the nearest integer.
        recip  : bool  : If computing self/zVec, or its reciprocal self/zVec.
        
        Notes
        -----
        * Scaler requires a __mul__ operation compatible with integers.

        If rounding a value equidistant from multiple integers minimizing
        the norm is selected for. This can seem tricky if not careful about
        ambient space (and the norm being used).

        Example: (2.5, -2.5) -> (2, -2). If simply choosing the minimal
        value the result would be (2, -3) which has higher norm than
        (2, -2).

        """
        if zVec == (1,0):
            re_, im_ = self.vec
        else:
            norm_ = self.norm if recip else euclidean_normSq(zVec)
            sgn = -1 if recip else 1
            re_ = (self.real*zVec[0] + self.imag*zVec[1]) / norm_
            im_ = (self.imag*zVec[0] - self.real*zVec[1]) / (sgn*norm_)
            
        res = (scaler*re_, scaler*im_)
        return roundVec(res) if nearest else res

    def __floordiv__(self, other):
        """"Return self // other.

        Parameters
        ----------
        other : * : Point (element) being compared.

        Notes
        -----
        * Anything ameible with 'proj_to_2tuple' can be used as a parameter.
        
        """
        otherVec = proj_to_2tuple(other)
        return GaussianInteger(vec=self.__componentQuot(zVec=otherVec))

    def __rfloordiv__(self, other):
        """"Return other // self."""
        otherVec = proj_to_2tuple(other)
        return GaussianInteger(vec=self.__componentQuot(zVec=otherVec, recip=True))

    def __ifloordiv__(self, other):
        """"Return self // other."""
        return self.__floordiv__(other)

    def __truediv__(self, other):
        """"Return self / other.

        Parameters
        ----------
        other : * : Point (element) being compared.

        Notes
        -----
        * Anything ameible with 'proj_to_2tuple' can be used as a parameter.
        
        """
        otherVec = proj_to_2tuple(other)
        return self.__componentQuot(zVec=otherVec, nearest=False)

    def __rtruediv__(self, other):
        """"Return other / self."""
        otherVec = proj_to_2tuple(other)
        return self.__componentQuot(zVec=otherVec, nearest=False, recip=True)

    def __ifloordiv__(self, other):
        """"Return self / other."""
        return self.__truediv__(other)

    def __mod__(self, other):
        """"Return self % other.

        Parameters
        ----------
        other : * : Point (element) being compared.

        Notes
        -----
        * Anything ameible with 'proj_to_2tuple' can be used as a parameter.
        
        """
        k = self.__floordiv__(other)
        return self.__sub__(k*other)
            
    def __rmod__(self, other):
        """Return other % self."""
        k = self.__rfloordiv__(other)
        return (other - self.__mul__(k))
    
    def __imod__(self, other):
        """Return self % other."""
        return self.__mod__(other)

    # Euclidean Algorithm Over Z[i]

    def gcd(self, other):
        """"Return greatest common factor over Z[i].

        Parameters
        ----------
        other : * : Point (element) being compared.

        Notes
        -----
        * Anything ameible with 'proj_to_2tuple' can be used as a parameter.
        
        """
        zVec = proj_to_2tuple(other)
        z, w = GaussianInteger(vec=zVec), GaussianInteger(vec=self.vec)
        while w.norm > 0:
            z, w = w, z%w

        return z
