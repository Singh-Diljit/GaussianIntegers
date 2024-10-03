# GaussianIntegers

## Project Overview

This Python project implements a class for **Gaussian Integers** and includes a function to compute the **prime decomposition** and **all divisors** of Gaussian integers.

## Gaussian Integers

Gaussian integers are complex numbers where both the real and imaginary parts are integers. They form a **Euclidean domain**, which means that unique factorization into Gaussian primes is possible, similar to the fundamental theorem of arithmetic for the integers. The set of Gaussian integers is denoted $Z[i]$.

## Class Structure

### `GaussianInteger`

The core of the project is the `GaussianInteger` class, which represents and manipulates Gaussian integers.

#### Attributes:
- `real` (int): The real part, $a$, of the Gaussian integer.
- `imag` (int): The imaginary part, $b$, of the Gaussian integer.
- `vec` (tuple): A tuple of the form ($a$, $b$).

#### Select Methods:

1. The class features basic classification, representation, and accessing features such as:
	1. `isReal(self)`, and `isPurelyImag(self)`
	2. `isAdditiveID(self)`, and  `isMultiplicativeID(self)`
	3. `__bool__(self)`, `__hash__(self)`, `__repr__(self)`, and `__str__(self)`
	
2. Standard complex number features such as:
	1. `slope(self)`, and `slopeExact(self)`
	2. `conjugate(self)`, `norm(self)`, and `modulus(self)`
	
3. Various arithmetic related features (implemented to support cross class operations).
4.  Greatest common factor method: `gcd(self, other)` for computing the GCD over Z[i].
   
## Prime Decomposition

### Function: `primeDecomp(z: GaussianInteger)`

This function decomposes a Gaussian integer, $z$, into its prime factors over the Gaussian integers. 

For primes, $p \in Z$ a prime in $Z[i]$ can be constructed as follows:
1. If $p = 2$ then $1+i$ is a prime in $Z[i]$, namely you have: $2 = -i(1+i)^2$. This is the only prime to **ramify** in $Z[i]$. This can be proved by looking at the discriminant of $Z[i]$ and noting every prime that ramifies must divide that number (which is $-4$).

2. If $p = 3 \pmod 4$ then $p$ is also a prime in $Z[i]$. This is proved by considering norms and taking Fermat's theorem on sums of two squares into account. These primes are called **inert**.

3. If $p=1 \pmod 4$ then by Fermat's theorem on sums of two squares, there exist (unique) $a, b\in Z$ such that $p = a^2 + b^2$ and $p$ decomposed into two primes over $Z[i]$, namely $a+bi$ and $a-bi$. Again this can be proved by considering the multiplicative nature of norms. These primes are said to be **decomposed primes**.

The algorithm is designed to handle the factorization of numbers whose $Z$ - associated prime is under $10^8$. This can easily be extended by either banking more primes in the `primeDataZi` folder or including more sophisticated prime decomposition methods over $Z$. See `smallPrimes` project on my github to bypass this limitation.

#### Example Usage:

```python
# Example: Decomposing a Gaussian integer

el = GaussianInteger(26, 18)
print(primeDecomp(el))
>>> {(1, 1): 3, (2, 1): 3}

el = GaussianInteger(5, 5)
print([str(x) for x in divisorsZi(el)])
>>> ['2 - i', '5 + 5i', '2 + i', '3 + i', '1 + i', '5', '1', '1 + 3i']

print([str(x) for x in divisorsZi(el, proper=True)])
>>> ['2 - i', '2 + i', '3 + i', '1 + i', '5', '1 + 3i']

print([str(x) for x in divisorsZi(el, proper=True, sort=True)])
>>> ['1 + i', '2 - i', '2 + i', '3 + i', '1 + 3i', '5']
```

## Installation

This project requires **Python 3.x**. No additional dependencies are necessary.
