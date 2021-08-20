# Molar gas constant (also known as the gas constant, universal gas constant, ideal gas constant, or R)

## Unit is J*K^-1*mol^-1

from math import isclose

from .boltzmann import value as k
from .avogadro import value as Na

value = 8.3145
assert isclose(value, Na * k, rel_tol = 0.0001)