from sympy import Eq
from symplyphysics import (
    Quantity,
    Symbol,
    Function,
    symbols,
    units,
    validate_input,
    validate_output,
)

# Description
## Differental Joule-Thompson effect ...

# Law: (dT/dp)_H = (T * (dV/dT)_p - V)/C_p
## T - temperature
## p - pressure
## V - volume
## H - enthalpy
## C_p - isobaric heat capacity
