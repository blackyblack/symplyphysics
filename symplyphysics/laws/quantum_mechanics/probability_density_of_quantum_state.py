from sympy import Eq, sqrt
from symplyphysics import (
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    units,
)

# Description
## Probability density of a quantum system is the likelihood that the system will be in a particular
## position at a particular time.

# Law: rho(x, t) = Psi(x, t) * conj(Psi(x, t)) = abs(Psi(x, t))**2
## rho - probability density
## Psi - wave function
## x - position
## t - time
## abs - absolute value
## conj(z) - complex conjugate of `z`

# Links: Wikipedia <https://en.wikipedia.org/wiki/Wave_function#Position-space_wave_functions>

probability_density = Symbol("probability_density", 1 / units.length)
wave_function = Symbol("wave_function", 1 / sqrt(units.length))

law = Eq(probability_density, abs(wave_function)**2)


@validate_input(wave_function_value_=wave_function)
@validate_output(probability_density)
def calculate_probability_density(wave_function_value_: Quantity) -> Quantity:
    result = law.rhs.subs(wave_function, wave_function_value_)
    return Quantity(result)
