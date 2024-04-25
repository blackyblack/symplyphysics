from sympy import Eq, Rational
from symplyphysics import (
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    units,
    dimensionless,
    angle_type,
)

# Description
## As opposed to the classical harmonic oscillator, the energy levels of a quantum harmonic oscillator
## are quantized, meaning that its energy take a value out of a discrete range. These energy levels
## are equidistant, i.e. the difference between successive energy levels is the same for all levels.

# Law: E_n = (n + 1/2) * hbar * omega
## E_n - energy of n-th level
## n - quantum number of oscillator, which is any non-negative integer (0, 1, 2, ...)
## hbar - reduced Planck constant
## omega - angular frequency of oscillator

# Notes
## - This means that the energy of a quantum oscillator cannot be zero and the lowest it can be
##   is the zero-point energy `E_0 = hbar * omega / 2`.

energy_level = Symbol("energy_level", units.energy)
quantum_number = Symbol("quantum_number", dimensionless, integer=True, nonnegative=True)
angular_frequency = Symbol("angular_frequency", angle_type / units.time)

law = Eq(
    energy_level,
    (quantum_number + Rational(1, 2)) * units.hbar * angular_frequency,
)


@validate_input(
    quantum_number_=quantum_number,
    angular_frequency_=angular_frequency,
)
@validate_output(energy_level)
def calculate_energy_level(
    quantum_number_: int,
    angular_frequency_: int,
) -> Quantity:
    result = law.rhs.subs({
        quantum_number: quantum_number_,
        angular_frequency: angular_frequency_,
    })
    return Quantity(result)
