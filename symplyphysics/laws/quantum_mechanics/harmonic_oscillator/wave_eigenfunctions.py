from sympy import Eq, Rational, sqrt, exp, factorial, pi, pprint
from sympy.functions.special.polynomials import hermite
from sympy.physics.units import hbar
from symplyphysics import (
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    units,
    dimensionless,
    angle_type,
    symbols,
    clone_symbol,
    convert_to_float,
)

# Description
## ...

# Law: psi_n(x) = 1 / sqrt(2**n * factorial(n))
#               * ((m * w) / (pi * hbar))**(1 / 4)
#               * exp(-1 * (m * w) / (2 * hbar) * x**2)
#               * H_n(sqrt((m * w) / hbar) * x)
## psi_n - n-th wave eigenfunction
## n - mode number
## m - mass of oscillating quantum particle
## w - angular frequency of quantum oscillator
## hbar - reduced Planck constant
## H_n - n-th (physicists') Hermite polynomial

wave_function = Symbol("wave_function", 1 / sqrt(units.length))
mode_number = Symbol("mode_number", dimensionless, integer=True, nonnegative=True)
oscillator_mass = clone_symbol(symbols.basic.mass, "oscillator_mass")
angular_frequency = Symbol("angular_frequency", angle_type / units.time)
position = Symbol("position", units.length)

law = Eq(
    wave_function,
    (1 / (sqrt(2**mode_number * factorial(mode_number))))
    * ((oscillator_mass * angular_frequency) / (pi * hbar))**Rational(1, 4)
    * exp(-1 * (oscillator_mass * angular_frequency) / (2 * hbar) * position**2)
    * hermite(mode_number, sqrt(oscillator_mass * angular_frequency / hbar) * position)
)


@validate_input(
    mode_number_=mode_number,
    oscillator_mass_=oscillator_mass,
    angular_frequency_=angular_frequency,
    position_=position,
)
@validate_output(wave_function)
def calculate_wave_function(
    mode_number_: int,
    oscillator_mass_: Quantity,
    angular_frequency_: Quantity,
    position_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        mode_number: mode_number_,
        oscillator_mass: oscillator_mass_,
        angular_frequency: angular_frequency_,
        position: position_,
    })
    return Quantity(result)
