from sympy import (Eq, solve)
from sympy.physics.units import planck as planck_constant
from sympy.physics.units import speed_of_light
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## Photon is the elementary part of any electromagnetical radiation which has no mass and always moves with speed of light.
## In despite of no mass, it carries some momentum. The amount of this momentum depends only on the frequency of the photon.
## Law is: p = h * nu / c, where
## p is momentum of photon,
## h is Planck constant,
## nu is frequency of photon,
## c is speed of light.

photon_momentum = Symbol("photon_momentum", units.momentum)
photon_frequency = Symbol("frequency", units.frequency)

law = Eq(photon_momentum, planck_constant * photon_frequency / speed_of_light)


def print_law() -> str:
    return print_expression(law)


@validate_input(photon_frequency_=photon_frequency)
@validate_output(photon_momentum)
def calculate_momentum(photon_frequency_: Quantity) -> Quantity:
    result_momentum_expr = solve(law, photon_momentum, dict=True)[0][photon_momentum]
    result_expr = result_momentum_expr.subs({photon_frequency: photon_frequency_})
    return Quantity(result_expr)
