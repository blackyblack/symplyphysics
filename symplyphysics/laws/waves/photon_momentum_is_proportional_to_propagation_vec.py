from sympy import (Eq, solve, pi)
from sympy.physics.units import planck as planck_constant
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Photon is the elementary part of any electromagnetical radiation which has no mass and always moves with speed of light.
## In despite of no mass, it carries some momentum. The amount of this momentum depends only on the propagation vector of the photon.
## Law is: p = (h/2pi) * k, where
## p is momentum of photon,
## h is Planck constant,
## k is number of wave cycles per length.

photon_momentum = Symbol("photon_momentum", units.momentum)
wavenumber = Symbol("wavenumber", 1/units.length)

law = Eq(photon_momentum, (planck_constant/(2*pi)) * wavenumber)


def print_law() -> str:
    return print_expression(law)


@validate_input(wavenumber_=wavenumber)
@validate_output(photon_momentum)
def calculate_momentum(wavenumber_: Quantity) -> Quantity:
    result_momentum_expr = solve(law, photon_momentum, dict=True)[0][photon_momentum]
    result_expr = result_momentum_expr.subs(
        {wavenumber: wavenumber_})
    return Quantity(result_expr)
