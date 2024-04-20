from sympy import (Eq, solve, sqrt)
from sympy.physics.units import speed_of_light
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input, validate_output,
    dimensionless,)

# Description
## The ring resonator is a microstrip line in the shape of a circle. The effective permeability of a transmission line is such
## a permeability that a wave would propagate in an unlimited medium in the same way as it propagates in a transmission line.

## Law is: f = m * c / (l * eps), where
## f - the resonant frequency of the ring resonator,
## m - order interference,
## l - the length of the ring resonator,
## eps - effective permittivity of the resonator.

frequency = Symbol("frequency", units.frequency)

ring_length = Symbol("ring_length", units.length)
order_interference = Symbol("order_interference", dimensionless)
permittivity = Symbol("permittivity", dimensionless)

law = Eq(frequency, order_interference * speed_of_light / (ring_length * sqrt(permittivity)))


def print_law() -> str:
    return print_expression(law)


@validate_input(ring_length_=ring_length,
    order_interference_=order_interference,
    permittivity_=permittivity)
@validate_output(frequency)
def calculate_frequency(ring_length_: Quantity, order_interference_: int, permittivity_: float) -> Quantity:
    result_expr = solve(law, frequency, dict=True)[0][frequency]
    result_expr = result_expr.subs({
        ring_length: ring_length_,
        order_interference: order_interference_,
        permittivity: permittivity_,
    })
    return Quantity(result_expr)
