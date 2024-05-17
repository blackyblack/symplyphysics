from sympy import (Eq, solve, sqrt)
from sympy.physics.units import speed_of_light
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
)

# Description
## The ring resonator is a microstrip line in the shape of a circle.
## When a wave propagates along a microstrip line, part of the field goes out, since the microstrip line does
## not have metal borders on all sides, unlike, for example, rectangular waveguides. Then imagine an environment
## in which the field will have the same magnitude as the field of a microstrip line. The permittivity of such a
## medium will be called the effective permittivity of the line.
## A wave traveling through an ring resonator acquires a phase shift and interacts with a wave incident on the
## resonator. If the phase shift is a multiple of 2*pi*m, then these waves add up in phase. m is the interference order.

## Law is: f = m * c / (l * eps), where
## f - the resonant frequency of the ring resonator,
## m - order interference,
## c - speed of light,
## l - the length of circumference of the ring resonator,
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
def calculate_frequency(ring_length_: Quantity, order_interference_: int,
    permittivity_: float) -> Quantity:
    result_expr = solve(law, frequency, dict=True)[0][frequency]
    result_expr = result_expr.subs({
        ring_length: ring_length_,
        order_interference: order_interference_,
        permittivity: permittivity_,
    })
    return Quantity(result_expr)
