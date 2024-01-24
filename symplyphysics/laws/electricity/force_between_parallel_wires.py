from sympy import (Eq, solve, pi)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless
)
from sympy.physics.units import magnetic_constant

# Description
## Two parallel wires through which current flows interact with each other.
## The intensity of their interaction depends on the magnitude of the current,
## the distance from the wires, as well as on the material and length of the wire.

## Law is: F = mu0 * mu * I1 * I2 * l / (2 * pi * r), where
## F - force of interaction of wires,
## mu0 - magnetic constant,
## mu - relative permeability of medium,
## I1 - current in the first wire,
## I2 - current in the second wire,
## l - length of wires,
## r - distance.

force = Symbol("force", units.force)

relative_permeability = Symbol("relative_permeability", dimensionless)
first_wire_current = Symbol("first_wire_current", units.current)
second_wire_current = Symbol("second_wire_current", units.current)
length = Symbol("length", units.length)
distance = Symbol("distance", units.length)

law = Eq(force,
         magnetic_constant * relative_permeability * first_wire_current * second_wire_current * length
         / (2 * pi * distance))


def print_law() -> str:
    return print_expression(law)


@validate_input(relative_permeability_=relative_permeability,
    first_wire_current_=first_wire_current,
    second_wire_current_=second_wire_current,
    length_=length,
    distance_=distance)
@validate_output(force)
def calculate_force(relative_permeability_: float, first_wire_current_: Quantity,
    second_wire_current_: Quantity, length_: Quantity,
    distance_: Quantity) -> Quantity:
    result_expr = solve(law, force, dict=True)[0][force]
    result_expr = result_expr.subs({
        relative_permeability: relative_permeability_,
        first_wire_current: first_wire_current_,
        second_wire_current: second_wire_current_,
        length: length_,
        distance: distance_
    })
    return Quantity(result_expr)
