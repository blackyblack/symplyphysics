from sympy import (Eq, solve)
from sympy.physics.units import magnetic_constant
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input, validate_output, Dimensionless)

# Description
## The basic characteristic of a coil is its inductance - the ability of the coil to accumulate energy as magnetic field.
## Law: L = mu * mu_0 * N**2 * S / d, where
## L is the inductance of a coil,
## mu is magnetic permeability of core,
## mu_0 is magnetic constant (magnetic permeability of vacuum),
## N is number of turns,
## S is area of each turn,
## d is coil length.

coil_inductance = Symbol("coil_inductance", units.inductance)
magnetic_permeability = Symbol("magnetic_permeability", Dimensionless)
number_of_turns = Symbol("number_of_turns", Dimensionless)
turn_area = Symbol("turn_area", units.area)
coil_length = Symbol("coil_length", units.length)

law = Eq(coil_inductance,
    magnetic_constant * magnetic_permeability * number_of_turns**2 * turn_area / coil_length)


def print_law() -> str:
    return print_expression(law)


@validate_input(turn_area_=turn_area, coil_length_=coil_length)
@validate_output(coil_inductance)
def calculate_inductance(magnetic_permeability_: float, number_of_turns_: float,
    turn_area_: Quantity, coil_length_: Quantity) -> Quantity:
    result_inductance_expr = solve(law, coil_inductance, dict=True)[0][coil_inductance]
    result_expr = result_inductance_expr.subs({
        magnetic_permeability: magnetic_permeability_,
        number_of_turns: number_of_turns_,
        turn_area: turn_area_,
        coil_length: coil_length_
    })
    return expr_to_quantity(result_expr)
