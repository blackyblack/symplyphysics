from sympy import (Eq, solve)
from symplyphysics import (clone_symbol, symbols, units, Quantity, Symbol, print_expression,
    validate_input, validate_output)

# Description
## Coulomb's law states that the force F between two point charges, q1 and q2, in a vacuum is proportional to their product
## and inversely proportional to the square of their distance apart r.

# Law: F = k_e * ((q_1 * q_2)/r_{12}^2), where
## F is force between two electrically charged particles at rest (in vacuum),
## k_e is a Coulomb constant,
## q_1 and q_2 are the quantities of each charge,
## r is the distance between the charges.

# Note
## The symbols Q1 and Q2 in the Coulomb's law equation represent the quantities of charge on the two interacting objects.
## When using the "+" and "-" signs in the calculation of force,
## the result will be that a "-" value for force is a sign of an attractive force
## and a "+" value for force signifies a repulsive force. Mathematically, the force value
## would be found to be positive when Q1 and Q2 are of like charge - either both "+" or both "-".
## And the force value would be found to be negative when Q1 and Q2 are of opposite charge - one is "+"
## and the other is "-". This is consistent with the concept that oppositely charged objects
## have an attractive interaction and like charged objects have a repulsive interaction.

electrostatic_force = clone_symbol(symbols.dynamics.force)
first_charge = Symbol("first_charge", units.charge)
second_charge = Symbol("second_charge", units.charge)
distance = Symbol("distance", units.length)

law = Eq(electrostatic_force,
    (units.coulomb_constant * first_charge * second_charge) / (distance**2))


def print_law() -> str:
    return print_expression(law)


@validate_input(first_charge_=first_charge, second_charge_=second_charge, distance_=distance)
@validate_output(electrostatic_force)
def calculate_force(first_charge_: Quantity, second_charge_: Quantity,
    distance_: Quantity) -> Quantity:
    solved = solve(law, electrostatic_force, dict=True)[0][electrostatic_force]
    result_expr = solved.subs({
        first_charge: first_charge_,
        second_charge: second_charge_,
        distance: distance_,
    })
    return Quantity(result_expr)
