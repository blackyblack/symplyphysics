from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input, validate_output)

# Description
## is an experimental law of physics that calculates the amount of force between two electrically charged particles at rest.
## This electric force is conventionally called electrostatic force or Coulomb force

# Law: F = k_e * ((q_1 * q_2)/r_{12}^2), where
## F is force between two electrically charged particles at rest (in vacuum),
## k_e is a Coulomb constant,
## q_1 and q_2 are the quantities of each charge,
## r is the distance between the charges.

force_btw_charges = Symbol("force", units.force)
first_charge = Symbol("first_charge", units.charge)
second_charge = Symbol("second_charge", units.charge)
distance_btw_charges = Symbol("distance_btw_charges", units.length)

law = Eq(force_btw_charges, (units.coulomb_constant * first_charge * second_charge)/(distance_btw_charges**2))

def print_law() -> str:
    return print_expression(law)

@validate_input(first_charge_ = first_charge, second_charge_ = second_charge, distance_btw_charges_ = distance_btw_charges)
@validate_output(force_btw_charges)
def calculate_force_btw_two_charged_particles(first_charge_: Quantity, second_charge_: Quantity, distance_btw_charges_: Quantity) -> Quantity:
    solved = solve(law, force_btw_charges, dict=True)[0][force_btw_charges]
    result_expr = solved.subs({
        first_charge: first_charge_,
        second_charge: second_charge_,
        distance_btw_charges: distance_btw_charges_,
    })
    return Quantity(result_expr)
