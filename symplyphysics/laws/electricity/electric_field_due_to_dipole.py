from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The value of the electric field set up by a dipole at a distant point on the dipole axis 
## (which runs through both particles) is proportional to the inverse cube of the distance 
## to the dipole and the value of the dipole moment.

# Law:
## E = 2 * k_e * p / z^3
## E - value of electric field of dipole
## k_e - Coulomb's constant
## p = ql - value of the dipole moment
## z - distance to dipole on the dipole axis

# Condition
## z/l << 1 - the point of measuring the electric field should be far enough from the dipole itself

electric_field = Symbol("electric_field", units.force / units.charge)
dipole_moment = Symbol("dipole_moment", units.charge * units.length)
distance_to_dipole = Symbol("distance_to_dipole", units.length)

law = Eq(electric_field, 2 * units.coulomb_constant * dipole_moment / distance_to_dipole**3)


def print_law() -> str:
    return print_expression(law)


@validate_input(dipole_moment_=dipole_moment, distance_to_dipole_=distance_to_dipole)
@validate_output(electric_field)
def calculate_electric_field(dipole_moment_: Quantity, distance_to_dipole_: Quantity) -> Quantity:
    result = solve(law, electric_field)[0]
    result_field = result.subs({
        dipole_moment: dipole_moment_,
        distance_to_dipole: distance_to_dipole_,
    })
    return Quantity(result_field)
