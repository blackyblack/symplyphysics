from sympy import (Eq, solve)
from sympy.physics.units import electric_constant
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input, validate_output, Dimensionless)

# Description
## The basic characteristic of a capacitor is its capacitance - the ability of the capacitor to accumulate an electric charge.
## Law: C = epsilon * epsilon_0 * S / d, where
## C is capacitance of capacitor,
## epsilon is dielectric permeability of capacitor insulator,
## epsilon_0 is dielectric constant,
## S is area of each plate,
## d is distance between plates.

# Conditions
## - The whole electric field is condensed in space between plates, no losses.

capacitor_capacitance = Symbol("capacitor_capacitance", units.capacitance)
dielectric_permeability = Symbol("dielectric_permeability", Dimensionless)
plate_area = Symbol("plate_area", units.area)
distance_between_plates = Symbol("distance_between_plates", units.length)

law = Eq(capacitor_capacitance, electric_constant * dielectric_permeability * plate_area / distance_between_plates)


def print_law() -> str:
    return print_expression(law)


@validate_input(plate_area_=plate_area, distance_between_plates_=distance_between_plates)
@validate_output(capacitor_capacitance)
def calculate_capacitance(dielectric_permeability_: float, plate_area_: Quantity, distance_between_plates_: Quantity) -> Quantity:
    result_capacitance_expr = solve(law, capacitor_capacitance, dict=True)[0][capacitor_capacitance]
    result_expr = result_capacitance_expr.subs({
        dielectric_permeability: dielectric_permeability_,
        plate_area: plate_area_,
        distance_between_plates: distance_between_plates_
    })
    return expr_to_quantity(result_expr)
