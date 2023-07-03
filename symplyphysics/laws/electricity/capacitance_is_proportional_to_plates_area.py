from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input, validate_output, Dimensionless)
from sympy.physics.units import electric_constant

# Description
## The basic characteristic of a capacitor is its capacitance - the ability of the capacitor to accumulate an electric charge.
## Law: C = epsilon * epsilon_0 * S / d, where
## C is capasitance of capacitor,
## epsilon is dielectric permeability of capacitor insulator,
## epsilon_0 is dielectric constant,
## S is area of each plate,
## d is plate-to-plate clearance.

# Conditions
## The whole electric field is condensed in capacitor between plates, no losses.

capacitor_capacitance = Symbol("capacitor_capacitance", units.capacitance)
dielectric_permeability = Symbol("dielectric_permeability", Dimensionless)
plate_area = Symbol("plate_area", units.area)
clearance = Symbol("clearance", units.length)

law = Eq(capacitor_capacitance, electric_constant * dielectric_permeability * plate_area / clearance)


def print_law() -> str:
    return print_expression(law)


@validate_input(plate_area_=plate_area, clearance_=clearance)
@validate_output(capacitor_capacitance)
def calculate_capacitance(dielectric_permeability_: float, plate_area_: Quantity, clearance_: Quantity) -> Quantity:
    result_capacitance_expr = solve(law, capacitor_capacitance, dict=True)[0][capacitor_capacitance]
    result_expr = result_capacitance_expr.subs({
        dielectric_permeability: dielectric_permeability_,
        plate_area: plate_area_,
        clearance: clearance_
    })
    return expr_to_quantity(result_expr)
