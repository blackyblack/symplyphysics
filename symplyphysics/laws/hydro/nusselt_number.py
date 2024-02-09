from sympy import Eq, solve, S
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
                           validate_output, dimensionless, convert_to)


# Description
# Nusselt number is the ratio of convective to conductive heat transfer at
# a boundary in a fluid. There is a characteristic length in
# the formula. The characteristic length is the dimension
# that defines the length scale of a physical system. A characteristic length
# is usually the volume of a system divided by its surface: L = V / A,
# where V is the volume of the body, and A is the cross-sectional area.
# For example, it is used to calculate flow through circular and non-circular
# tubes in order to examine flow conditions. D = 4 * A / p, where
# D is characteristic diameter, A is the cross-sectional are, p is wetted perimeter.
# Law: Nu = h * L / k, where
# h is the heat transfer coefficient,
# L is the characteristic length,
# k is the thermal conductivity,
# Nu is Nusselt number.


heat_transfer_coefficient = Symbol(
    "heat_transfer_coefficient", units.power / (units.area * units.temperature)
)
characteristic_length = Symbol("characteristic_length", units.length)
thermal_conductivity = Symbol(
    "thermal_conductivity", units.power / (units.length * units.temperature)
)

nusselt_number = Symbol("nusselt_number", dimensionless)

law = Eq(nusselt_number, heat_transfer_coefficient *
         characteristic_length / thermal_conductivity)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    heat_transfer_coefficient_=heat_transfer_coefficient,
    characteristic_length_=characteristic_length,
    thermal_conductivity_=thermal_conductivity,
)
@validate_output(nusselt_number)
def calculate_nusselt_number(
    heat_transfer_coefficient_: Quantity,
    characteristic_length_: Quantity,
    thermal_conductivity_: Quantity
) -> float:
    result_expr = solve(law, nusselt_number, dict=True)[0][nusselt_number]
    result_applied = result_expr.subs({
        heat_transfer_coefficient: heat_transfer_coefficient_,
        characteristic_length: characteristic_length_,
        thermal_conductivity: thermal_conductivity_
    })
    result = Quantity(result_applied)
    return float(convert_to(result, S.One).evalf())
