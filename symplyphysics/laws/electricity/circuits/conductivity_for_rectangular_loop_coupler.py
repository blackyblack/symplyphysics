from sympy import Eq, solve, Matrix, sqrt
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
)

## Description
## The rectangular loop coupler consists of four sections. The conductivity of each section can be calculated by calculating
## the conductivity of the transmission line to which the coupler is connected and the power ratio at the outputs.

## Law is: Matrix([Y1, Y2, Y3, Y4]) = Matrix([Y0 / sqrt(k), Y0 * sqrt((k + 1) / k), Y0 * sqrt((k + 1) / k), Y0 / sqrt(k)]), where
## Y1 - first conductivity,
## Y2 - second conductivity,
## Y3 - third conductivity,
## Y4 - fourth conductivity,
## Y0 - transmission line conductivity,
## k - ratio coefficient of the power at the outputs of the coupler.

first_conductivity = Symbol("first_conductivity", units.conductance)
second_conductivity = Symbol("second_conductivity", units.conductance)
third_conductivity = Symbol("third_conductivity", units.conductance)
fourth_conductivity = Symbol("fourth_conductivity", units.conductance)

transmission_line_conductivity = Symbol("transmission_line_conductivity", units.conductance)
ratio_of_power = Symbol("ratio_of_power", dimensionless)


law = Eq(Matrix([first_conductivity, second_conductivity, third_conductivity, fourth_conductivity]),
         Matrix([transmission_line_conductivity / sqrt(ratio_of_power), transmission_line_conductivity * sqrt((ratio_of_power + 1) / ratio_of_power),
                 transmission_line_conductivity * sqrt((ratio_of_power + 1) / ratio_of_power), transmission_line_conductivity / sqrt(ratio_of_power)]))

def print_law() -> str:
    return print_expression(law)


@validate_input(transmission_line_conductivity_=transmission_line_conductivity, ratio_of_power_=ratio_of_power)
@validate_output(units.conductance)
def calculate_conductivities(transmission_line_conductivity_: Quantity, ratio_of_power_: float) -> tuple[Quantity, Quantity, Quantity, Quantity]:
    result = solve(law, [first_conductivity, second_conductivity, third_conductivity, fourth_conductivity], dict=True)[0]
    result_Y1 = result[first_conductivity]
    result_Y2 = result[second_conductivity]
    result_Y3 = result[third_conductivity]
    result_Y4 = result[fourth_conductivity]
    substitutions = {
        transmission_line_conductivity: transmission_line_conductivity_,
        ratio_of_power: ratio_of_power_,
    }
    result_Y1 = Quantity(result_Y1.subs(substitutions))
    result_Y2 = Quantity(result_Y2.subs(substitutions))
    result_Y3 = Quantity(result_Y3.subs(substitutions))
    result_Y4 = Quantity(result_Y4.subs(substitutions))
    return (result_Y1, result_Y2, result_Y3, result_Y4)
