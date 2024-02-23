from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The bar stiffness coefficient depends on the Young's module, cross-sectional area and length.
## Young's module is a tabular value that differs for each material.

## Law is: k = E * S / l, where
## k - coefficient of stiffness,
## E - Young's module,
## S - area,
## l - length.

coefficient_of_stiffness = Symbol("coefficient_of_stiffness", units.force / units.length)

module_of_young = Symbol("module_of_young", units.pressure)
area = Symbol("area", units.area)
length = Symbol("length", units.length)

law = Eq(coefficient_of_stiffness, module_of_young * area / length)


def print_law() -> str:
    return print_expression(law)


@validate_input(module_of_young_=module_of_young,
    area_=area,
    length_=length)
@validate_output(coefficient_of_stiffness)
def calculate_coefficient_of_stiffness(module_of_young_: Quantity, area_: Quantity,
    length_: Quantity) -> Quantity:
    result_expr = solve(law, coefficient_of_stiffness, dict=True)[0][coefficient_of_stiffness]
    result_expr = result_expr.subs({
        module_of_young: module_of_young_,
        area: area_,
        length: length_
    })
    return Quantity(result_expr)
