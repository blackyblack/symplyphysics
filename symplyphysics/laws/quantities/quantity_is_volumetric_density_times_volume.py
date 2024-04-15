from sympy import Eq, Symbol as SymSymbol
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
)
from symplyphysics.core.dimensions import assert_equivalent_dimension

# Description
## An extensive quantity of interest can be obtained by multiplying the corresponding volumetric
## density by volume

# Law: X = rho_X * V
## X - extensive quantity
## rho_X - (intensive) volumetric density
## V - volume

extensive_quantity = SymSymbol("extensive_quantity")
volumetric_density = SymSymbol("volumetric_density")
volume = Symbol("volume", units.volume)

law = Eq(extensive_quantity, volumetric_density * volume)


def print_law() -> str:
    return print_expression(law)


@validate_input(volume_=volume)
def calculate_extensive_quantity(
    volumetric_density_: Quantity,
    volume_: Quantity,
) -> Quantity:
    result_expr = law.rhs.subs({
        volumetric_density: volumetric_density_,
        volume: volume_,
    })
    result = Quantity(result_expr)

    assert_equivalent_dimension(
        result,
        "return",
        "calculate_extensive_quantity",
        volumetric_density_.dimension * units.volume,
    )

    return result
