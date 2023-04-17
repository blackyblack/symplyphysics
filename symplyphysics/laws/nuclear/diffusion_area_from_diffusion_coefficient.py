from sympy import Expr
from symplyphysics import (
    Eq, pretty, solve, units, expr_to_quantity
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

# Description
## The physical meaning of the diffusion length can be seen by calculating the mean square distance that
## a neutron travels in the one direction from the plane source to its absorption point.

## Law: L^2 = D / Σa
## Where:
## Σa - macroscopic absorption cross-section of the fuel.
##   See [macroscopic transport cross-section](./macroscopic_transport_cross_section.py) implementation.
## D - diffusion coefficient.
##   See [diffusion coefficient](./neutron_diffusion_coefficient_from_scattering_cross_section.py) implementation.
## L^2 - diffusion area.
## L - diffusion length.

diffusion_coefficient = Symbol("diffusion_coefficient", units.length)
macroscopic_absorption_cross_section = Symbol("macroscopic_absorption_cross_section", 1 / units.length)
diffusion_area = Symbol("diffusion_area", units.length**2)

law = Eq(diffusion_area, diffusion_coefficient / macroscopic_absorption_cross_section)

def print(expr: Expr) -> str:
    symbols = [diffusion_coefficient, macroscopic_absorption_cross_section, diffusion_area]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(diffusion_coefficient_=diffusion_coefficient, macroscopic_absorption_cross_section_=macroscopic_absorption_cross_section)
@validate_output_symbol(diffusion_area)
def calculate_diffusion_area(
    diffusion_coefficient_: Quantity,
    macroscopic_absorption_cross_section_: Quantity) -> Quantity:
    result_diffusion_expr = solve(law, diffusion_area, dict=True)[0][diffusion_area]
    result_expr = result_diffusion_expr.subs({
        diffusion_coefficient: diffusion_coefficient_,
        macroscopic_absorption_cross_section: macroscopic_absorption_cross_section_})
    return expr_to_quantity(result_expr)
