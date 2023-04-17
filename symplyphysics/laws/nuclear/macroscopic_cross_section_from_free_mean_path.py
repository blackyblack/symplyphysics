from sympy import Expr
from symplyphysics import (
    Eq, pretty, solve, units, expr_to_quantity
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

# Description
## Macroscopic cross-section - represents the effective target area of all of the nuclei contained
## in the volume of the material (such as fuel pellet). It is the probability of neutron-nucleus interaction per
## centimeter of neutron travel.

## Macroscopic cross-section: Σ = 1 / λ
## Where:
## λ (mean free path) is equal to the average value of x, the distance traveled by a neutron without any
## interaction, over the interaction probability distribution.
## Σ is the macroscopic cross-section.

mean_free_path = Symbol("mean_free_path", units.length)
macroscopic_cross_section = Symbol("macroscopic_cross_section", 1 / units.length)

law = Eq(macroscopic_cross_section, 1 / mean_free_path)

def print(expr: Expr) -> str:
    symbols = [mean_free_path, macroscopic_cross_section]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(mean_free_path_=mean_free_path)
@validate_output_symbol(macroscopic_cross_section)
def calculate_cross_section(mean_free_path_: Quantity) -> Quantity:
    result_cross_section_expr = solve(law, macroscopic_cross_section, dict=True)[0][macroscopic_cross_section]
    result_expr = result_cross_section_expr.subs(mean_free_path, mean_free_path_)
    return expr_to_quantity(result_expr)
