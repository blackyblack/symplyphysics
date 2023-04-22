from sympy import (Eq, solve)
from symplyphysics import (
    units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol,
)

# Description
## Macroscopic cross-section - represents the effective target area of all of the nuclei contained
## in the volume of the material (such as fuel pellet). It is the probability of neutron-nucleus interaction per
## centimeter of neutron travel.

## Macroscopic cross-section: Σ = σ * N
## Where:
## σ (microscopic cross-section) is the effective target area that a nucleus interacts with an incident neutron.
## N (atomic number density) is the number of atoms of a given type per unit volume of the material.
##   See [atomic number density](../chemistry/atomic_number_density_from_material_density_atomic_weight.py) implementation.
## Σ is the macroscopic cross-section.

microscopic_cross_section = Symbol("microscopic_cross_section", units.length**2)
atomic_number_density = Symbol("atomic_number_density", 1 / units.length**3)
macroscopic_cross_section = Symbol("macroscopic_cross_section", 1 / units.length)

law = Eq(macroscopic_cross_section, microscopic_cross_section * atomic_number_density)

def print() -> str:
    return print_expression(law)

@validate_input_symbols(microscopic_cross_section_=microscopic_cross_section, atomic_number_density_=atomic_number_density)
@validate_output_symbol(macroscopic_cross_section)
def calculate_cross_section(microscopic_cross_section_: Quantity, atomic_number_density_: Quantity) -> Quantity:
    result_cross_section_expr = solve(law, macroscopic_cross_section, dict=True)[0][macroscopic_cross_section]
    result_expr = result_cross_section_expr.subs({
        microscopic_cross_section: microscopic_cross_section_,
        atomic_number_density: atomic_number_density_})
    return expr_to_quantity(result_expr)
