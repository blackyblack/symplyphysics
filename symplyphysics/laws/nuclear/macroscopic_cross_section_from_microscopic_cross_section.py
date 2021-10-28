from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Macroscopic cross-section - represents the effective target area of all of the nuclei contained
## in the volume of the material (such as fuel pellet). It is the probability of neutron-nucleus interaction per
## centimeter of neutron travel.

## Macroscopic cross-section: Σ = σ * N
## Where:
## σ (microscopic cross-section) is the effective `target area` that a nucleus interacts with an incident neutron.
## N (atomic number density) is the number of atoms of a given type per unit volume of the material.
## Σ is the macroscopic cross-section.

microscopic_cross_section = symbols('microscopic_cross_section')
atomic_number_density = symbols('atomic_number_density')
macroscopic_cross_section = symbols('macroscopic_cross_section')

law = Eq(macroscopic_cross_section, microscopic_cross_section * atomic_number_density)

def print():
    return pretty(law, use_unicode=False)

@validate_input(microscopic_cross_section_=units.length ** 2, atomic_number_density_=(1 / units.length ** 3))
@validate_output(1 / units.length)
def calculate_cross_section(microscopic_cross_section_: Quantity, atomic_number_density_: Quantity) -> Quantity:
    result_cross_section_expr = solve(law, macroscopic_cross_section, dict=True)[0][macroscopic_cross_section]
    result_expr = result_cross_section_expr.subs({
        microscopic_cross_section: microscopic_cross_section_,
        atomic_number_density: atomic_number_density_})
    return expr_to_quantity(result_expr, 'macro-cross-section')
