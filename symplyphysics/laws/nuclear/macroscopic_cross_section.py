from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Macroscopic cross-section - represents the effective target area of all of the nuclei contained
## in the volume of the material (such as fuel pellet). It is the probability of neutron-nucleus interaction per
## centimeter of neutron travel.

## Macroscopic cross-section: Σ = 1 / λ
## Where:
## λ (mean free path) is equal to the average value of x, the distance traveled by a neutron without any
## interaction, over the interaction probability distribution.
## Σ is the macroscopic cross-section.

mean_free_path = symbols('mean_free_path')
macroscopic_cross_section = symbols('macroscopic_cross_section')

law = Eq(macroscopic_cross_section, 1 / mean_free_path)

def print():
    return pretty(law, use_unicode=False)

@validate_input(mean_free_path_=units.length)
@validate_output(1 / units.length)
def calculate_cross_section(mean_free_path_: Quantity) -> Quantity:
    result_cross_section_expr = solve(law, macroscopic_cross_section, dict=True)[0][macroscopic_cross_section]
    result_expr = result_cross_section_expr.subs({
        mean_free_path: mean_free_path_})
    return expr_to_quantity(result_expr, 'macro-cross-section')
