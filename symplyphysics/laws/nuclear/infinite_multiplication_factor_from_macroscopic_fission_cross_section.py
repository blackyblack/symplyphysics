from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units, S,
    validate_input, expr_to_quantity, convert_to
)

# Description
## Infinite multiplication factor: k_infinite = v * Σf / Σa
## Where:
## v - average number of neutrons produced per fission.
## Σf - overall macroscopic fission cross-section.
##   See [macroscopic cross-section](./macroscopic_cross_section_from_free_mean_path.py) implementation.
## Σa - overall macroscopic absorption cross-section.
## k_infinite - infinite multiplication factor.
##   See [infinite multiplication factor](./infinite_multiplication_factor.py)

neutrons_per_fission = symbols('neutrons_per_fission')
macroscopic_fission_cross_section = symbols('macroscopic_fission_cross_section')
macroscopic_absorption_cross_section = symbols('macroscopic_absorption_cross_section')
infinite_multiplication_factor = symbols('infinite_multiplication_factor')

law = Eq(infinite_multiplication_factor,
    neutrons_per_fission * macroscopic_fission_cross_section / macroscopic_absorption_cross_section)

def print():
    return pretty(law, use_unicode=False)

@validate_input(macroscopic_fission_cross_section_=(1 / units.length), macroscopic_absorption_cross_section_=(1 / units.length))
def calculate_multiplication_factor(neutrons_per_fission_: float,
    macroscopic_fission_cross_section_: Quantity,
    macroscopic_absorption_cross_section_: Quantity) -> float:

    result_factor_expr = solve(law, infinite_multiplication_factor, dict=True)[0][infinite_multiplication_factor]
    result_expr = result_factor_expr.subs({
        neutrons_per_fission: neutrons_per_fission_,
        macroscopic_fission_cross_section: macroscopic_fission_cross_section_,
        macroscopic_absorption_cross_section: macroscopic_absorption_cross_section_})
    result_factor = expr_to_quantity(result_expr, 'inf_multiplication_factor')
    return convert_to(result_factor, S.One).n()
