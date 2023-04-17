from sympy import Expr
from sympy.physics.units.dimensions import Dimension
from symplyphysics import (
    Eq, pretty, solve, units, S, Probability, expr_to_quantity, convert_to
)
from symplyphysics.core.quantity_decorator import validate_input_symbols
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

# Description
## Thermal neutron utilization factor (f), is the ratio of the number of neutrons absorbed in the fuel
## versus the total number of absorptions everywhere in the reactor, i.e., in the fuel, moderator, cladding,
## and other reactor materials

## Law: f = Σa_fuel / Σa_total
## Where:
## Σa_fuel - macroscopic absorption cross-section of the fuel.
##   See [macroscopic cross-section](./macroscopic_cross_section_from_free_mean_path.py) implementation.
## Σa_total - macroscopic absorption cross-section of the fuel, moderator, cladding, etc, combined.
## f - thermal neutron utilisation factor

macroscopic_fuel_absorption_cross_section = Symbol("macroscopic_fuel_absorption_cross_section", 1 / units.length)
macroscopic_total_absorption_cross_section = Symbol("macroscopic_total_absorption_cross_section", 1 / units.length)
thermal_utilisation_factor = Symbol("thermal_utilisation_factor", Dimension(S.One))

law = Eq(thermal_utilisation_factor,
    macroscopic_fuel_absorption_cross_section / macroscopic_total_absorption_cross_section)

def print(expr: Expr) -> str:
    symbols = [macroscopic_fuel_absorption_cross_section, macroscopic_total_absorption_cross_section, thermal_utilisation_factor]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(
    macroscopic_fuel_absorption_cross_section_=macroscopic_fuel_absorption_cross_section,
    macroscopic_total_absorption_cross_section_=macroscopic_total_absorption_cross_section)
def calculate_utilisation_factor(
    macroscopic_fuel_absorption_cross_section_: Quantity,
    macroscopic_total_absorption_cross_section_: Quantity) -> Probability:
    if macroscopic_fuel_absorption_cross_section_.scale_factor > macroscopic_total_absorption_cross_section_.scale_factor:
        raise ValueError(f"macroscopic_fuel_absorption_cross_section_ ({macroscopic_fuel_absorption_cross_section_.scale_factor}) should be <= "
        f"macroscopic_total_absorption_cross_section_ ({macroscopic_total_absorption_cross_section_.scale_factor})")

    result_factor_expr = solve(law, thermal_utilisation_factor, dict=True)[0][thermal_utilisation_factor]
    result_expr = result_factor_expr.subs({
        macroscopic_fuel_absorption_cross_section: macroscopic_fuel_absorption_cross_section_,
        macroscopic_total_absorption_cross_section: macroscopic_total_absorption_cross_section_})
    result_factor = expr_to_quantity(result_expr)
    return Probability(convert_to(result_factor, S.One).evalf())
