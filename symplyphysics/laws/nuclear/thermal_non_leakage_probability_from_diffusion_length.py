"""
Thermal non-leakage probability from diffusion area and geometric buckling
==========================================================================

Thermal non-leakage probability can be found from the diffusion area of thermal neutrons
and geometric buckling.

**Links:**

#. `NuclearPower <https://www.nuclear-power.com/nuclear-power/reactor-physics/nuclear-fission-chain-reaction/thermal-non-leakage-probability/>`__.
#. `Wikipedia, last row in table <https://en.wikipedia.org/wiki/Six_factor_formula>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    convert_to_float,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.symbols.probability import Probability

thermal_diffusion_area = clone_as_symbol(symbols.neutron_diffusion_area, display_symbol="L_th^2", display_latex="L_\\text{th}^2")
"""
Thermal :symbols:`neutron_diffusion_area`.
"""

geometric_buckling = symbols.geometric_buckling
"""
:symbols:`geometric_buckling`.
"""

thermal_non_leakage_probability = symbols.thermal_non_leakage_probability
"""
:symbols:`thermal_non_leakage_probability`.
"""

law = Eq(thermal_non_leakage_probability, 1 / (1 + thermal_diffusion_area * geometric_buckling))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(thermal_diffusion_area_=thermal_diffusion_area,
    geometric_buckling_=geometric_buckling)
@validate_output(thermal_non_leakage_probability)
def calculate_probability(thermal_diffusion_area_: Quantity,
    geometric_buckling_: Quantity) -> Probability:
    result_probability_expr = solve(law, thermal_non_leakage_probability,
        dict=True)[0][thermal_non_leakage_probability]
    result_expr = result_probability_expr.subs({
        thermal_diffusion_area: thermal_diffusion_area_,
        geometric_buckling: geometric_buckling_
    })
    return Probability(convert_to_float(result_expr))
