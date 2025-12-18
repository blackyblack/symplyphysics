"""
Fast non-leakage probability from Fermi age and geometric buckling
==================================================================

The fast non-leakage probability can be calculated from the Fermi age of neutrons and
the geometric buckling associated with the given nuclear reactor.

**Links:**

#. `Wikipedia, fifth row in table <https://en.wikipedia.org/wiki/Six_factor_formula>`__.
"""

from sympy import Eq, solve, exp
from symplyphysics import (
    Quantity,
    convert_to_float,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.symbols.probability import Probability

geometric_buckling = symbols.geometric_buckling
"""
:symbols:`geometric_buckling`.
"""

fermi_age = symbols.neutron_fermi_age
"""
:symbols:`neutron_fermi_age`.
"""

fast_non_leakage_probability = symbols.fast_non_leakage_probability
"""
:symbols:`fast_non_leakage_probability`.
"""

law = Eq(fast_non_leakage_probability, exp(-1 * geometric_buckling * fermi_age))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(geometric_buckling_=geometric_buckling, neutron_fermi_age_=fermi_age)
@validate_output(fast_non_leakage_probability)
def calculate_probability(geometric_buckling_: Quantity,
    neutron_fermi_age_: Quantity) -> Probability:
    result_probability_expr = solve(law, fast_non_leakage_probability,
        dict=True)[0][fast_non_leakage_probability]
    result_expr = result_probability_expr.subs({
        geometric_buckling: geometric_buckling_,
        fermi_age: neutron_fermi_age_
    })
    return Probability(convert_to_float(result_expr))
