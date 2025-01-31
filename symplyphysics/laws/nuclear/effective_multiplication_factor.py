"""
Effective multiplication factor from infinite multiplication factor and probabilities
=====================================================================================

The effective multiplication factor can be found from the infinite multiplication factor
as well as fast and thermal non-leakage probabilities.

**Links:**

#. `Wikipedia, formula above table <https://en.wikipedia.org/wiki/Six_factor_formula>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import symbols
from symplyphysics.core.symbols.probability import Probability

fast_non_leakage_probability = symbols.fast_non_leakage_probability
"""
:symbols:`fast_non_leakage_probability`.
"""

thermal_non_leakage_probability = symbols.thermal_non_leakage_probability
"""
:symbols:`thermal_non_leakage_probability`.
"""

infinite_multiplication_factor = symbols.infinite_multiplication_factor
"""
:symbols:`infinite_multiplication_factor`.
"""

effective_multiplication_factor = symbols.effective_multiplication_factor
"""
:symbols:`effective_multiplication_factor`.
"""

law = Eq(
    effective_multiplication_factor,
    infinite_multiplication_factor * fast_non_leakage_probability * thermal_non_leakage_probability)
"""
:laws:symbol::

:laws:latex::
"""


def calculate_multiplication_factor(infinite_multiplication_factor_: float,
    fast_non_leakage_probability_: Probability,
    thermal_non_leakage_probability_: Probability) -> float:

    result_factor_expr = solve(law, effective_multiplication_factor,
        dict=True)[0][effective_multiplication_factor]
    result_expr = result_factor_expr.subs({
        infinite_multiplication_factor: infinite_multiplication_factor_,
        fast_non_leakage_probability: fast_non_leakage_probability_,
        thermal_non_leakage_probability: thermal_non_leakage_probability_
    })
    return result_expr.evalf()
