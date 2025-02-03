"""
Infinite multiplication factor formula
======================================

Infinite multiplication factor can be found using thermal fission factor, fast fission
factor, resonance escape probability, and thermal utilization factor.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Nuclear_chain_reaction#Six-factor_formula>`__.
# `NuclearPower <https://www.nuclear-power.com/nuclear-power/reactor-physics/nuclear-fission-chain-reaction/four-factor-formula-infinite-multiplication-factor/>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import symbols
from symplyphysics.core.symbols.probability import Probability

thermal_fission_factor = symbols.thermal_fission_factor
"""
:symbols:`thermal_fission_factor`.
"""

thermal_utilization_factor = symbols.thermal_utilization_factor
"""
:symbols:`thermal_utilization_factor`.
"""

resonance_escape_probability = symbols.resonance_escape_probability
"""
:symbols:`resonance_escape_probability`.
"""

fast_fission_factor = symbols.fast_fission_factor
"""
:symbols:`fast_fission_factor`.
"""

infinite_multiplication_factor = symbols.infinite_multiplication_factor
"""
:symbols:`infinite_multiplication_factor`.
"""

law = Eq(infinite_multiplication_factor,
    thermal_fission_factor * fast_fission_factor * resonance_escape_probability * thermal_utilization_factor)
"""
:laws:symbol::

:laws:latex::
"""


def calculate_multiplication_factor(neutron_reproduction_: float, fast_fission_: float,
    resonance_escape_probability_: Probability, thermal_utilisation_: Probability) -> float:

    result_factor_expr = solve(law, infinite_multiplication_factor,
        dict=True)[0][infinite_multiplication_factor]
    result_expr = result_factor_expr.subs({
        thermal_fission_factor: neutron_reproduction_,
        fast_fission_factor: fast_fission_,
        resonance_escape_probability: resonance_escape_probability_,
        thermal_utilization_factor: thermal_utilisation_
    })
    return result_expr.evalf()
