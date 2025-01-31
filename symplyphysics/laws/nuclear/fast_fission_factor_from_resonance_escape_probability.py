"""
Fast fission factor from resonance escape probability
=====================================================

The fast fission factor is the ratio of the fast neutrons produced by fissions at all
energies to the number of fast neutrons produced in thermal fission. Unfortunately
proper fast fission factor formula is too complex for actual use. Approximate formula is
not very useful since there is no actual data on :math:`P_\\text{FAF}`,
:math:`P_\\text{TAF}`, and :math:`u_\\text{f}` values.

**Links:**

#. `Wikipedia, fourth row in table <https://en.wikipedia.org/wiki/Six_factor_formula>`__.
"""

from sympy import Eq
from symplyphysics import symbols, clone_as_symbol

fast_fission_factor = symbols.fast_fission_factor
"""
:symbols:`fast_fission_factor`.
"""

resonance_escape_probability = symbols.resonance_escape_probability
"""
:symbols:`resonance_escape_probability`.
"""

thermal_utilization_factor = symbols.thermal_utilization_factor
"""
:symbols:`thermal_utilization_factor`.
"""

fast_utilization = symbols.fast_utilization
"""
:symbols:`fast_utilization`.
"""

fast_neutrons_per_fission = clone_as_symbol(symbols.particle_count, display_symbol="nu_f", display_latex="\\nu_\\text{f}")
"""
Number of fast neutrons produced per fission. See :symbols:`particle_count`.
"""

thermal_neutrons_per_fission = clone_as_symbol(symbols.particle_count, display_symbol="nu_t", display_latex="\\nu_\\text{t}")
"""
Number of thermal neutrons produced per fission. See :symbols:`particle_count`.
"""

fast_absorption_fission_probability = symbols.fast_absorption_fission_probability
"""
:symbols:`fast_absorption_fission_probability`.
"""

thermal_absorption_fission_probability = symbols.thermal_absorption_fission_probability
"""
:symbols:`thermal_absorption_fission_probability`.
"""

thermal_non_leakage_probability = symbols.thermal_non_leakage_probability
"""
:symbols:`thermal_non_leakage_probability`.
"""

_first_ratio = (1 - resonance_escape_probability) / resonance_escape_probability
_second_ratio = (
    (fast_utilization * fast_neutrons_per_fission * fast_absorption_fission_probability)
    / (thermal_utilization_factor * thermal_neutrons_per_fission * thermal_absorption_fission_probability * thermal_non_leakage_probability)
)

law = Eq(fast_fission_factor, 1 + _first_ratio * _second_ratio,)
"""
:laws:symbol::

:laws:latex::
"""

# There is no calculate() method for the lack of data and use of this formula.
