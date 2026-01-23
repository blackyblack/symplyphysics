"""
Quality factor via damping ratio
================================

In case of a damped oscillating system where the damping force is linearly proportional
to the oscillator's frequency, the :math:`Q` factor of the oscillator is inversely
proportional to the system's damping ratio. Hence greater values of the :math:`Q` factor
correspond to lower values of the damping ratio and to lower damping.

**Links:**

#. `Wikipedia, derivable from here <https://en.wikipedia.org/wiki/Damping#Q_factor_and_decay_rate>`__.
"""

from sympy import Eq
from symplyphysics import convert_to_float, validate_input, validate_output, symbols

quality_factor = symbols.quality_factor
"""
:symbols:`quality_factor` of the oscillator.
"""

damping_ratio = symbols.damping_ratio
"""
:symbols:`damping_ratio` of the oscillator. Also see :ref:`Damped harmonic oscillator equation`.
"""

law = Eq(quality_factor, 1 / (2 * damping_ratio))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(damping_ratio_=damping_ratio)
@validate_output(quality_factor)
def calculate_quality_factor(damping_ratio_: float) -> float:
    result = law.rhs.subs(damping_ratio, damping_ratio_)
    return convert_to_float(result)
