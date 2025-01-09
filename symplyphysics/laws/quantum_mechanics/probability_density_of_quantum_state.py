"""
Probability density of quantum state
====================================

Probability density of a quantum system is the likelihood that the system will be in a particular
position at a particular time.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Wave_function#Position-space_wave_functions>`__.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols

probability_density = symbols.probability_density
"""
:symbols:`probability_density` of the system.
"""

wave_function = symbols.wave_function
"""
:symbols:`wave_function` of the system.
"""

law = Eq(probability_density, abs(wave_function)**2)
"""
:laws:symbols::

:laws:latex::
"""


@validate_input(wave_function_value_=wave_function)
@validate_output(probability_density)
def calculate_probability_density(wave_function_value_: Quantity) -> Quantity:
    result = law.rhs.subs(wave_function, wave_function_value_)
    return Quantity(result)
