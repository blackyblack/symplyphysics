r"""
Statistical weight of macrostate
================================

If a physical system can be described as having several states which can be occupied by
different numbers of particles but with the total number of particles being conserved and
a condition that all allowed microstates of the closed system are equiprobable, the formula
for the statistical weight of the system can be found in combinatorics.

**Notes:**

#. Law can also be represented in form :math:`\\Omega = \frac{N!}{\prod_i (N_i!)}`
   (:code:`Omega = factorial(N) / Product(factorial(N_i), i)`)

**Links:**

#. `Chemistry LibreTexts, formula 1.5.1 <https://chem.libretexts.org/Courses/Western_Washington_University/Biophysical_Chemistry_(Smirnov_and_McCarty)/01%3A_Biochemical_Thermodynamics/1.05%3A_The_Boltzmann_Distribution_and_the_Statistical_Definition_of_Entropy>`__.
"""

from typing import Sequence
from sympy import Eq, Idx, factorial
from symplyphysics import (
    validate_input,
    validate_output,
    global_index,
    ProductIndexed,
    IndexedSum,
    symbols,
)
from symplyphysics.core.symbols.symbols import clone_as_indexed

statistical_weight = symbols.statistical_weight
"""
:symbols:`statistical_weight` of the system's macrostate.
"""

particle_count_in_state = clone_as_indexed(symbols.particle_count)
"""
:symbols:`particle_count` in state :math:`i`.
"""

law = Eq(
    statistical_weight,
    factorial(IndexedSum(particle_count_in_state[global_index], global_index)) /
    ProductIndexed(factorial(particle_count_in_state[global_index]), global_index),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(particle_counts_=particle_count_in_state)
@validate_output(statistical_weight)
def calculate_statistical_weight(particle_counts_: Sequence[int]) -> int:
    local_index = Idx("local_index", (1, len(particle_counts_)))
    result = law.rhs.subs(global_index, local_index).doit()
    for idx_, particle_count_ in enumerate(particle_counts_, 1):
        result = result.subs(particle_count_in_state[idx_], particle_count_)
    return int(result)
