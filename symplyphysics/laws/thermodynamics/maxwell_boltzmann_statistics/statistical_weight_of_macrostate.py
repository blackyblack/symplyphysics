"""
Statistical weight of macrostate
================================

If a physical system can be described as having several states which can be occupied by
different numbers of particles but with the total number of particles being conserved and
a condition that all allowed microstates of the closed system are equiprobable, the formula
for the statistical weight of the system can be found in combinatorics.
"""

from typing import Sequence
from sympy import Eq, Idx, factorial
from symplyphysics import (
    Symbol,
    dimensionless,
    validate_input,
    validate_output,
    global_index,
    ProductIndexed,
    SumIndexed,
    SymbolIndexed,
)

statistical_weight = Symbol("statistical_weight", dimensionless, integer=True)
"""
Statistical weight of the system's macrostate.

Symbol:
    :code:`W`
"""

particle_count_in_state = SymbolIndexed("particle_count_in_state", dimensionless, integer=True)
r"""
Number of particles in state :math:`i`.

Symbol:
    :code:`N_i`

Latex:
    :math:`N_i`
"""

law = Eq(
    statistical_weight,
    factorial(SumIndexed(particle_count_in_state[global_index], global_index)) /
    ProductIndexed(factorial(particle_count_in_state[global_index]), global_index),
)
r"""
:code:`W = factorial(N) / Product(factorial(N_i), i) = factorial(Sum(N_i, i)) / Product(factorial(N_i), i)`

Latex:
    .. math::
        W = \frac{N}{\prod_i (N_i!)} = \frac{(\sum_i N_i)!}{\prod_i (N_i!)}
"""


@validate_input(particle_counts_=particle_count_in_state)
@validate_output(statistical_weight)
def calculate_statistical_weight(particle_counts_: Sequence[int]) -> int:
    local_index = Idx("local_index", (1, len(particle_counts_)))
    result = law.rhs.subs(global_index, local_index).doit()
    for idx_, particle_count_ in enumerate(particle_counts_, 1):
        result = result.subs(particle_count_in_state[idx_], particle_count_)
    return int(result)
