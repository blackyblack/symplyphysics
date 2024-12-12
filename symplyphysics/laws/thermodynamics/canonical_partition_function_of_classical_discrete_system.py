"""
Canonical partition function of a classical discrete system
===========================================================

Let us assume a canonical ensemble, i.e. a thermodynamically large system that is in thermal contact
with the environment, with a temperature T and whose volume and number of constituent particles remain
constant. For a classical discrete system the partition function is the sum of the Boltzmann factors
of all the possible energy states:

**Links:**

#. `Wikipedia, derivable from here <https://en.wikipedia.org/wiki/Partition_function_(statistical_mechanics)#Classical_discrete_system>`__.

..
    NOTE replce `boltzmann_factor` with actual exponent?
"""

from typing import Sequence
from sympy import Eq, Idx, solve
from symplyphysics import (
    convert_to_float,
    validate_input,
    validate_output,
    SumIndexed,
    global_index,
    symbols,
)
from symplyphysics.core.symbols.symbols import clone_as_indexed

partition_function = symbols.partition_function
"""
:symbols:`partition_function` of the system.
"""

boltzmann_factor = clone_as_indexed(symbols.boltzmann_factor)
"""
:symbols:`boltzmann_factor` of energy state :math:`i`. See :doc:`Boltzmann factor
<definitions.boltzmann_factor_via_state_energy_and_temperature>`.
"""

law = Eq(partition_function, SumIndexed(boltzmann_factor[global_index], global_index))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(boltzmann_factors_=boltzmann_factor)
@validate_output(partition_function)
def calculate_partition_function(boltzmann_factors_: Sequence[float]) -> float:
    local_index = Idx("index_local", (1, len(boltzmann_factors_)))
    partition_function_law = law.subs(global_index, local_index)
    partition_function_law = partition_function_law.doit()
    solved = solve(partition_function_law, partition_function, dict=True)[0][partition_function]
    for i, v in enumerate(boltzmann_factors_):
        solved = solved.subs(boltzmann_factor[i + 1], v)
    return convert_to_float(solved)
