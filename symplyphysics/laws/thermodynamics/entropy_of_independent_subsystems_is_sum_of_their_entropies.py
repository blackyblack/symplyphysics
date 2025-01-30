"""
Entropy of independent subsystems is sum of their entropies
===========================================================

If a thermodynamic system can be decomposed into several subsystems which are all statistically
independent, the total entropy of the system can be calculated as the sum of the entropies of all
the subsystems. Mathematically speaking, this is a representation of such a property of entropy
known as `subadditivity <https://en.wikipedia.org/wiki/Subadditivity>`_.

**Conditions:**

#. The subsystems must be (approximately) independent in the statistical sense.

**Links:**

#. `ScienceDirect <https://www.sciencedirect.com/science/article/abs/pii/S0378437103002620>`__.
"""

from typing import Sequence
from sympy import Eq, Idx
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    SumIndexed,
    global_index,
    symbols,
)
from symplyphysics.core.symbols.symbols import clone_as_indexed

total_entropy = symbols.entropy
"""
Total :symbols:`entropy` of the system as a whole.
"""

subsystem_entropy = clone_as_indexed(symbols.entropy)
"""
:symbols:`entropy` of the :math:`i`-th subsystem.
"""

law = Eq(total_entropy, SumIndexed(subsystem_entropy[global_index], global_index))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(subsystem_entropies_=subsystem_entropy)
@validate_output(total_entropy)
def calculate_total_entropy(subsystem_entropies_: Sequence[Quantity]) -> Quantity:
    local_index = Idx("local_index", (1, len(subsystem_entropies_)))
    result = law.rhs.subs(global_index, local_index).doit()
    for idx, subsystem_entropy_ in enumerate(subsystem_entropies_, 1):
        result = result.subs(subsystem_entropy[idx], subsystem_entropy_)
    return Quantity(result)
