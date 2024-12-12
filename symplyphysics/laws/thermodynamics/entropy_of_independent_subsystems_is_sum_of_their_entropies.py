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

#. `Wikipedia <https://www.sciencedirect.com/science/article/abs/pii/S0378437103002620>`__.
"""

from typing import Sequence
from sympy import Eq, Idx
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    SumIndexed,
    SymbolIndexed,
    global_index,
)

total_entropy = Symbol("total_entropy", units.energy / units.temperature)
"""
Total entropy of the system as a whole.

Symbol:
    :code:`S`
"""

subsystem_entropy = SymbolIndexed("subsystem_entropy", units.energy / units.temperature)
r"""
Entropy of the :math:`i`-th subsystem.

Symbol:
    :code:`S_i`

Latex:
    :math:`S_i`
"""

law = Eq(total_entropy, SumIndexed(subsystem_entropy[global_index], global_index))
r"""
:code:`S = Sum(S_i, i)`

Latex:
    .. math::
        S = \sum_i S_i
"""


@validate_input(subsystem_entropies_=subsystem_entropy)
@validate_output(total_entropy)
def calculate_total_entropy(subsystem_entropies_: Sequence[Quantity]) -> Quantity:
    local_index = Idx("local_index", (1, len(subsystem_entropies_)))
    result = law.rhs.subs(global_index, local_index).doit()
    for idx, subsystem_entropy_ in enumerate(subsystem_entropies_, 1):
        result = result.subs(subsystem_entropy[idx], subsystem_entropy_)
    return Quantity(result)
