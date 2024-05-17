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

# Description
## If a thermodynamic system can be decomposed into several subsystems which are all statistically
## independent, the total entropy of the system can be calculated as the sum of the entropies of all
## the subsystems. Mathematically speaking, this is a representation of such a property of entropy
## known as [subadditivity](https://en.wikipedia.org/wiki/Subadditivity).

# Law: S = Sum(S_i, i)
## S - total entropy of system
## S_i - entropy of i-th subsystem

# Conditions
## - The subsystems are (approximately) independent in the statistical sense

total_entropy = Symbol("total_entropy", units.energy / units.temperature)
subsystem_entropy = SymbolIndexed("subsystem_entropy", units.energy / units.temperature)

law = Eq(total_entropy, SumIndexed(subsystem_entropy[global_index], global_index))


@validate_input(subsystem_entropies_=subsystem_entropy)
@validate_output(total_entropy)
def calculate_total_entropy(subsystem_entropies_: Sequence[Quantity]) -> Quantity:
    local_index = Idx("local_index", (1, len(subsystem_entropies_)))
    result = law.rhs.subs(global_index, local_index).doit()
    for idx, subsystem_entropy_ in enumerate(subsystem_entropies_, 1):
        result = result.subs(subsystem_entropy[idx], subsystem_entropy_)
    return Quantity(result)
