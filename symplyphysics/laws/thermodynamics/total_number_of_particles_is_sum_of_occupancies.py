from typing import Sequence
from sympy import Eq, Idx
from symplyphysics import (
    dimensionless,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    global_index,
    SymbolIndexed,
    SumIndexed,
)

# Description
## Maxwell-Boltzmann, Fermi-Dirac and Bose-Einstein distributions describe the average number of particles
## in some energy level or particle state. They are normalized by the condition that the sum of the number
## of particles in each energy level or particle state should be equal to the total number of particles
## in the system.

# Law: N = Sum(n_i, i)
## N - total number of particles
## n_i - occupancy of energy level or particle state of index `i`
## i - index

# Note
## This can be used to express chemical potential `mu` as a function of temperature `T` and total number
## of particles `N` for Fermi-Dirac and Bose-Einstein distributions.

total_particle_count = Symbol("total_particle_count", dimensionless)
occupancy = SymbolIndexed("occupancy", dimensionless)

law = Eq(total_particle_count, SumIndexed(occupancy[global_index], global_index))


def print_law() -> str:
    return print_expression(law)


@validate_input(occupancies_=occupancy)
@validate_output(total_particle_count)
def calculate_total_particle_count(occupancies_: Sequence[float]) -> int:
    local_index_ = Idx("local_index_", (1, len(occupancies_)))
    result = law.rhs.subs(global_index, local_index_).doit()
    for idx_, count_ in enumerate(occupancies_, 1):
        result = result.subs(occupancy[idx_], count_)
    return int(result)
