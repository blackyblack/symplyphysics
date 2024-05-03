from typing import Iterable
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

# Description
## If a physical system can be described as having several states which can be occupied by
## different numbers of particles but with the total number of particles being conserved and
## a condition that all allowed microstates of the closed system are equiprobable, the formula
## for the statistical weight of the system can be found in combinatorics:

# Law: W = N! / Product((N_i)!, i) = Sum(N_i, i)! / Product((N_i)!, i)
## W - statistical weight of system's macrostate
## N - total number of particles
## N_i - number of particles of i-th state
## a! - factorial of `a`

statistical_weight = Symbol("statistical_weight", dimensionless, integer=True)
particle_count = SymbolIndexed("particle_count", dimensionless, integer=True)

law = Eq(
    statistical_weight,
    factorial(SumIndexed(particle_count[global_index], global_index))
    / ProductIndexed(factorial(particle_count[global_index]), global_index),
)


@validate_input(particle_counts_=particle_count)
@validate_output(statistical_weight)
def calculate_statistical_weight(particle_counts_: Iterable[int]) -> int:
    local_index = Idx("local_index", (1, len(particle_counts_)))
    result = law.rhs.subs(global_index, local_index).doit()
    for idx_, particle_count_ in enumerate(particle_counts_, 1):
        result = result.subs(particle_count[idx_], particle_count_)
    return int(result)
