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
## The Fermi-Dirac distribution describes the average number of particles in some energy state.
## It is normalized by the condition that the sum of the number of particles in each energy state
## should be equal to the total number of particles in the system.

# Law: N = Sum(N_i, i)
## N - total number of particles
## N_i - average number of fermions in single-particle state `i`
## i - energy state index

# Note
## This can be used to express chemical potential `mu` as a function of temperature `T` and total number
## of particles `N`.

total_fermion_count = Symbol("total_fermion_count", dimensionless)
average_fermion_count = SymbolIndexed("average_fermion_count", dimensionless)
energy_state_index = global_index

law = Eq(total_fermion_count, SumIndexed(average_fermion_count[energy_state_index], energy_state_index))


def print_law() -> str:
    return print_expression(law)


@validate_input(average_fermion_counts_=average_fermion_count)
@validate_output(total_fermion_count)
def calculate_total_fermion_count(average_fermion_counts_: Sequence[float]) -> int:
    energy_state_index_ = Idx("state_index", (1, len(average_fermion_counts_)))
    result = law.rhs.subs(energy_state_index, energy_state_index_).doit()
    for idx_, count_ in enumerate(average_fermion_counts_, 1):
        result = result.subs(average_fermion_count[idx_], count_)
    int_result = int(result)
    if int_result != result:
        raise ValueError(f"The total number of particles should be an integer, not {result}")
    return int_result
