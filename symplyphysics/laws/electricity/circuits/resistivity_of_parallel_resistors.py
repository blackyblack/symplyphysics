from sympy import Eq, Idx, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    SymbolIndexed,
    SumIndexed,
    global_index,
)
from symplyphysics.core.expr_comparisons import expr_equals
import symplyphysics.laws.electricity.circuits.conductivity_of_parallel_resistors as parallel_conductivity
import symplyphysics.definitions.electrical_conductivity_is_inversed_resistance as conductance_definition

# Description
## If resistors are connected in parallel, total resistance is the inverse sum of conductance of each resistor.
## Law: R_parallel = 1/sum(sigma_i), where
## R_parallel - total resistance,
## sigma_i - conductance of i-th resistor.

parallel_resistance = Symbol("parallel_resistance", units.impedance)
conductance = SymbolIndexed("conductance", units.conductance)
law = Eq(parallel_resistance, 1 / SumIndexed(conductance[global_index], global_index))

# Derive the law from the conductivity law for parallel resistors

# Deriving the law for two resistors. Unfortunately, it is not possible to formally derive this for an arbitrary
# number of resistors due to the proof limitation of sympy. But it is possible to follow the following technique
# for any other number of resistors to prove the given equivalence.

resistance1 = Symbol("resistance1", units.impedance)
resistance2 = Symbol("resistance2", units.impedance)

conductance1 = solve(
    conductance_definition.definition.subs(conductance_definition.object_resistance, resistance1),
    conductance_definition.object_conductivity)[0]

conductance2 = solve(
    conductance_definition.definition.subs(conductance_definition.object_resistance, resistance2),
    conductance_definition.object_conductivity)[0]

local_index_ = Idx("index_local_", (1, 2))
parallel_conductance_law = parallel_conductivity.law.subs(global_index, local_index_)
parallel_conductance_law = parallel_conductance_law.doit()
parallel_conductance_law = parallel_conductance_law.subs({
    parallel_conductivity.conductance[1]: conductance1,
    parallel_conductivity.conductance[2]: conductance2,
})
parallel_conductance = solve(parallel_conductance_law,
    parallel_conductivity.parallel_conductance)[0]

parallel_resistance_from_conductivity_law = solve(
    conductance_definition.definition.subs(conductance_definition.object_conductivity,
    parallel_conductance), conductance_definition.object_resistance)[0]

parallel_resistance_law = law.subs(global_index, local_index_)
parallel_resistance_law = parallel_resistance_law.doit()
parallel_resistance_law = parallel_resistance_law.subs({
    conductance[1]: conductance1,
    conductance[2]: conductance2,
})
parallel_resistance_from_law_in_question = solve(parallel_resistance_law, parallel_resistance)[0]

assert expr_equals(parallel_resistance_from_conductivity_law,
    parallel_resistance_from_law_in_question)


def print_law() -> str:
    return print_expression(law)


@validate_input(resistances_=units.impedance)
@validate_output(parallel_resistance)
def calculate_parallel_resistance(resistances_: list[Quantity]) -> Quantity:
    conductances_ = tuple(
        conductance_definition.calculate_conductivity(resistance) for resistance in resistances_)
    local_index = Idx("index_local", (1, len(resistances_)))
    resistances_law = law.subs(global_index, local_index)
    resistances_law = resistances_law.doit()
    solved = solve(resistances_law, parallel_resistance, dict=True)[0][parallel_resistance]
    for i, v in enumerate(conductances_):
        solved = solved.subs(conductance[i + 1], v)
    return Quantity(solved)
