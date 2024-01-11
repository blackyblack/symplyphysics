from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.symbols.symbols import tuple_of_symbols
from symplyphysics.core.operations.sum_array import SumArray
from symplyphysics.core.symbols.symbols import tuple_of_symbols
import symplyphysics.laws.electricity.circuits.conductivity_of_parallel_resistors as parallel_conductivity
import symplyphysics.definitions.electrical_conductivity_is_inversed_resistance as conductance_definition

# Description
## If resistors are connected in parallel, total resistance is the inverse sum of conductance of each resistor.
## Law: R_parallel = 1/sum(sigma_i), where
## R_parallel - total resistance,
## sigma_i - conductance of i-th resistor.

conductances = Symbol("conductances", units.conductance)
parallel_resistance = Symbol("parallel_resistance", units.impedance)
law = Eq(parallel_resistance, 1 / SumArray(conductances), evaluate=False)


# Derive the law from the conductivity law for parallel resistors

# Deriving the law for two resistors. Unfortunately, it is not possible to formally derive this for an arbitrary
# number of resistors due to the proof limitation of sympy. But it is possible to follow the following technique
# for any other number of resistors to prove the given equivalence.

resistance1, resistance2 = tuple_of_symbols("resistance", units.impedance, 2)

conductance1 = solve(
    conductance_definition.definition.subs(
        conductance_definition.object_resistance,
        resistance1
    ), 
    conductance_definition.object_conductivity
)[0]

conductance2 = solve(
    conductance_definition.definition.subs(
        conductance_definition.object_resistance,
        resistance2
    ),
    conductance_definition.object_conductivity
)[0]

parallel_conductance = solve(
    parallel_conductivity.law.subs(
        parallel_conductivity.conductances,
        (conductance1, conductance2)
    ),
    parallel_conductivity.parallel_conductance
)[0]

parallel_resistance_from_conductivity_law = solve(
    conductance_definition.definition.subs(
        conductance_definition.object_conductivity,
        parallel_conductance
    ),
    conductance_definition.object_resistance
)[0]

parallel_resistance_from_law_in_question = solve(
    law.subs(conductances, (conductance1, conductance2)),
    parallel_resistance
)[0]

assert expr_equals(
    parallel_resistance_from_conductivity_law, 
    parallel_resistance_from_law_in_question
)


def print_law() -> str:
    return print_expression(law)


@validate_input(resistances_=units.impedance)
@validate_output(units.impedance)
def calculate_parallel_resistance(resistances_: list[Quantity]) -> Quantity:
    conductances_ = tuple(
        conductance_definition.calculate_conductivity(resistance) 
        for resistance in resistances_
    )
    conductance_symbols = tuple_of_symbols("conductance", units.conductance, len(conductances_))
    resistances_law = law.subs(conductances, conductance_symbols).doit()
    solved = solve(resistances_law, parallel_resistance, dict=True)[0][parallel_resistance]
    for symbol, value in zip(conductance_symbols, conductances_):
        solved = solved.subs(symbol, value)
    return Quantity(solved)
