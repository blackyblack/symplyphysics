from sympy import (Eq, Idx, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, global_index)
from symplyphysics.definitions import electrical_conductivity_is_inversed_resistance as conductance_definition
from symplyphysics.laws.electricity.circuits import conductivity_of_parallel_resistors as parallel_resistors_law

# Description
## If two resistors are connected in parallel, total conductance is a sum of conductances of each resistor.
## See [conductivity_of_parallel_resistors](./conductivity_of_parallel_resistors.py) for additional information.
## Law: sigma_parallel = sigma1 + sigma2, where
## sigma_parallel is total conductance,
## sigma1 and sigma2 are conductances of two parallel resistors.

parallel_conductance = Symbol("parallel_conductance", units.conductance)
first_conductance = Symbol("first_conductance", units.conductance)
second_conductance = Symbol("second_conductance", units.conductance)

law = Eq(parallel_conductance, first_conductance + second_conductance)

# Derive the same law from more general law for any number of resistors

local_index_ = Idx("local_index_", (1, 2))
two_resistors_law = parallel_resistors_law.law.subs(global_index, local_index_).doit()
two_resistors_law = two_resistors_law.subs({
    parallel_resistors_law.conductance[1]: first_conductance,
    parallel_resistors_law.conductance[2]: second_conductance,
})
assert two_resistors_law.rhs == law.rhs


def print_law() -> str:
    return print_expression(law)


@validate_input(first_resistance_=units.impedance, second_resistance_=units.impedance)
@validate_output(units.impedance)
def calculate_resistance(first_resistance_: Quantity, second_resistance_: Quantity) -> Quantity:
    first_resistance = Symbol("first_resistance", units.impedance)
    second_resistance = Symbol("second_resistance", units.impedance)
    conductance1 = solve(conductance_definition.definition,
        conductance_definition.conductivity,
        dict=True)[0][conductance_definition.conductivity].subs(
        {conductance_definition.resistance: first_resistance})
    conductance2 = solve(conductance_definition.definition,
        conductance_definition.conductivity,
        dict=True)[0][conductance_definition.conductivity].subs(
        {conductance_definition.resistance: second_resistance})
    result_conductance_expr = solve(law, parallel_conductance,
        dict=True)[0][parallel_conductance].subs({
        first_conductance: conductance1,
        second_conductance: conductance2
        })
    result_resistance = solve(conductance_definition.definition,
        conductance_definition.resistance,
        dict=True)[0][conductance_definition.resistance].subs(
        {conductance_definition.conductivity: result_conductance_expr})

    result_expr = result_resistance.subs({
        first_resistance: first_resistance_,
        second_resistance: second_resistance_
    })
    return Quantity(result_expr)
