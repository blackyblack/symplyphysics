from sympy import (Eq, solve)
from symplyphysics.definitions import electrical_conductivity_is_inversed_resistance as conductance_definition
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## If two resistors are connected in parallel, total conductance is sum of conductances of each resistor.
## Law: sigma_parallel = sigma1 + sigma2, where
## sigma_parallel is total conductance,
## sigma1 and sigma2 are conductances of two parallel resistors.

parallel_conductance = Symbol("parallel_conductance", units.conductance)
first_conductance = Symbol("first_conductance", units.conductance)
second_conductance = Symbol("second_conductance", units.conductance)

law = Eq(parallel_conductance, first_conductance + second_conductance)

def print_law() -> str:
    return print_expression(law)

first_resistance = Symbol("first_resistance", units.impedance)
second_resistance = Symbol("second_resistance", units.impedance)
parallel_resistance = Symbol("parallel_resistance", units.impedance)

@validate_input(first_resistance_=first_resistance, second_resistance_=second_resistance)
@validate_output(parallel_resistance)
def calculate_resistance(first_resistance_: Quantity, second_resistance_: Quantity) -> Quantity:
    conductance1 = solve(conductance_definition.definition, conductance_definition.object_conductivity, dict=True)[0][conductance_definition.object_conductivity].subs(
        {conductance_definition.object_resistance: first_resistance})
    conductance2 = solve(conductance_definition.definition, conductance_definition.object_conductivity, dict=True)[0][conductance_definition.object_conductivity].subs(
        {conductance_definition.object_resistance: second_resistance})
    result_conductance_expr = solve(law, parallel_conductance, dict=True)[0][parallel_conductance].subs({first_conductance: conductance1, second_conductance: conductance2})
    result_resistance = solve(conductance_definition.definition, conductance_definition.object_resistance, dict=True)[0][conductance_definition.object_resistance].subs(
        {conductance_definition.object_conductivity: result_conductance_expr})

    result_expr = result_resistance.subs({
        first_resistance: first_resistance_,
        second_resistance: second_resistance_
    })
    return Quantity(result_expr)
