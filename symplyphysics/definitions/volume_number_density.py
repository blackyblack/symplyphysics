from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Volume number density is the number of specified objects per unit volume.

## Definition: N = n / V
## Where:
## n is the total number of objects
## V is volume

number_density, objects, volume = symbols('number_density objects volume')
definition = Eq(number_density, objects / volume)

definition_dimension_SI = 1 / units.meter**3

def print():
    return pretty(definition, use_unicode=False)

def print_dimension():
    return pretty(definition_dimension_SI, use_unicode=False)

@validate_input(volume_=units.volume)
@validate_output(1 / units.volume)
def calculate_number_density(objects_: int, volume_: Quantity) -> Quantity:
    solved = solve(definition, number_density, dict=True)[0][number_density]
    result_expr = solved.subs({
        objects: objects_,
        volume: volume_})
    return expr_to_quantity(result_expr, 'volume_number_density')