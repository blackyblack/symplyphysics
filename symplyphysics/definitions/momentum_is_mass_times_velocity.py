from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Momentum is the multiplication of velocity and mass

## Definition: M = m * V
## Where:
## m is a mass of the object
## V is it's velocity

momentum, mass, velocity = symbols('momentum mass velocity')
definition = Eq(momentum, mass * velocity)

definition_dimension_SI = units.kilogram * units.meter / units.second

def print():
    return pretty(definition, use_unicode=False)

def print_dimension():
    return pretty(definition_dimension_SI, use_unicode=False)

@validate_input(velocity_=units.velocity)
@validate_input(mass_=units.mass)
@validate_output(units.kilogram * units.meter / units.second)

def calculate_momuntum(mass_: Quantity, velocity_: Quantity) -> Quantity:
    solved = solve(definition, momentum, dict=True)[0][momentum]
    result_expr = solved.subs({
        mass: mass_,
        velocity: velocity_})
    return expr_to_quantity(result_expr, 'momentum')