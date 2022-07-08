from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Momentum is the multiplication of velocity and mass
## Definition: P = m * V
## Where:
## m is a mass of the object
## V is it's velocity
## P is momentum

momentum, mass, velocity = symbols('momentum mass velocity')
definition = Eq(momentum, mass * velocity)

definition_dimension_SI = units.kilogram * units.meter / units.second

def print():
    return pretty(definition, use_unicode=False)

def print_dimension():
    return pretty(definition_dimension_SI, use_unicode=False)

@validate_input(velocity_=units.velocity)
@validate_input(mass_=units.mass)
@validate_output(units.momentum)
def calculate_momentum(mass_: Quantity, velocity_: Quantity) -> Quantity:
    solved = solve(definition, momentum, dict=True)[0][momentum]
    result_expr = solved.subs({
        mass: mass_,
        velocity: velocity_})
    return expr_to_quantity(result_expr, 'momentum')

@validate_input(velocity_ =units.velocity)
@validate_input(momentum_ =units.momentum)
@validate_output(units.mass)
def calculate_mass(momentum_: Quantity, velocity_: Quantity) -> Quantity:
    solved = solve(definition, mass, dict=True)[0][mass]
    result_expr = solved.subs({
        momentum: momentum_,
        velocity: velocity_})
    return expr_to_quantity(result_expr, 'mass')

@validate_input(momentum_ = units.momentum)
@validate_input(mass_ = units.mass)
@validate_output(units.velocity)
def calculate_velocity(momentum_: Quantity, mass_: Quantity) -> Quantity:
    solved = solve(definition, velocity, dict=True)[0][velocity]
    result_expr = solved.subs({
        momentum: momentum_,
        mass: mass_})
    return expr_to_quantity(result_expr, 'velocity')