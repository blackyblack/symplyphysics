from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematic import distance_from_constant_velocity as velocity_definition

# Description
## Wavelength is the spatial period of a periodic wave — the distance over which the wave's shape repeats.
## It is the distance between consecutive corresponding points of the same phase on the wave.
## Law: λ = v * T, where
## λ is wavelength,
## v is wave propagation speed (phase speed),
## T is oscillation period.

# Conditions:
## - v is constant.

wavelength = Symbol("wavelength", units.length)
propagation_speed = Symbol("propagation_speed", units.velocity)
oscillation_period = Symbol("oscillation_period", units.time)

law = Eq(wavelength, propagation_speed * oscillation_period)

# Derive the same law from constant velocity motion, assuming wave length is a distance that wave travels
# during oscillation period.

# Prove that derived movement function equals to law.rhs, given initial position = 0
# and propagation_speed is constant_velocity
constant_velocity_movement_definition = velocity_definition.law.subs({
    velocity_definition.constant_velocity: propagation_speed,
    velocity_definition.movement_time: oscillation_period,
    velocity_definition.initial_position: 0
})
assert expr_equals(constant_velocity_movement_definition.rhs, law.rhs)


def print() -> str:
    return print_expression(law)


@validate_input_symbols(velocity_=propagation_speed, period_=oscillation_period)
@validate_output_symbol(wavelength)
def calculate_wavelength(velocity_: Quantity, period_: Quantity) -> Quantity:
    applied_definition = solve(law, wavelength, dict=True)[0][wavelength]
    result_expr = applied_definition.subs({
        propagation_speed: velocity_,
        oscillation_period: period_
    })
    return expr_to_quantity(result_expr)
