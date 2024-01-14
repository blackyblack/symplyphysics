from sympy import Eq, solve, sqrt

from symplyphysics import (Quantity, Symbol, dimensionless, print_expression,
                           units, validate_input, validate_output)

# Speed of souns for ideal gases
# c = sqrt( gamma * R * T / M ), where
# gamma is heat capacity ratio (adiabatic index),
# T is temperature,
# M is mass of one mole of this gas.
# R is ideal gas constant,
# c is speed of sound.

temperature = Symbol("temperature", units.temperature)
heat_capacity_ratio = Symbol("heat_capacity_ratio", dimensionless)
mole_mass = Symbol("mole_mass", units.mass / units.amount_of_substance)
speed_of_sound = Symbol("speed_of_sound", units.velocity)


law = Eq(
    speed_of_sound,
    sqrt(heat_capacity_ratio * units.molar_gas_constant * temperature / mole_mass)
)


def print_law():
    print_expression(law)


@validate_input(
    temperature_=temperature,
    heat_capacity_ratio_=heat_capacity_ratio,
    mole_mass_=mole_mass
)
@validate_output(speed_of_sound)
def calculate_speed_of_sound(temperature_, heat_capacity_ratio_, mole_mass_):

    result_expr = solve(law, speed_of_sound, dict=True)[0][speed_of_sound]

    result_applied = result_expr.subs(
        {
            temperature: temperature_,
            heat_capacity_ratio: heat_capacity_ratio_,
            mole_mass: mole_mass_,
        }
    )

    return Quantity(result_applied)
