from sympy import Eq, sqrt, Derivative
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Derived by Laplace, the formula for the speed of sound in fluids uses the pressure-density dependence
## and the fact that the oscillations in a sound wave happen so fast, and the thermal conductivity of air
## is so small, that there is no heat transfer in the sound wave, i.e. it is an adiabatic, and therefore
## isentropic, process. This is in contrast with the Newton's formula, who thought that the sound propagation
## is an isothermal process in the assumption that the temperature differences between different parts of
## the sound wave immediately level out, which eventually turned out to be inconsistent with experimental data.

# Law: c = sqrt((dp/d(rho))_S)
## c - speed of sound in fluid
## p - pressure
## rho - density
## (d/d(rho)_S - derivative w.r.t. density at constant entropy (i.e. adiabatic conditions)

speed_of_sound = Symbol("speed_of_sound", units.velocity)
pressure = Function("pressure", units.pressure)
density = Symbol("density", units.mass / units.volume)
entropy = Symbol("entropy", units.energy / units.temperature)

law = Eq(speed_of_sound, sqrt(Derivative(pressure(density, entropy), density)))


def print_law() -> str:
    return print_expression(law)


@validate_input(
    density_change_=density,
    pressure_change_=pressure,
)
@validate_output(speed_of_sound)
def calculate_speed_of_sound(
    density_change_: Quantity,
    pressure_change_: Quantity,
) -> Quantity:
    pressure_ = (pressure_change_ / density_change_) * density
    result = law.rhs.subs(pressure(density, entropy), pressure_).doit()
    return Quantity(result)
