from sympy import Eq, sin, exp, pi
from symplyphysics import (
    dimensionless,
    units,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
)

# Description
## ...

# Law: T_n(x, t) = B_n * sin(n * pi * x / L) * exp(-chi * (n * pi / L)**2 * t)
## T_n - temperature function, solution to the heat equation
## B_n - real-valued coefficient
## x - position, spatial coordinate
## t - time
## L - length of rod
## chi - thermal diffusivity
## n - any positive integer (1, 2, 3, ...)

# Notes
## - `x` can take a value in the range `[0, L]`.
## - The total solution is `T(x, t) = sum(T_n(x, t), n)`
## - To find the values of B_n one must refer to the boundary condition `f(x) = sum(T_n(x, 0), n)`.

temperature = Symbol("temperature", units.temperature)
coefficient = Symbol("coefficient", units.temperature)
thermal_diffusivity = Symbol("thermal_diffusivity", units.area / units.time)
mode_number = Symbol("mode_number", dimensionless, integer=True, positive=True)
maximum_position = Symbol("maximum_position", units.length)
position = Symbol("position", units.length)
time = Symbol("time", units.time)

law = Eq(
    temperature,
    coefficient
    * sin(mode_number * pi * position / maximum_position)
    * exp(-1 * thermal_diffusivity * (mode_number * pi / maximum_position)**2 * time)
)


@validate_input(
    coefficient_=coefficient,
    thermal_diffusivity_=thermal_diffusivity,
    mode_number_=mode_number,
    maximum_position_=maximum_position,
    position_=position,
    time_=time,
)
@validate_output(temperature)
def calculate_temperature(
    coefficient_: Quantity,
    thermal_diffusivity_: Quantity,
    mode_number_: int,
    maximum_position_: Quantity,
    position_: Quantity,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        coefficient: coefficient_,
        thermal_diffusivity: thermal_diffusivity_,
        mode_number: mode_number_,
        maximum_position: maximum_position_,
        position: position_,
        time: time_,
    })
    return Quantity(result)
