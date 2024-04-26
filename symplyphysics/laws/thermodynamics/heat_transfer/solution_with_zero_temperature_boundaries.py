from sympy import Eq, sin, exp, pi, solve, dsolve, symbols, Function as SymFunction
from symplyphysics import (
    dimensionless,
    units,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics.heat_transfer import (
    equation_in_homogeneous_medium_in_one_dimension as heat_equation,
)

# Description
## The heat equation coupled with a boundary condition can be solved to get a unique solution.

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
maximum_position = Symbol("maximum_position", units.length, positive=True)
position = Symbol("position", units.length)
time = Symbol("time", units.time)

law = Eq(
    temperature,
    coefficient
    * sin(mode_number * pi * position / maximum_position)
    * exp(-1 * thermal_diffusivity * (mode_number * pi / maximum_position)**2 * time)
)

# Derive law from heat transfer equation. See
# [this site](https://tutorial.math.lamar.edu/classes/de/solvingheatequation.aspx) for reference.

# We are aiming to solve the heat equation with initial boundary conditions `T(0, t) = 0`
# and `T(L, t) = 0`.

# We are using a method of solving partial differential equations called the separation of variables.
# The solution function is separated into two functions, one of the spatial variable and another of
# the time variable.
_position_factor, _time_factor = symbols("position_factor time_factor", cls=SymFunction)

_temperature_expr = _position_factor(position) * _time_factor(time)

_heat_transfer_eqn = heat_equation.law.subs({
    heat_equation.position: position,
    heat_equation.time: time,
    heat_equation.thermal_diffusivity: thermal_diffusivity,
}).subs({
    heat_equation.temperature(position, time): _temperature_expr,
}).doit()

_separation_constant = symbols("separation_constant", positive=True)

_time_factor_via_constant_eqn = Eq(
    _time_factor(time).diff(time) / (_time_factor(time) * thermal_diffusivity),
    -1 * _separation_constant,
)

_time_factor_expr = dsolve(
    _time_factor_via_constant_eqn,
    _time_factor(time),
    ics={_time_factor(0): 1},  # not interested in the coefficient
).rhs

_position_factor_via_constant_expr = solve(
    (_heat_transfer_eqn, _time_factor_via_constant_eqn),
    (_position_factor(position), _time_factor(time)),
    dict=True,
)[0][_position_factor(position)]

_position_factor_via_constant_eqn = Eq(
    _position_factor(position),
    _position_factor_via_constant_expr,
)

_position_factor_expr = dsolve(
    _position_factor_via_constant_eqn,
    _position_factor(position),
    ics={_position_factor(0): 0}  # see initial boundary conditions at left end
).rhs.subs("C1", 1)  # not interested in the coefficient

# Initial boundary condition at right end can only be achieved with certain values of the
# separation constant. Since sympy does not give infinite solutions to trigonometric equations
# this is the correct value of the separation constant for the right-end boundary condition.
_separation_constant_expr = (pi * mode_number / maximum_position)**2
_position_factor_at_right_end = _position_factor_expr.subs(position, maximum_position)
assert expr_equals(_position_factor_at_right_end.subs(_separation_constant, _separation_constant_expr), 0)

_position_factor_expr = _position_factor_expr.subs(_separation_constant, _separation_constant_expr)
_time_factor_expr = _time_factor_expr.subs(_separation_constant, _separation_constant_expr)

# The solution can be multiplied by any coefficient independent of time and position.
_temperature_derived = coefficient * _temperature_expr.subs({
    _position_factor(position): _position_factor_expr,
    _time_factor(time): _time_factor_expr
})

assert expr_equals(_temperature_derived, law.rhs)


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
