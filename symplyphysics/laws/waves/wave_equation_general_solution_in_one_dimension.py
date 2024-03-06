from sympy import Eq, symbols, Function as SymFunction, cos
from symplyphysics import (
    units,
    angle_type,
    Symbol,
    Quantity,
    print_expression,
    validate_input,
)
from symplyphysics.core.quantity_decorator import validate_output_same
from symplyphysics.core.symbols.quantities import scale_factor

# Description
## The solution of the 1D wave equation are sums of a right-traveling function and a left-traveling function.
## Here, "traveling" means that the shape of these individual arbitrary functions with respect to x stays constant, 
## however, the functions are translated left and right with time with a certain velocity called phase velocity.

# Law: u(x, t) = u_minus(k*x - w*t) + u_plus(k*x + w*t)
## u - solution of the 1D wave equation
## u_minus - solution moving in the positive (right) direction of x-axis
## u_plus - solution moving in the negative (left) direction of x-axis
## k - angular wavenumber
## x - spatial variable
## w - angular frequency
## t - time

general_solution = symbols("general_solution", cls=SymFunction)
left_traveling_solution = symbols("left_traveling_solution", cls=SymFunction)
right_traveling_solution = symbols("right_traveling_solution", cls=SymFunction)
angular_wavenumber = Symbol("angular_wavenumber", angle_type / units.length)
position = Symbol("position", units.length)
angular_frequency = Symbol("angular_frequency", angle_type / units.time)
time = Symbol("time", units.time)

law = Eq(
    general_solution(time),
    right_traveling_solution(angular_wavenumber * position - angular_frequency * time)
    + left_traveling_solution(angular_wavenumber * position + angular_frequency * time)
)

# TODO: prove it's a solution of the wave equation


def print_law() -> str:
    return print_expression(law)


@validate_input(
    angular_wavenumber_=angular_wavenumber,
    position_=position,
    angular_frequency_=angular_frequency,
    time_=time,
    phase_lag_=angle_type,
)
@validate_output_same("amplitude_")
def calculate_displacement(
    amplitude_: Quantity,
    angular_wavenumber_: Quantity,
    position_: Quantity,
    angular_frequency_: Quantity,
    time_: Quantity,
    phase_lag_: Quantity | float,
) -> Quantity:
    """Solution of the wave equation in the form `A * cos(k*x + w*t + phi)`"""

    solution = law.rhs.subs({
        right_traveling_solution(angular_wavenumber * position - angular_frequency * time): (
            amplitude_ * cos(angular_wavenumber * position - angular_frequency * time + scale_factor(phase_lag_))
        ),
        left_traveling_solution(angular_wavenumber * position + angular_frequency * time): 0,
    })
    result = solution.subs({
        angular_wavenumber: angular_wavenumber_,
        position: position_,
        angular_frequency: angular_frequency_,
        time: time_,
    })
    return Quantity(result)
