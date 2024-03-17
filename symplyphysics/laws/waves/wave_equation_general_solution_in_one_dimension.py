from sympy import Eq, symbols, Function as SymFunction, cos
from symplyphysics import (
    units,
    angle_type,
    Symbol,
    Function,
    Quantity,
    print_expression,
    validate_input,
)
from symplyphysics.core.quantity_decorator import validate_output_same
from symplyphysics.core.symbols.quantities import scale_factor

# Description
## The solution of the 1D wave equation ...

# Law: u(x, t) = f(psi(x, t))
## u - solution of the 1D wave equation
## f - unary function that accepts angle-type input
## psi - [phase of the wave](./phase_of_traveling_wave.py)
## x - spatial variable
## t - time

general_solution = symbols("general_solution", cls=SymFunction)
solution_function = symbols("solution_function", cls=SymFunction)
wave_phase = Function("wave_phase", angle_type)
position = Symbol("position", units.length)
time = Symbol("time", units.time)

law = Eq(general_solution(position, time), solution_function(wave_phase(position, time)))

# TODO: prove it's a solution of the wave equation


def print_law() -> str:
    return print_expression(law)


@validate_input(wave_phase_=wave_phase)
@validate_output_same("amplitude_")
def calculate_displacement(
    amplitude_: Quantity,
    wave_phase_: Quantity | float,
) -> Quantity:
    """Solution of the wave equation in the form `A * cos(psi)`"""

    wave_phase_ = scale_factor(wave_phase_)
    result = law.rhs.subs(
        solution_function(wave_phase(position, time)),
        amplitude_ * cos(wave_phase_),
    )
    return Quantity(result)
