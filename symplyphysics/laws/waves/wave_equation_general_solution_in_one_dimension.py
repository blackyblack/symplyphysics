from sympy import Eq, symbols, Function as SymFunction, cos, S
from symplyphysics import (
    units,
    angle_type,
    Symbol,
    Function,
    Quantity,
    print_expression,
    validate_input,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.quantity_decorator import validate_output_same
from symplyphysics.core.symbols.quantities import scale_factor
from symplyphysics.definitions import wave_equation_in_one_dimension as wave_eqn
from symplyphysics.laws.waves import (
    phase_of_traveling_wave as phase_law,
    phase_velocity_from_angular_frequency_and_wavenumber as phase_velocity_law,
)

# Description
## Any function of a single variable can be a solution of the wave equation if
## it depends on the phase of the wave, i.e. the position and time variables
## only appear in its expression in the form of the [phase of the wave](./phase_of_traveling_wave.py)

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

# Prove this is a solution of the wave equation

_angular_frequency, _angular_wavenumber = symbols("angular_frequency angular_wavenumber")

_phase_velocity = phase_velocity_law.law.rhs.subs({
    phase_velocity_law.angular_frequency: _angular_frequency,
    phase_velocity_law.angular_wavenumber: _angular_wavenumber,
})

_phase = phase_law.law.rhs.subs({
    phase_law.angular_wavenumber: _angular_wavenumber,
    phase_law.position: position,
    phase_law.angular_frequency: _angular_frequency,
    phase_law.time: time,
})

_solution = law.rhs.subs(
    wave_phase(position, time), _phase
)

_eqn = wave_eqn.definition.subs({
    wave_eqn.position: position,
    wave_eqn.time: time,
    wave_eqn.phase_velocity: _phase_velocity,
})

_eqn_lhs_subs = _eqn.lhs.subs(
    wave_eqn.displacement(position, time), _solution
).doit()

_eqn_rhs_subs = _eqn.rhs.subs(
    wave_eqn.displacement(position, time), _solution
).doit()

assert expr_equals(_eqn_lhs_subs, _eqn_rhs_subs)


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
