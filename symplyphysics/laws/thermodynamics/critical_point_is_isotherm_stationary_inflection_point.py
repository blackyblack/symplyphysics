from sympy import Eq, Derivative, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    symbols,
)
from symplyphysics.laws.thermodynamics.equations_of_state import (
    van_der_waals_state_equation as van_der_waals_law,
)

# Description
## Critical point (in the thermodynamic sense) is such values of volume, pressure, and temperature at 
## which only one phase exists and at the vicinity of which the physical property of the phases of the 
## substance change dramatically. Algebraically, the critical point is the stationary inflection point
## of the isothermal pressure-volume dependency line.

# Law: (dp/dV)_T = 0, (d**2(p)/dV**2)_T = 0
## p - pressure
## V - volume
## T - temperature
## (d/dV)_T - partial derivative with respect to volume at constant temperature

# Note
## - These equations need to be solved together with the equation of state.

pressure = Function("pressure", units.pressure)
volume = Symbol("volume", units.volume)
temperature = symbols.thermodynamics.temperature

inflection_point_eqn = Eq(Derivative(pressure(volume), volume, 1), 0)
flat_tangent_eqn = Eq(Derivative(pressure(volume), volume, 2), 0)

law = inflection_point_eqn, flat_tangent_eqn

# Find critical point for van der Waals equation of state

_equation_of_state = van_der_waals_law.law.subs({
    van_der_waals_law.volume: volume,
    van_der_waals_law.pressure: pressure(volume),
    van_der_waals_law.temperature: temperature,
})

_pressure_expr = solve(_equation_of_state, pressure(volume))[0]

_inflection_eqn = inflection_point_eqn.subs(
    pressure(volume), _pressure_expr
).doit().simplify()

_tangent_eqn = flat_tangent_eqn.subs(
    pressure(volume), _pressure_expr
).doit().simplify()

_critical_point = solve(
    eqns :=(_equation_of_state, _inflection_eqn, _tangent_eqn),
    (volume, pressure(volume), temperature),
    dict=True,
)[0]

_critical_volume = _critical_point[volume]
_critical_pressure = _critical_point[pressure(volume)]
_critical_temperature = _critical_point[temperature]


def print_law() -> str:
    return print_expression(law)


@validate_input(
    bonding_forces_parameter_=van_der_waals_law.bonding_forces_parameter,
    molecules_volume_parameter_=van_der_waals_law.molecules_volume_parameter,
    amount_of_substance_=van_der_waals_law.amount_of_substance,
)
# FIXME @validate_multiple_output(volume, pressure, temperature)
def calculate_critical_point(
    bonding_forces_parameter_: Quantity,
    molecules_volume_parameter_: Quantity,
    amount_of_substance_: Quantity,
) -> tuple[Quantity, Quantity, Quantity]:
    # Calculate critical point for van der Waals equation of state

    subs_ = {
        van_der_waals_law.bonding_forces_parameter: bonding_forces_parameter_,
        van_der_waals_law.molecules_volume_parameter: molecules_volume_parameter_,
        van_der_waals_law.amount_of_substance: amount_of_substance_,
    }
    volume_ = _critical_volume.subs(subs_)
    pressure_ = _critical_pressure.subs(subs_)
    temperature_ = _critical_temperature.subs(subs_)
    return volume_, pressure_, temperature_
