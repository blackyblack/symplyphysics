from sympy import Eq, Rational, solve
from symplyphysics import (
    dimensionless,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    convert_to_float,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics.equations_of_state.van_der_waals import (
    equation, critical_molar_volume, critical_pressure, critical_temperature,
)

# Description
## The dimensionless form of the van der Waals equation of state features reduced quantities,
## which are simply the usual thermodynamic quantities divided by their value at the critical
## point. One notable property of the dimensionless equation of state is that it contains no
## substance-specific quantities, i.e. all van der Waals fluids will plot on the same reduced
## pressure-volume curve at the same reduced temperature.

# Law: (p* + 3 / (V*)**2) * (V* - 1/3) = 8/3 * T*
## p* = p/p_c - reduced pressure, see [critical pressure](./critical_pressure.py)
## V* = V/V_c - reduced volume, see [critical volume](./critical_volume.py)
## T* = T/T_c - reduced temperature, see [critical temperature](./critical_temperature.py)

reduced_pressure = Symbol("reduced_pressure", dimensionless)
reduced_volume = Symbol("reduced_volume", dimensionless)
reduced_temperature = Symbol("reduced_temperature", dimensionless)

law = Eq(
    (reduced_pressure + 3 / reduced_volume**2) * (reduced_volume - Rational(1, 3)),
    Rational(8, 3) * reduced_temperature,
)

# Derive from van der Waals equation of state

_critical_pressure_expr = critical_pressure.law.rhs.subs({
    critical_pressure.bonding_forces_parameter: equation.bonding_forces_parameter,
    critical_pressure.molecules_volume_parameter: equation.molecules_volume_parameter,
})

_critical_molar_volume_expr = critical_molar_volume.law.rhs.subs({
    critical_molar_volume.molecules_volume_parameter: equation.molecules_volume_parameter,
})

_critical_temperature_expr = critical_temperature.law.rhs.subs({
    critical_temperature.bonding_forces_parameter: equation.bonding_forces_parameter,
    critical_temperature.molecules_volume_parameter: equation.molecules_volume_parameter,
})

_reduced_eqn = equation.law.subs({
    equation.pressure: reduced_pressure * _critical_pressure_expr,
    equation.volume: reduced_volume * _critical_molar_volume_expr * equation.amount_of_substance,
    equation.temperature: reduced_temperature * _critical_temperature_expr,
})

_reduced_temperature_derived = solve(_reduced_eqn, reduced_temperature)[0]

_reduced_temperature_from_law = solve(law, reduced_temperature)[0]

assert expr_equals(_reduced_temperature_derived, _reduced_temperature_from_law)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    reduced_volume_=reduced_volume,
    reduced_temperature_=reduced_temperature,
)
@validate_output(reduced_pressure)
def calculate_reduced_pressure(
    reduced_volume_: float,
    reduced_temperature_: float,
) -> float:
    expr = solve(law, reduced_pressure)[0]
    result = expr.subs({
        reduced_volume: reduced_volume_,
        reduced_temperature: reduced_temperature_,
    })
    return convert_to_float(result)
