"""
Dimensionless equation
======================

The *dimensionless form* of the van der Waals equation of state features :ref:`reduced quantities
<vdw_reduced_units_def>`. One notable property of the dimensionless equation of state is that it
contains no substance-specific quantities, i.e. all van der Waals fluids will plot on the same
reduced pressure-volume curve at the same reduced temperature.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Van_der_Waals_equation#Critical_point_and_corresponding_states>`__.
"""

from sympy import Eq, Rational, solve
from symplyphysics import (
    dimensionless,
    Symbol,
    validate_input,
    validate_output,
    convert_to_float,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics.equations_of_state.van_der_waals import (
    critical_molar_volume,
    critical_pressure,
    critical_temperature,
    reduced_pressure as reduced_pressure_law,
    reduced_volume as reduced_volume_law,
    reduced_temperature as reduced_temperature_law,
    van_der_vaals_equation,
)

reduced_pressure = Symbol("p_r", dimensionless)
"""
See :doc:`laws.thermodynamics.equations_of_state.van_der_waals.reduced_pressure`.
"""

reduced_volume = Symbol("V_r", dimensionless)
"""
See :doc:`laws.thermodynamics.equations_of_state.van_der_waals.reduced_volume`.
"""

reduced_temperature = Symbol("T_r", dimensionless)
"""
See :doc:`laws.thermodynamics.equations_of_state.van_der_waals.reduced_temperature`.
"""

law = Eq(
    (reduced_pressure + 3 / reduced_volume**2) * (reduced_volume - Rational(1, 3)),
    Rational(8, 3) * reduced_temperature,
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive from van der Waals equation of state

_critical_pressure_expr = critical_pressure.law.rhs.subs({
    critical_pressure.attractive_forces_parameter: van_der_vaals_equation.attractive_forces_parameter,
    critical_pressure.excluded_volume_parameter: van_der_vaals_equation.excluded_volume_parameter,
})

_critical_molar_volume_expr = critical_molar_volume.law.rhs.subs({
    critical_molar_volume.excluded_volume_parameter: van_der_vaals_equation.excluded_volume_parameter,
})

_critical_temperature_expr = critical_temperature.law.rhs.subs({
    critical_temperature.attractive_forces_parameter: van_der_vaals_equation.attractive_forces_parameter,
    critical_temperature.excluded_volume_parameter: van_der_vaals_equation.excluded_volume_parameter,
})

_pressure = solve(reduced_pressure_law.law, reduced_pressure_law.pressure)[0].subs({
    reduced_pressure_law.reduced_pressure: reduced_pressure,
    reduced_pressure_law.critical_pressure: _critical_pressure_expr,
})

_molar_volume = solve(reduced_volume_law.law, reduced_volume_law.volume)[0].subs({
    reduced_volume_law.reduced_volume: reduced_volume,
    reduced_volume_law.critical_volume: _critical_molar_volume_expr,
})

_temperature = solve(reduced_temperature_law.law, reduced_temperature_law.temperature)[0].subs({
    reduced_temperature_law.reduced_temperature: reduced_temperature,
    reduced_temperature_law.critical_temperature: _critical_temperature_expr,
})

_reduced_eqn = van_der_vaals_equation.law.subs({
    van_der_vaals_equation.pressure: _pressure,
    van_der_vaals_equation.molar_volume: _molar_volume,
    van_der_vaals_equation.temperature: _temperature,
})

_reduced_temperature_derived = solve(_reduced_eqn, reduced_temperature)[0]

_reduced_temperature_from_law = solve(law, reduced_temperature)[0]

assert expr_equals(_reduced_temperature_derived, _reduced_temperature_from_law)


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
