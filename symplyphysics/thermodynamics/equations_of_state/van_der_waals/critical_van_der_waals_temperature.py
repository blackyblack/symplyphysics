"""
Critical van der Waals temperature
==================================

Critical temperature of a van der Waals fluid depends on the parameters :math:`a` and :math:`b`
of the van der Waals equation and the molar gas constant :math:`R`. See :ref:`Critical parameters <vdw_critical_parameters_def>`.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Van_der_Waals_equation#Critical_point_and_corresponding_states>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
    quantities,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.thermodynamics.equations_of_state.van_der_waals import van_der_vaals_equation as waals_law
from symplyphysics.thermodynamics.phase_transitions.critical_point import critical_point_is_isotherm_stationary_inflection_point as critical_point

critical_temperature = clone_as_symbol(symbols.temperature,
    display_symbol="T_c",
    display_latex="T_\\text{c}")
"""
Critical :symbols:`temperature` of the van der Waals fluid.
"""

attractive_forces_parameter = symbols.attractive_forces_parameter
"""
:symbols:`attractive_forces_parameter`.
"""

excluded_volume_parameter = symbols.excluded_volume_parameter
"""
:symbols:`excluded_volume_parameter`.
"""

law = Eq(
    critical_temperature,
    (8 * attractive_forces_parameter) /
    (27 * quantities.molar_gas_constant * excluded_volume_parameter),
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive the law from the van der Waals equation of state and the conditions on the critical
# point, which is the stationary inflection point of the (p, V) isotherm.

## Express pressure on the isotherm as a function of molar volume.
_isotherm_pressure = solve(waals_law.law, waals_law.pressure)[0].subs(waals_law.molar_volume,
    critical_point.volume)

## At the critical point both the first and the second derivatives of pressure with respect
## to volume vanish.
_inflection_point_eqn = critical_point.inflection_point_condition.subs(
    critical_point.pressure(critical_point.volume), _isotherm_pressure).doit()
_flat_tangent_eqn = critical_point.flat_tangent_condition.subs(
    critical_point.pressure(critical_point.volume), _isotherm_pressure).doit()

## Solve the two conditions for the critical molar volume and the critical temperature.
_critical_point_solved = solve(
    (_inflection_point_eqn, _flat_tangent_eqn),
    (critical_point.volume, waals_law.temperature),
    dict=True,
)[0]

_critical_temperature_derived = _critical_point_solved[waals_law.temperature]

assert expr_equals(_critical_temperature_derived, law.rhs)


@validate_input(
    bonding_forces_parameter_=attractive_forces_parameter,
    molecules_volume_parameter_=excluded_volume_parameter,
)
@validate_output(critical_temperature)
def calculate_critical_temperature(
    bonding_forces_parameter_: Quantity,
    molecules_volume_parameter_: Quantity,
) -> Quantity:
    temperature_ = law.rhs.subs({
        attractive_forces_parameter: bonding_forces_parameter_,
        excluded_volume_parameter: molecules_volume_parameter_,
    })
    return Quantity(temperature_)
