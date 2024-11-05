"""
Temperature derivative via volume derivative
============================================

The **Joule—Thompson effect** describes the change in temperature that accompanies the expansion of
a gas without production of work or transfer of heat, which is in effect an isenthalpic process.

**Notes:**

#. The left-hand side of the equation is also called the **Joule—Thompson coefficient**.

**Conditions:**

#. Particle count is assumed to be constant.
#. Heat capacity is assumed to be independent of temperature.
"""

from sympy import Eq, Derivative, solve
from symplyphysics import (
    Quantity,
    symbols,
    units,
    validate_input,
    validate_output,
    clone_as_function,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.geometry.line import two_point_function, Point2D
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation

pressure = symbols.pressure
"""
:symbols:`pressure` inside the system.
"""

enthalpy = symbols.enthalpy
"""
:symbols:`enthalpy` of the system.
"""

temperature = clone_as_function(symbols.temperature, [pressure, enthalpy])
"""
:symbols:`temperature` of the system as a function of :attr:`~pressure` and
:attr:`~enthalpy`.
"""

volume = clone_as_function(
    symbols.volume,
    [temperature(pressure, enthalpy), pressure],
)
"""
:symbols:`volume` of the system as a function of :attr:`~temperature` and
:attr:`~pressure`.
"""

isobaric_heat_capacity = clone_as_symbol(symbols.heat_capacity, subscript="p")
r"""
:symbols:`heat_capacity` of the system at constant pressure.
"""

law = Eq(Derivative(temperature(pressure, enthalpy), pressure), (1 / isobaric_heat_capacity) *
    (temperature(pressure, enthalpy) * Derivative(volume(temperature(pressure, enthalpy), pressure),
    temperature(pressure, enthalpy)) - volume(temperature(pressure, enthalpy), pressure)))
"""
:laws:symbol::

:laws:latex::
"""

# TODO: derive from enthalpy differential and Maxwell relations.

# Note that the Joule-Thompson effect can only happen in real gases, it is non-existent in the ideal gas.

_volume_expr = solve(ideal_gas_equation.law, ideal_gas_equation.volume)[0].subs({
    ideal_gas_equation.pressure: pressure,
    ideal_gas_equation.temperature: temperature(pressure, enthalpy),
})

_joule_thompson_coefficient = law.rhs.subs(volume(temperature(pressure, enthalpy), pressure),
    _volume_expr).doit()

assert expr_equals(_joule_thompson_coefficient, 0)


@validate_input(
    volume_before_=volume,
    volume_after_=volume,
    temperature_before_=temperature,
    temperature_after_=temperature,
    isobaric_heat_capacity_=isobaric_heat_capacity,
)
@validate_output(units.temperature / units.pressure)
def calculate_temperature_derivative(
    volume_before_: Quantity,
    volume_after_: Quantity,
    temperature_before_: Quantity,
    temperature_after_: Quantity,
    isobaric_heat_capacity_: Quantity,
) -> Quantity:
    temperature_sym = symbols.temperature
    expr = law.rhs.subs(temperature(pressure, enthalpy), temperature_sym)

    volume_ = two_point_function(
        Point2D(temperature_before_, volume_before_),
        Point2D(temperature_after_, volume_after_),
        temperature_sym,
    )

    result = expr.subs(volume(temperature_sym, pressure),
        volume_).doit().subs(isobaric_heat_capacity, isobaric_heat_capacity_).simplify()

    # Result does not depend on temperature
    assert expr_equals(result.diff(temperature_sym), 0)

    return Quantity(result)
