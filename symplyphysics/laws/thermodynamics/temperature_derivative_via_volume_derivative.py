from sympy import Eq, Derivative, solve
from symplyphysics import (
    Quantity,
    Symbol,
    Function,
    symbols,
    units,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.geometry.line import two_point_function, Point2D
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation

# Description
## The Joule-Thompson effect describes the change in temperature that accompanies the expansion of
## a gas without production of work or transfer of heat, in effect this is an isenthalpic process.

# Law: (dT/dp)_H = (T * (dV/dT)_p - V)/C_p
## T - temperature
## p - pressure
## V - volume
## H - enthalpy
## C_p - isobaric heat capacity

# Notes
## - The left-hand side is also called the Joule-Thompson coefficient.

# Conditions
## - Changes in particle count are not taken into account and it is assumed to be constant.
## - Heat capacity is assumed to be independent of temperature.

temperature = Function("temperature", units.temperature)
pressure = Symbol("pressure", units.pressure)
volume = Function("volume", units.volume)
enthalpy = Symbol("enthalpy", units.energy)
isobaric_heat_capacity = Symbol("isobaric_heat_capacity", units.energy / units.temperature)

law = Eq(
    Derivative(temperature(pressure, enthalpy), pressure),
    ((temperature(pressure, enthalpy)
      * Derivative(volume(temperature(pressure, enthalpy), pressure), temperature(pressure, enthalpy))
      - volume(temperature(pressure, enthalpy), pressure)
     ) / isobaric_heat_capacity)
)

# TODO: derive from enthalpy differential and Maxwell relations.

# Note that the Joule-Thompson effect can only happen in real gases, it is non-existent in the ideal gas.

_volume_expr = solve(ideal_gas_equation.law, ideal_gas_equation.volume)[0].subs({
    ideal_gas_equation.pressure: pressure,
    ideal_gas_equation.temperature: temperature(pressure, enthalpy),
})

_joule_thompson_coefficient = law.rhs.subs(
    volume(temperature(pressure, enthalpy), pressure), _volume_expr
).doit()

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
    temperature_sym = symbols.thermodynamics.temperature
    expr = law.rhs.subs(temperature(pressure, enthalpy), temperature_sym)

    volume_ = two_point_function(
        Point2D(temperature_before_, volume_before_),
        Point2D(temperature_after_, volume_after_),
        temperature_sym,
    )

    result = expr.subs(
        volume(temperature_sym, pressure), volume_
    ).doit().subs(
        isobaric_heat_capacity, isobaric_heat_capacity_
    ).simplify()

    # Result does not depend on temperature
    assert expr_equals(result.diff(temperature_sym), 0)

    return Quantity(result)
