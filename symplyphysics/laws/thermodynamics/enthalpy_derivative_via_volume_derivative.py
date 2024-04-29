from sympy import Eq, Derivative
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.geometry.line import two_point_function, Point2D

# Description
## Isothermal derivative of enthalpy w.r.t. pressure can be found as a function of volume and temperature.

# (dH/dp)_T = V(T, p) - T * (dV/dT)_p
## H - enthalpy
## T - absolute temperature
## p - pressure
## V - volume
## (d/dp)_T - partial derivative w.r.t. pressure at constant temperature
## (d/dT)_p - partial derivative w.r.t. temperature at constant pressure

# Conditions
## - Assuming an infinitesimal quasi-static isothermal process.

enthalpy = Function("enthalpy", units.energy)
temperature = symbols.thermodynamics.temperature
pressure = Symbol("pressure", units.pressure)
volume = Function("volume", units.volume)

law = Eq(
    Derivative(enthalpy(temperature, pressure), pressure),
    volume(temperature, pressure) - temperature * Derivative(volume(temperature, pressure), temperature)
)

# TODO: Derive from differential of enthalpy and Maxwell relation


@validate_input(
    volume_before_=volume,
    volume_after_=volume,
    temperature_before_=temperature,
    temperature_after_=temperature,
)
@validate_output(units.energy / units.pressure)
def calculate_enthalpy_derivative(
    volume_before_: Quantity,
    volume_after_: Quantity,
    temperature_before_: Quantity,
    temperature_after_: Quantity,
) -> Quantity:
    volume_ = two_point_function(
        Point2D(temperature_before_, volume_before_),
        Point2D(temperature_after_, volume_after_),
        temperature,
    )
    result = law.rhs.subs(volume(temperature, pressure), volume_).doit()

    # the result does not depend on temperature
    assert expr_equals(result.diff(temperature), 0)

    return Quantity(result)
