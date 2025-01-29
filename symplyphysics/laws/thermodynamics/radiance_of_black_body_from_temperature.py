"""
Radiance of black body from temperature
=======================================

The Stefan—Boltzmann law, also known as Stefan's law, states that the total energy radiated per
unit surface area of a black body per unit time is directly proportional to the fourth
power of the black body's thermodynamic temperature.

**Notation:**

#. :quantity_notation:`stefan_boltzmann_constant`.

**Conditions:**

#. The body is completely black, i.e. it absorbs all energy.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Stefan%E2%80%93Boltzmann_law>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
    quantities,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import (
    radiant_exitance_is_radiant_flux_emitted_per_unit_area as exitance_def,)
from symplyphysics.laws.thermodynamics import radiation_power_via_temperature as radiation_law

radiance = symbols.radiant_exitance
"""
:symbols:`radiant_exitance` of the body.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the body.
"""

law = Eq(radiance, quantities.stefan_boltzmann_constant * temperature**4)
"""
:laws:symbol::

:laws:latex::
"""

# Derive from law of thermal radiation power

_thermal_radiation_power = radiation_law.law.rhs.subs({
    radiation_law.emissivity: 1,  # see note, epsilon = 1 for idealized black body
    radiation_law.temperature: temperature,
    radiation_law.surface_area: exitance_def.area,
})

_radiant_exitance_derived = exitance_def.definition.rhs.subs(
    exitance_def.radiant_flux(exitance_def.area),
    _thermal_radiation_power,
).doit()

assert expr_equals(_radiant_exitance_derived, law.rhs)


@validate_input(temperature_=temperature)
@validate_output(radiance)
def calculate_radiance(temperature_: Quantity) -> Quantity:
    solved = solve(law, radiance, dict=True)[0][radiance]
    result_expr = solved.subs(temperature, temperature_)
    return Quantity(result_expr)
