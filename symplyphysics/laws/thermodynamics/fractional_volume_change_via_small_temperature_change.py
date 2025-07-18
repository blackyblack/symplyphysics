"""
Fractional volume change via small temperature change
=====================================================

Thermal expansion is the phenomenon when a body increases its dimensions in response to an increase
in temperature. For small temperature changes, the expansion coefficient is approximately constant
and the fractional change in the body's volume is proportional to the change in the body's temperature.

**Conditions:**

#. The temperature change :math:`\\Delta T` is small enough for the expansion coefficient
   :math:`\\alpha_V` to be constant.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Thermal_expansion#Volume>`__.

..
    TODO move to `dynamics/deformations`?
"""

from sympy import Eq, symbols as sym_symbols, solve
from symplyphysics import (
    Quantity,
    convert_to_float,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import (
    volumetric_coefficient_of_thermal_expansion as coefficient_def,)
from symplyphysics.laws.quantities import (
    fractional_change_is_change_over_initial_value as fractional_change_law,)

fractional_volume_change = clone_as_symbol(symbols.fractional_change, subscript="V")
"""
Change in the body's :symbols:`volume` divided by its initial volume. See :symbols:`fractional_change`.
"""

volumetric_expansion_coefficient = clone_as_symbol(symbols.thermal_expansion_coefficient,
    subscript="V")
"""
Volumetric :symbols:`thermal_expansion_coefficient`.
"""

temperature_change = clone_as_symbol(symbols.temperature,
    display_symbol="Delta(T)",
    display_latex="\\Delta T")
"""
Change in body's :symbols:`temperature`.
"""

law = Eq(fractional_volume_change, volumetric_expansion_coefficient * temperature_change)
"""
:laws:symbol::

:laws:latex::
"""

# Derive from definition of the coefficient of thermal expansion.

_volume_change = sym_symbols("volume_change")
_parameters = coefficient_def.parameters
_temperature = coefficient_def.temperature
_volume = coefficient_def.volume(_temperature, _parameters)

# Approximate the derivative with a fraction.
_coefficient_eqn = coefficient_def.definition.subs({
    coefficient_def.volumetric_expansion_coefficient: volumetric_expansion_coefficient,
    _volume.diff(_temperature): _volume_change / temperature_change,
})

_volume_change_eqn = fractional_change_law.law.subs({
    fractional_change_law.fractional_change: fractional_volume_change,
    fractional_change_law.change: _volume_change,
    fractional_change_law.initial_value: _volume,
})

_fractional_change_expr = solve(
    (_coefficient_eqn, _volume_change_eqn),
    (_volume_change, fractional_volume_change),
    dict=True,
)[0][fractional_volume_change]

assert expr_equals(_fractional_change_expr, law.rhs)


@validate_input(
    volumetric_expansion_coefficient_=volumetric_expansion_coefficient,
    temperature_change_=temperature_change,
)
@validate_output(fractional_volume_change)
def calculate_fractional_volume_change(
    volumetric_expansion_coefficient_: Quantity,
    temperature_change_: Quantity,
) -> float:
    result = law.rhs.subs({
        volumetric_expansion_coefficient: volumetric_expansion_coefficient_,
        temperature_change: temperature_change_,
    })
    return convert_to_float(result)
