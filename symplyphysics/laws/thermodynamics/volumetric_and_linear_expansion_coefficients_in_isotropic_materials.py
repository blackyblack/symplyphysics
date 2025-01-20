"""
Volumetric and linear expansion coefficients in isotropic materials
===================================================================

Coefficients of thermal expansion describe how the size of an object changes with a change in temperature
at a constant pressure. In isotropic materials, the volumetric coefficient is three times the linear one.

**Conditions:**

#. The material is isotropic.
#. The pressure is constant during the volume change.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Thermal_expansion#Isotropic_materials>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    clone_as_symbol,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import (
    linear_coefficient_of_thermal_expansion as linear_def,
    volumetric_coefficient_of_thermal_expansion as volumetric_def,
)

volumetric_expansion_coefficient = clone_as_symbol(symbols.thermal_expansion_coefficient,
    subscript="V")
"""
Volumetric :symbols:`thermal_expansion_coefficient` of the material. Also see
:doc:`Volumetric expansion coefficient <definitions.volumetric_coefficient_of_thermal_expansion>`.
"""

linear_expansion_coefficient = clone_as_symbol(symbols.thermal_expansion_coefficient, subscript="l")
r"""
Linear :symbols:`thermal_expansion_coefficient` of the material. Also see
:doc:`Linear expansion coefficient <definitions.linear_coefficient_of_thermal_expansion>`.
"""

law = Eq(volumetric_expansion_coefficient, 3 * linear_expansion_coefficient)
"""
:laws:symbol::

:laws:latex::
"""

# Derive from their definitions for a cube-sized object

_volume_change = clone_as_symbol(symbols.volume, display_symbol="dV")
_initial_length = clone_as_symbol(symbols.length, subscript="0")
_length_change = clone_as_symbol(symbols.length, display_symbol="dl")
_temperature_change = clone_as_symbol(symbols.temperature, display_symbol="dT")

_final_length_expr = _initial_length + _length_change

_volume = clone_as_symbol(symbols.volume)
_length = clone_as_symbol(symbols.length)
_volume_eqn = Eq(_volume, _length**3)
_volume_expr = solve(_volume_eqn, _volume)[0]
_initial_volume_expr = _volume_expr.subs(_length, _initial_length)
_final_volume_expr = _volume_expr.subs(_length, _final_length_expr)

_volumetric_expr = volumetric_def.definition.rhs.subs({
    volumetric_def.volume(volumetric_def.temperature).diff(volumetric_def.temperature):
    _volume_change / _temperature_change,
    _volume_change:
    _final_volume_expr - _initial_volume_expr,
    volumetric_def.volume(volumetric_def.temperature):
        _initial_volume_expr,
})

# Note: if dimension changes are not small, the third and fourth terms of the following series expansion
# can be taken into account.
_volumetric_expr_series = _volumetric_expr.series(_length_change, 0, 2).removeO()

_linear_eqn = linear_def.definition.subs({
    linear_def.length(linear_def.temperature, linear_def.pressure).diff(linear_def.temperature):
    _length_change / _temperature_change,
    linear_def.length(linear_def.temperature, linear_def.pressure):
        _initial_length,
    linear_def.linear_expansion_coefficient:
        linear_expansion_coefficient,
})

_volumetric_expr_derived = solve(
    (Eq(volumetric_expansion_coefficient, _volumetric_expr_series), _linear_eqn),
    (volumetric_expansion_coefficient, _length_change),
    dict=True,
)[0][volumetric_expansion_coefficient]

assert expr_equals(_volumetric_expr_derived, law.rhs)


@validate_input(linear_expansion_coefficient_=linear_expansion_coefficient)
@validate_output(volumetric_expansion_coefficient)
def calculate_volumetric_expansion_coefficient(linear_expansion_coefficient_: Quantity) -> Quantity:
    result = law.rhs.subs(linear_expansion_coefficient, linear_expansion_coefficient_)
    return Quantity(result)
