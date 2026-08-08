"""
Voltage is electric field times distance
========================================

Voltage between two points in space can be found as the negative line integral of the electric
field across the path that connects those points. In the case of a constant electric field it can
be simplified to the product of the electric field strength and the distance between the points,
multiplying by the necessary sign.

**Conditions:**

#. The electric field is constant between the two points. This might be achieved by
   choosing a small enough distance between the points.

#. The electric field must be conservative, i.e. this law only applies in the electrostatic case.

**Links:**

#. `Physics LibreTexts <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_II_-_Thermodynamics_Electricity_and_Magnetism_(OpenStax)/07%3A_Electric_Potential/7.03%3A_Electric_Potential_and_Potential_Difference>`__.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols, clone_as_symbol
from symplyphysics.core.expr_comparisons import expr_equals_abs
from symplyphysics.electromagnetism.electrostatics.voltage import voltage_is_line_integral_of_electric_field as voltage_integral_law

voltage = symbols.voltage
"""
:symbols:`voltage` between two points.
"""

electric_field_strength = symbols.electric_field_strength
"""
:symbols:`electric_field_strength`.
"""

distance = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` between two points.
"""

law = Eq(voltage, electric_field_strength * distance)
"""
:laws:symbol::

:laws:latex::
"""

# Derive the law from the line integral of the electric field.

## Assumption: the electric field is constant between the two points, so the tangent component
## of the electric field along the integration path is the constant field strength.
_integrand = voltage_integral_law.law.rhs.subs(
    voltage_integral_law.electric_field_component(voltage_integral_law.distance),
    electric_field_strength,
)

## The integration path is a straight line connecting the two points, so the final distance
## exceeds the initial one by the distance between the points.
_initial_distance = clone_as_symbol(symbols.distance, subscript="0")
_voltage_derived = _integrand.subs({
    voltage_integral_law.initial_distance: _initial_distance,
    voltage_integral_law.final_distance: _initial_distance + distance,
}).doit()

## The sign of the voltage depends on the mutual orientation of the electric field and the
## integration path, hence only the magnitudes are compared.
assert expr_equals_abs(_voltage_derived, law.rhs)


@validate_input(
    electric_field_strength_=electric_field_strength,
    distance_=distance,
)
@validate_output(voltage)
def calculate_voltage(
    electric_field_strength_: Quantity,
    distance_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        electric_field_strength: electric_field_strength_,
        distance: distance_,
    })
    return Quantity(result)
