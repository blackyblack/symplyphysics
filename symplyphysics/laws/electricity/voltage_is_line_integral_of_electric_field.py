"""
Voltage is line integral of electric field
==========================================

In electrostatics, voltage change between two points can be found as the line integral
of the electric field along any path connecting the two points.

**Conditions:**

#. Applies to electrostatic fields, although electrostatic approximation might be
   used for electromagnetic fields at lower frequencies.
"""

from sympy import Eq, Integral
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    clone_as_function,
    clone_as_symbol,
    symbols,
)
from symplyphysics.core.geometry.line import two_point_function, Point2D

voltage = SymbolNew("V", units.voltage)
"""
Voltage between two points.
"""

distance = SymbolNew("s", units.length)
"""
Distance traveled.
"""

electric_field_component = clone_as_function(
    symbols.electric_field_strength,
    [distance],
    subscript="s",
)
"""
Component of the electric field vector tangent to the integration path.
See :symbols:`electric_field_strength`.
"""

initial_distance = SymbolNew("s_0", units.length)
"""
Initial distance.
"""

final_distance = SymbolNew("s_1", units.length)
"""
Final distance.
"""

law = Eq(
    voltage,
    -1 * Integral(electric_field_component(distance), (distance, initial_distance, final_distance)),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    initial_electric_field_component_=electric_field_component,
    final_electric_field_component=electric_field_component,
    initial_distance_=initial_distance,
    final_distance_=final_distance,
)
@validate_output(voltage)
def calculate_voltage(
    initial_electric_field_component_: Quantity,
    final_electric_field_component: Quantity,
    initial_distance_: Quantity,
    final_distance_: Quantity,
) -> Quantity:
    electric_field_component_ = two_point_function(
        Point2D(initial_distance_, initial_electric_field_component_),
        Point2D(final_distance_, final_electric_field_component),
        distance,
    )
    _voltage = law.rhs.subs({
        electric_field_component(distance): electric_field_component_,
        initial_distance: initial_distance_,
        final_distance: final_distance_,
    }).doit()
    return Quantity(_voltage)
