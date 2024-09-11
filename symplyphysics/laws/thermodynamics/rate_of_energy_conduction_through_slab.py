"""
Rate of energy conduction through slab
======================================

The rate at which energy is conducted through a slab for which one face is maintained at a higher
temperature than the other face is proportional to the temperature difference of the faces and
its face area and inversely proportional to its length.
"""

from sympy import Eq
from symplyphysics import (
    clone_symbol,
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

energy_conduction_rate = Symbol("energy_conduction_rate", units.power, positive=True)
"""
Rate of energy conductivity through the slab.

Symbol:
    :code:`P`
"""

thermal_conductivity = Symbol("thermal_conductivity",
    units.power / (units.length * units.temperature),
    positive=True)
"""
Thermal conductivity of the slab's material.

Symbol:
    :code:`k`
"""

face_area = Symbol("face_area", units.area, positive=True)
"""
Area of the face of the slab.

Symbol:
    :code:`A`
"""

slab_thickness = Symbol("slab_thickness", units.length, positive=True)
"""
Distance between the two faces of the slab.

Symbol:
    :code:`l`
"""

temperature_difference = clone_symbol(symbols.temperature,
    display_symbol="dT",
    display_latex="\\Delta T",
    real=True)
"""
:attr:`~symplyphysics.symbols.temperature` difference between the two faces of the slab.
"""

law = Eq(
    energy_conduction_rate,
    thermal_conductivity * face_area * abs(temperature_difference) / slab_thickness,
)
r"""
:code:`P = k * A * abs(dT) / l`

Latex:
    .. math::
        P = k A \frac{|\Delta T|}{l}
"""


@validate_input(
    thermal_conductivity_=thermal_conductivity,
    face_area_=face_area,
    slab_thickness_=slab_thickness,
    temperature_difference_=temperature_difference,
)
@validate_output(energy_conduction_rate)
def calculate_energy_conduction_rate(
    thermal_conductivity_: Quantity,
    face_area_: Quantity,
    slab_thickness_: Quantity,
    temperature_difference_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        thermal_conductivity: thermal_conductivity_,
        face_area: face_area_,
        slab_thickness: slab_thickness_,
        temperature_difference: temperature_difference_,
    })
    return Quantity(result)
