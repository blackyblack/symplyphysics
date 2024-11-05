"""
Rate of energy conduction through slab
======================================

The rate at which energy is conducted through a slab for which one face is maintained at a higher
temperature than the other face is proportional to the temperature difference of the faces and
its face area and inversely proportional to its length.
"""

from sympy import Eq
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
)

energy_conduction_rate = symbols.power
"""
Rate of energy conductivity through the slab. See :symbols:`power`.
"""

thermal_conductivity = symbols.thermal_conductivity
"""
:symbols:`thermal_conductivity` of the slab's material.
"""

face_area = symbols.area
"""
:symbols:`area` of the face of the slab.
"""

slab_thickness = symbols.thickness
"""
Distance between the two faces of the slab. See :symbols:`thickness`.
"""

temperature_difference = clone_as_symbol(
    symbols.temperature,
    display_symbol="Delta(T)",
    display_latex="\\Delta T",
    real=True,
)
"""
:symbols:`temperature` difference between the two faces of the slab.
"""

law = Eq(
    energy_conduction_rate,
    thermal_conductivity * face_area * abs(temperature_difference) / slab_thickness,
)
"""
:laws:symbol::

:laws:latex::
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
