from sympy import Eq
from symplyphysics import (
    clone_symbol,
    symbols,
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The rate at which energy is conducted through a slab for which one face is maintained at a higher
## temperature than the other face is proportional to the temperature difference of the faces and
## its face area and inversely proportional to its length.

# Law: P_cond = k * A * |delta_T| / L
## P_cond - rate of energy conductivity through slab
## k - thermal conductivity of material of slab
## A - area of slab face
## L - thickness of slab (distance between the two faces)
## delta_T - difference between the temperatures of the faces
## |x| - the absolute value of x

energy_conduction_rate = Symbol("energy_conduction_rate", units.power, positive=True)
thermal_conductivity = Symbol("thermal_conductivity",
    units.power / (units.length * units.temperature),
    positive=True)
face_area = Symbol("face_area", units.area, positive=True)
slab_thickness = Symbol("slab_thickness", units.length, positive=True)
temperature_difference = clone_symbol(symbols.thermodynamics.temperature,
    "temperature_difference",
    real=True)

law = Eq(
    energy_conduction_rate,
    thermal_conductivity * face_area * abs(temperature_difference) / slab_thickness,
)


def print_law() -> str:
    return print_expression(law)


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
