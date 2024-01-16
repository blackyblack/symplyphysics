from symplyphysics import Symbol, units, validate_input, validate_output, Quantity, scale_vector, list_of_quantities
from symplyphysics.core.vectors.vectors import QuantityVector, Vector

# Description
## The electic field E is defined at any point in terms of the electrostatic force F 
## that would be exerted on a test charge q0 placed there.

# Definition: E = F / q0
## E - electic field
## F - electrostatic force
## q0 - test charge

test_charge = Symbol("test_charge", units.charge)


def electric_field_law(electrostatic_force_: Vector) -> Vector:
    return scale_vector(1 / test_charge, electrostatic_force_)


def electrostatic_force_law(electric_field_: Vector) -> Vector:
    return scale_vector(test_charge, electric_field_)


@validate_input(electrostatic_force_=units.force, test_charge_=test_charge)
@validate_output(units.force / units.charge)
def calculate_electric_field(electrostatic_force_: QuantityVector, test_charge_: Quantity) -> QuantityVector:
    result_electric_field = electric_field_law(electrostatic_force_)
    electric_field_components = list_of_quantities(
        result_electric_field.components, 
        {test_charge: test_charge_}
    )
    return QuantityVector(electric_field_components, electrostatic_force_.coordinate_system)


@validate_input(electric_field_=units.force / units.charge, test_charge_=test_charge)
@validate_output(units.force)
def calculate_electrostatic_force(electric_field_: QuantityVector, test_charge_: Quantity) -> QuantityVector:
    result_force = electrostatic_force_law(electric_field_)
    force_components = list_of_quantities(
        result_force.components,
        {test_charge: test_charge_}
    )
    return QuantityVector(force_components, electric_field_.coordinate_system)
