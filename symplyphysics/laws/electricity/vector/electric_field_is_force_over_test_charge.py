from symplyphysics import Symbol, units, validate_input, validate_output, Quantity, scale_vector
from symplyphysics.core.vectors.vectors import QuantityVector, Vector

# Description
## The electric field E is defined at any point in terms of the electrostatic force F
## that would be exerted on a test charge q0 placed there.

# Definition: E = F / q0
## E - electric field vector
## F - electrostatic force vector
## q0 - test charge

test_charge = Symbol("test_charge", units.charge)


def electric_field_law(electrostatic_force_: Vector) -> Vector:
    return scale_vector(1 / test_charge, electrostatic_force_)


def electrostatic_force_law(electric_field_: Vector) -> Vector:
    return scale_vector(test_charge, electric_field_)


@validate_input(electrostatic_force_=units.force, test_charge_=test_charge)
@validate_output(units.force / units.charge)
def calculate_electric_field(electrostatic_force_: QuantityVector,
    test_charge_: Quantity) -> QuantityVector:
    result_vector = electric_field_law(electrostatic_force_.to_base_vector())
    return QuantityVector.from_base_vector(result_vector, subs={test_charge: test_charge_})


@validate_input(electric_field_=units.force / units.charge, test_charge_=test_charge)
@validate_output(units.force)
def calculate_electrostatic_force(electric_field_: QuantityVector,
    test_charge_: Quantity) -> QuantityVector:
    result_vector = electrostatic_force_law(electric_field_.to_base_vector())
    return QuantityVector.from_base_vector(result_vector, subs={test_charge: test_charge_})
