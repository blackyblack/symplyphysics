from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    units,
    errors,
    Quantity,
    QuantityVector,
)
from symplyphysics.laws.electricity.vector import electric_field_is_force_over_test_charge as electric_field

# The electric field at a point is [1, 0, 2] N/C if the force exerted on a test charge of 0.5 C is [0.5, 0, 1] N


@fixture(name="test_args")
def test_args_fixture():
    q0 = Quantity(0.5 * units.coulomb)
    F = QuantityVector(
        [Quantity(0.5 * units.newton),
        Quantity(0 * units.newton),
        Quantity(1 * units.newton)])
    E = QuantityVector([
        Quantity(1 * units.newton / units.coulomb),
        Quantity(0 * units.newton / units.coulomb),
        Quantity(2 * units.newton / units.coulomb)
    ])
    Args = namedtuple("Args", "q0 F E")
    return Args(q0=q0, F=F, E=E)


def test_basic_electric_field_definition(test_args):
    result = electric_field.calculate_electric_field(test_args.F, test_args.q0)
    assert len(result.components) == 3
    for result_component, correct_component in zip(result.components, test_args.E.components):
        assert_equal(result_component, correct_component)


def test_basic_electrostatic_force_law(test_args):
    result = electric_field.calculate_electrostatic_force(test_args.E, test_args.q0)
    assert len(result.components) == 3
    for result_component, correct_component in zip(result.components, test_args.F.components):
        assert_equal(result_component, correct_component)


def test_bad_force_in_electric_field(test_args):
    Fb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        electric_field.calculate_electric_field(Fb, test_args.q0)
    with raises(TypeError):
        electric_field.calculate_electric_field(100, test_args.q0)


def test_bad_charge_in_electric_field(test_args):
    q0b = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        electric_field.calculate_electric_field(test_args.F, q0b)
    with raises(TypeError):
        electric_field.calculate_electric_field(test_args.F, 100)


def test_bad_electric_field_in_electrostatic_force(test_args):
    Eb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        electric_field.calculate_electrostatic_force(Eb, test_args.q0)
    with raises(TypeError):
        electric_field.calculate_electrostatic_force(100, test_args.q0)


def test_bad_charge_in_electrostatic_force(test_args):
    q0b = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        electric_field.calculate_electrostatic_force(test_args.E, q0b)
    with raises(TypeError):
        electric_field.calculate_electrostatic_force(test_args.E, 100)
