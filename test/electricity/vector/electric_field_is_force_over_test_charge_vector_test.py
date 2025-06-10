from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import units, errors, Quantity
from symplyphysics.laws.electricity.vector import electric_field_is_force_over_test_charge as electric_field

from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector, CARTESIAN
from symplyphysics.core.experimental.approx import assert_equal_vectors

# The electric field at a point is [1, 0, 2] N/C if the force exerted on a test charge of 0.5 C is [0.5, 0, 1] N

Args = namedtuple("Args", "q0 F E")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    q0 = Quantity(0.5 * units.coulomb)
    F = QuantityCoordinateVector([
        Quantity(0.5 * units.newton),
        Quantity(0 * units.newton),
        Quantity(1 * units.newton),
    ], CARTESIAN)
    E = QuantityCoordinateVector([
        Quantity(1 * units.newton / units.coulomb),
        Quantity(0 * units.newton / units.coulomb),
        Quantity(2 * units.newton / units.coulomb)
    ], CARTESIAN)
    return Args(q0=q0, F=F, E=E)


def test_basic_electric_field_definition(test_args: Args) -> None:
    result = electric_field.calculate_electric_field(test_args.F, test_args.q0)
    assert_equal_vectors(result, test_args.E)


def test_basic_electrostatic_force_law(test_args: Args) -> None:
    result = electric_field.calculate_electrostatic_force(test_args.E, test_args.q0)
    assert_equal_vectors(result, test_args.F)


def test_bad_force_in_electric_field(test_args: Args) -> None:
    f_bad_vector = QuantityCoordinateVector([units.second, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        electric_field.calculate_electric_field(f_bad_vector, test_args.q0)

    Fb = Quantity(units.newton)
    with raises(ValueError):
        electric_field.calculate_electric_field(Fb, test_args.q0)

    Fb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        electric_field.calculate_electric_field(Fb, test_args.q0)
    with raises(TypeError):
        electric_field.calculate_electric_field(100, test_args.q0)


def test_bad_charge_in_electric_field(test_args: Args) -> None:
    q0b = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        electric_field.calculate_electric_field(test_args.F, q0b)
    with raises(TypeError):
        electric_field.calculate_electric_field(test_args.F, 100)


def test_bad_electric_field_in_electrostatic_force(test_args: Args) -> None:
    e_bad_vector = QuantityCoordinateVector([units.second, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        electric_field.calculate_electrostatic_force(e_bad_vector, test_args.q0)

    Eb = Quantity(units.newton)
    with raises(ValueError):
        electric_field.calculate_electrostatic_force(Eb, test_args.q0)

    Eb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        electric_field.calculate_electrostatic_force(Eb, test_args.q0)
    with raises(TypeError):
        electric_field.calculate_electrostatic_force(100, test_args.q0)


def test_bad_charge_in_electrostatic_force(test_args: Args) -> None:
    q0b = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        electric_field.calculate_electrostatic_force(test_args.E, q0b)
    with raises(TypeError):
        electric_field.calculate_electrostatic_force(test_args.E, 100)
