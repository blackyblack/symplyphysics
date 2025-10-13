from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import units, errors
from symplyphysics.laws.electricity.vector import electric_field_vector_due_to_dipole as law

from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector, CARTESIAN
from symplyphysics.core.experimental.approx import assert_equal_vectors

Args = namedtuple("Args", "p r e")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    p = QuantityCoordinateVector([
        2e-13 * units.coulomb * units.meter,
        -2e-13 * units.coulomb * units.meter,
        1e-13 * units.coulomb * units.meter,
    ], CARTESIAN)

    r = QuantityCoordinateVector([
        1 * units.meter,
        -1 * units.meter,
        2 * units.meter,
    ], CARTESIAN)

    e = QuantityCoordinateVector([
        6.12e-5 * units.volt / units.meter,
        -6.12e-5 * units.volt / units.meter,
        30.6e-5 * units.volt / units.meter,
    ], CARTESIAN)

    return Args(p=p, r=r, e=e)


def test_electric_field_law(test_args: Args) -> None:
    result = law.calculate_electric_field(test_args.p, test_args.r)
    assert_equal_vectors(result, test_args.e)


def test_electric_dipole_moment(test_args: Args) -> None:
    result = law.calculate_electric_dipole_moment(test_args.e, test_args.r)

    assert_equal_vectors(result, test_args.p)


def test_bad_electric_dipole_moment(test_args: Args) -> None:
    pb_vector = QuantityCoordinateVector([units.candela, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_electric_field(pb_vector, test_args.r)

    pb_scalar = units.coulomb * units.meter
    with raises(ValueError):
        law.calculate_electric_field(pb_scalar, test_args.r)

    with raises(TypeError):
        law.calculate_electric_field(100, test_args.r)
    with raises(TypeError):
        law.calculate_electric_field([100], test_args.r)


def test_bad_position_vector(test_args: Args) -> None:
    rb_vector = QuantityCoordinateVector([units.candela, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_electric_field(test_args.p, rb_vector)
    with raises(errors.UnitsError):
        law.calculate_electric_dipole_moment(test_args.e, rb_vector)

    rb_scalar = units.meter
    with raises(ValueError):
        law.calculate_electric_field(test_args.p, rb_scalar)
    with raises(ValueError):
        law.calculate_electric_dipole_moment(test_args.e, rb_scalar)

    with raises(TypeError):
        law.calculate_electric_field(test_args.p, 100)
    with raises(TypeError):
        law.calculate_electric_field(test_args.p, [100])
    with raises(TypeError):
        law.calculate_electric_dipole_moment(test_args.e, 100)
    with raises(TypeError):
        law.calculate_electric_dipole_moment(test_args.e, [100])


def test_bad_electric_field(test_args: Args) -> None:
    eb_vector = QuantityCoordinateVector([units.candela, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_electric_dipole_moment(eb_vector, test_args.r)

    eb_scalar = units.volt / units.meter
    with raises(ValueError):
        law.calculate_electric_dipole_moment(eb_scalar, test_args.r)

    with raises(TypeError):
        law.calculate_electric_dipole_moment(100, test_args.r)
    with raises(TypeError):
        law.calculate_electric_dipole_moment([100], test_args.r)
