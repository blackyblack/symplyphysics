from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import errors, units, Quantity
from symplyphysics.laws.dynamics.vector import normal_force_via_pressure_and_vector_area as law

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

Args = namedtuple("Args", "p a")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    p = Quantity(5 * units.pascal)
    a = QuantityCoordinateVector([
        1e-5 * units.meter**2,
        0,
        -2e-5 * units.meter**2,
    ], CARTESIAN)
    return Args(p=p, a=a)


def test_law(test_args: Args) -> None:
    result = law.calculate_normal_force(test_args.p, test_args.a)

    expected = QuantityCoordinateVector([
        5e-5 * units.newton,
        0,
        -1e-4 * units.newton,
    ], CARTESIAN)

    assert_equal_vectors(result, expected)


def test_bad_pressure(test_args: Args) -> None:
    bad_scalar = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_normal_force(bad_scalar, test_args.a)
    with raises(TypeError):
        law.calculate_normal_force(100, test_args.a)


def test_bad_area(test_args: Args) -> None:
    bad_vector = QuantityCoordinateVector([units.candela, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_normal_force(test_args.p, bad_vector)

    bad_scalar = units.meter**2
    with raises(ValueError):
        law.calculate_normal_force(test_args.p, bad_scalar)

    with raises(TypeError):
        law.calculate_normal_force(test_args.p, 100)
    with raises(TypeError):
        law.calculate_normal_force(test_args.p, [100])
