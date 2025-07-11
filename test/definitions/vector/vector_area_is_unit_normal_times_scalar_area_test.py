from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import units, Quantity, errors
from symplyphysics.definitions.vector import vector_area_is_unit_normal_times_scalar_area as law

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

Args = namedtuple("Args", "n a")


@fixture(name="test_args")
def test_args_fixture() -> None:
    n = QuantityCoordinateVector([0.1, -0.4, 0.911], CARTESIAN)

    a = Quantity(3e-7 * units.meter**2)

    return Args(n=n, a=a)


def test_law(test_args: Args) -> None:
    result = law.calculate_vector_area(test_args.n, test_args.a)

    expected = QuantityCoordinateVector([
        3.0e-8 * units.meter**2,
        -1.2e-7 * units.meter**2,
        2.73e-7 * units.meter**2,
    ], CARTESIAN)

    assert_equal_vectors(result, expected, relative_tolerance=2e-3)


def test_bad_normal(test_args: Args) -> None:
    bad_vector = QuantityCoordinateVector([units.candela, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_vector_area(bad_vector, test_args.a)

    bad_scalar = 1.0
    with raises(ValueError):
        law.calculate_vector_area(bad_scalar, test_args.a)
    with raises(ValueError):
        law.calculate_vector_area([bad_scalar], test_args.a)


def test_bad_area(test_args: Args) -> None:
    bad_scalar = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_vector_area(test_args.n, bad_scalar)
    with raises(TypeError):
        law.calculate_vector_area(test_args.n, 100)
    
